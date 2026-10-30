# ============================================================
# 2×2 因子对照实验：冗余 × 泄漏
# 目的：分解 SHAP 崩溃中，冗余和泄漏各自的贡献
# 条件：
#   A：无冗余 + 无泄漏（基线）
#   B：有冗余 + 无泄漏（隔离冗余效应）
#   C：无冗余 + 有泄漏（隔离泄漏效应）
#   D：有冗余 + 有泄漏（完整现象）
# ============================================================

# ---------------------------- 0. 依赖 ----------------------------
required_pkgs <- c("tidyverse", "ranger", "fastshap", "pROC", "glmnet")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(tidyverse)
library(ranger)
library(fastshap)
library(pROC)
library(glmnet)

set.seed(2026)

# ---------------------------- 1. 辅助函数 ----------------------------

empirical_wasserstein1 <- function(x, y, n_quantiles = 80) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 5 || length(y) < 5) return(0)
  probs <- seq(0.02, 0.98, length.out = n_quantiles)
  mean(abs(quantile(x, probs, na.rm = TRUE) - quantile(y, probs, na.rm = TRUE)))
}

spearman_cor <- function(a, b) {
  s <- suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
  if (is.na(s)) 0 else s
}

jaccard_set <- function(a, b) {
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) 0 else inter / uni
}

compute_shap_importance <- function(rf, X, nsim = 25, sample_size = 150) {
  pfun <- function(object, newdata) {
    predict(object, data = newdata)$predictions[, 2]
  }
  if (nrow(X) > sample_size) {
    idx <- sample(nrow(X), sample_size)
    X_sub <- X[idx, , drop = FALSE]
  } else {
    X_sub <- X
  }
  shap_mat <- tryCatch({
    explain(rf, X = X_sub, pred_wrapper = pfun, nsim = nsim,
            adjust = TRUE, shap_only = TRUE)
  }, error = function(e) NULL)
  if (is.null(shap_mat)) return(NULL)
  colMeans(abs(shap_mat))
}

# ---------------------------- 2. 数据生成（支持 rho = 0） ----------------------------

simulate_data <- function(n = 600,
                          n_groups = 4,
                          features_per_group = 5,
                          rho_within = 0.85,
                          missing_rate = 0.1,
                          label_noise = 0.1,
                          seed = 123) {
  set.seed(seed)
  p <- n_groups * features_per_group
  
  g_labels <- rep(seq_len(n_groups), each = features_per_group)
  f_labels <- rep(seq_len(features_per_group), times = n_groups)
  feat_names <- paste0("G", g_labels, "_F", f_labels)
  
  if (rho_within > 0.05) {
    # ---------- 有冗余：组内高相关 ----------
    X_list <- list()
    base_list <- list()
    for (g in seq_len(n_groups)) {
      base <- rnorm(n)
      base_list[[g]] <- base
      G <- sapply(seq_len(features_per_group), function(j) {
        rho_within * base + sqrt(1 - rho_within^2) * rnorm(n)
      })
      X_list[[g]] <- G
    }
    X <- do.call(cbind, X_list)
    colnames(X) <- feat_names
    
    # 信号通过潜在因子 → 组内成员可互换
    group_beta <- rnorm(n_groups, mean = 1.3, sd = 0.3) *
      sample(c(-1, 1), n_groups, replace = TRUE)
    logit <- rep(0, n)
    for (g in seq_len(n_groups)) {
      logit <- logit + group_beta[g] * base_list[[g]]
    }
  } else {
    # ---------- 无冗余：独立特征，信号特征强度递减 ----------
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- feat_names
    
    # 每组第 1 个特征携带信号，强度按组递减
    signal_idx <- seq(1, p, by = features_per_group)
    betas <- sort(abs(rnorm(n_groups, 2.0, 0.2)), decreasing = TRUE) *
      sample(c(-1, 1), n_groups, replace = TRUE)
    logit <- rep(0, n)
    for (g in seq_len(n_groups)) {
      logit <- logit + betas[g] * X[, signal_idx[g]]
    }
  }
  
  logit <- logit + rnorm(n, 0, 0.5)
  prob <- plogis(logit)
  prob <- (1 - label_noise) * prob + label_noise * 0.5
  y <- rbinom(n, 1, prob)
  
  if (missing_rate > 0) {
    for (j in seq_len(p)) {
      X[sample(n, floor(n * missing_rate)), j] <- NA
    }
  }
  
  df <- as.data.frame(X)
  df$y <- y
  attr(df, "p") <- p
  df
}

# ---------------------------- 3. LASSO 选择 ----------------------------

lasso_select <- function(X, y, top_k, seed) {
  set.seed(seed)
  X_mat <- as.matrix(X)
  cvfit <- tryCatch(
    cv.glmnet(X_mat, y, family = "binomial", alpha = 1,
              nfolds = 5, standardize = TRUE),
    error = function(e) NULL
  )
  if (is.null(cvfit)) {
    cors <- sapply(X, function(col) abs(cor(col, y, use = "complete.obs")))
    return(names(sort(cors, decreasing = TRUE))[seq_len(min(top_k, ncol(X)))])
  }
  coefs <- as.numeric(coef(cvfit, s = "lambda.min"))[-1]
  names(coefs) <- colnames(X)
  ord <- order(abs(coefs), decreasing = TRUE)
  nz <- names(coefs)[abs(coefs) > 1e-8]
  if (length(nz) >= top_k) {
    return(names(coefs)[ord[seq_len(top_k)]])
  }
  cors <- sapply(X, function(col) abs(cor(col, y, use = "complete.obs")))
  backup <- setdiff(names(sort(cors, decreasing = TRUE)), nz)
  c(nz, backup[seq_len(top_k - length(nz))])
}

# ---------------------------- 4. 两种管道 ----------------------------

run_clean_iteration <- function(train, test, p, top_k, seed) {
  set.seed(seed)
  med <- sapply(train[, seq_len(p)], median, na.rm = TRUE)
  sds <- sapply(train[, seq_len(p)], sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds < 1e-8] <- 1
  
  prep <- function(d) {
    for (j in seq_len(p)) {
      d[[j]][is.na(d[[j]])] <- med[j]
      d[[j]] <- (d[[j]] - med[j]) / sds[j]
    }
    d
  }
  tr <- prep(train); te <- prep(test)
  sel <- lasso_select(tr[, seq_len(p)], tr$y, top_k, seed)
  list(train = tr[, c(sel, "y")],
       test  = te[, c(sel, "y")],
       selected = sel)
}

run_leak_iteration <- function(train, test, p, top_k, leak_strength, seed) {
  set.seed(seed)
  n_test <- nrow(test)
  n_leak <- round(leak_strength * n_test)
  if (n_leak > 0) {
    leak_idx <- sample(n_test, n_leak)
    fit_df <- rbind(train, test[leak_idx, , drop = FALSE])
  } else {
    fit_df <- train
  }
  
  med <- sapply(fit_df[, seq_len(p)], median, na.rm = TRUE)
  sds <- sapply(fit_df[, seq_len(p)], sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds < 1e-8] <- 1
  
  prep <- function(d) {
    for (j in seq_len(p)) {
      d[[j]][is.na(d[[j]])] <- med[j]
      d[[j]] <- (d[[j]] - med[j]) / sds[j]
    }
    d
  }
  tr <- prep(train); te <- prep(test); fp <- prep(fit_df)
  
  sel <- lasso_select(fp[, seq_len(p)], fp$y, top_k, seed)
  list(train = tr[, c(sel, "y")],
       test  = te[, c(sel, "y")],
       selected = sel)
}

# ---------------------------- 5. 管道稳定性汇总 ----------------------------

summarize_pipeline <- function(results, n_iter, shap_nsim, seed_offset = 0) {
  aucs <- numeric(n_iter)
  imps <- vector("list", n_iter)
  sels <- lapply(results, `[[`, "selected")
  
  for (i in seq_len(n_iter)) {
    set.seed(seed_offset + 1000 + i)
    rf <- ranger(y ~ ., data = results[[i]]$train,
                 probability = TRUE, num.trees = 200,
                 importance = "impurity")
    pred <- predict(rf, data = results[[i]]$test)$predictions[, 2]
    aucs[i] <- tryCatch(
      as.numeric(roc(results[[i]]$test$y, pred, quiet = TRUE)$auc),
      error = function(e) NA_real_
    )
    X_train <- results[[i]]$train %>% dplyr::select(-y)
    imps[[i]] <- compute_shap_importance(rf, X_train, nsim = shap_nsim)
  }
  
  all_feats <- unique(unlist(sels))
  aligned <- matrix(0, n_iter, length(all_feats),
                    dimnames = list(NULL, all_feats))
  for (i in seq_len(n_iter)) {
    if (!is.null(imps[[i]])) {
      aligned[i, names(imps[[i]])] <- imps[[i]]
    }
  }
  
  cors <- c()
  for (i in seq_len(n_iter - 1)) {
    for (j in (i + 1):n_iter) {
      cors <- c(cors, spearman_cor(aligned[i, ], aligned[j, ]))
    }
  }
  
  jac <- c()
  for (i in seq_len(n_iter - 1)) {
    for (j in (i + 1):n_iter) {
      jac <- c(jac, jaccard_set(sels[[i]], sels[[j]]))
    }
  }
  
  tibble(
    auc_mean = mean(aucs, na.rm = TRUE),
    auc_sd   = sd(aucs, na.rm = TRUE),
    shap_spearman_mean = mean(cors, na.rm = TRUE),
    shap_spearman_sd   = sd(cors, na.rm = TRUE),
    feature_jaccard_mean = mean(jac, na.rm = TRUE),
    n_unique_subsets = length(unique(sapply(sels, paste, collapse = ",")))
  )
}

# ---------------------------- 6. 单次模拟 ----------------------------

run_simulation <- function(sim_params, n_iter = 10, shap_nsim = 25) {
  data_params <- sim_params[c("n", "n_groups", "features_per_group",
                              "rho_within", "missing_rate",
                              "label_noise", "seed")]
  df <- do.call(simulate_data, data_params)
  p <- attr(df, "p")
  top_k <- sim_params$top_k
  
  set.seed(sim_params$seed + 1)
  n <- nrow(df)
  train_idx <- sample(n, floor(0.7 * n))
  train <- df[train_idx, ]
  test  <- df[-train_idx, ]
  
  clean_res <- vector("list", n_iter)
  leak_res  <- vector("list", n_iter)
  
  for (i in seq_len(n_iter)) {
    s <- sim_params$seed * 100 + i
    clean_res[[i]] <- run_clean_iteration(train, test, p, top_k, s)
    leak_res[[i]]  <- run_leak_iteration(train, test, p, top_k,
                                         sim_params$leak_strength, s)
  }
  
  clean_summary <- summarize_pipeline(clean_res, n_iter, shap_nsim)
  leak_summary  <- summarize_pipeline(leak_res,  n_iter, shap_nsim)
  
  tibble(
    leak_strength     = sim_params$leak_strength,
    rho_within        = sim_params$rho_within,
    n                 = sim_params$n,
    top_k             = top_k,
    seed              = sim_params$seed,
    clean_auc_mean    = clean_summary$auc_mean,
    clean_auc_sd      = clean_summary$auc_sd,
    clean_shap_spearman = clean_summary$shap_spearman_mean,
    clean_jaccard     = clean_summary$feature_jaccard_mean,
    leak_auc_mean     = leak_summary$auc_mean,
    leak_auc_sd       = leak_summary$auc_sd,
    leak_shap_spearman = leak_summary$shap_spearman_mean,
    leak_jaccard      = leak_summary$feature_jaccard_mean,
    delta_auc         = leak_summary$auc_mean - clean_summary$auc_mean,
    shap_drop         = clean_summary$shap_spearman_mean -
      leak_summary$shap_spearman_mean
  )
}

# ---------------------------- 7. 2×2 因子对照 ----------------------------

run_factorial <- function(n_sims = 30, n_iter = 10, shap_nsim = 25,
                          leak_range = c(0.3, 0.7),
                          seed_start = 8000) {
  set.seed(seed_start)
  all_results <- list()
  
  for (i in seq_len(n_sims)) {
    s <- seed_start + i
    n_i <- sample(c(400, 600), 1)
    top_k_i <- sample(c(3, 4), 1)
    leak_i <- runif(1, leak_range[1], leak_range[2])
    
    conditions <- list(
      list(rho = 0,    leak = 0,      label = "A_无冗余_无泄漏"),
      list(rho = 0.85, leak = 0,      label = "B_有冗余_无泄漏"),
      list(rho = 0,    leak = leak_i, label = "C_无冗余_有泄漏"),
      list(rho = 0.85, leak = leak_i, label = "D_有冗余_有泄漏")
    )
    
    for (cond in conditions) {
      sp <- list(n = n_i, n_groups = 4, features_per_group = 5,
                 rho_within = cond$rho,
                 missing_rate = 0.1, label_noise = 0.1,
                 leak_strength = cond$leak,
                 top_k = top_k_i, seed = s)
      res <- tryCatch(
        run_simulation(sp, n_iter = n_iter, shap_nsim = shap_nsim),
        error = function(e) NULL
      )
      if (!is.null(res)) {
        res$condition <- cond$label
        all_results[[length(all_results) + 1]] <- res
      }
    }
    
    if (i %% 5 == 0) cat("完成", i, "/", n_sims, "\n")
  }
  bind_rows(all_results)
}

# ---------------------------- 8. 执行 ----------------------------

cat("\n=== 2×2 对照实验开始 ===\n")
cat("条件数：4，每组", 30, "个模拟\n")

factorial_df <- run_factorial(
  n_sims     = 30,
  n_iter     = 10,
  shap_nsim  = 25,
  leak_range = c(0.3, 0.7),
  seed_start = 8000
)

saveRDS(factorial_df, "factorial_v4.rds")
cat("\n结果已保存至 factorial_v4.rds\n")

# ---------------------------- 9. 汇总 ----------------------------

cat("\n=== 按条件汇总 ===\n")
summary_factorial <- factorial_df %>%
  group_by(condition) %>%
  summarise(
    n                   = n(),
    mean_clean_auc      = mean(clean_auc_mean,      na.rm = TRUE),
    mean_leak_auc       = mean(leak_auc_mean,       na.rm = TRUE),
    mean_abs_delta_auc  = mean(abs(delta_auc),      na.rm = TRUE),
    mean_clean_shap     = mean(clean_shap_spearman, na.rm = TRUE),
    mean_leak_shap      = mean(leak_shap_spearman,  na.rm = TRUE),
    mean_shap_drop      = mean(shap_drop,           na.rm = TRUE),
    sd_shap_drop        = sd(shap_drop,             na.rm = TRUE),
    mean_clean_jaccard  = mean(clean_jaccard,       na.rm = TRUE),
    mean_leak_jaccard   = mean(leak_jaccard,        na.rm = TRUE),
    .groups = "drop"
  )
print(summary_factorial)

# ---------------------------- 10. 关键对比 ----------------------------

cat("\n=== 关键对比 ===\n")
get_val <- function(label, col) {
  summary_factorial %>%
    filter(condition == label) %>%
    pull(!!sym(col))
}

cat("\n[对比 1] 冗余对干净基线的影响\n")
A_shap <- get_val("A_无冗余_无泄漏", "mean_clean_shap")
B_shap <- get_val("B_有冗余_无泄漏", "mean_clean_shap")
cat("  无冗余 + 无泄漏 clean_shap:", round(A_shap, 3), "\n")
cat("  有冗余 + 无泄漏 clean_shap:", round(B_shap, 3), "\n")
cat("  → 纯冗余的贡献:", round(A_shap - B_shap, 3), "\n")

cat("\n[对比 2] 泄漏在无冗余下的影响\n")
C_drop <- get_val("C_无冗余_有泄漏", "mean_shap_drop")
cat("  无冗余 + 有泄漏 shap_drop:", round(C_drop, 3), "\n")

cat("\n[对比 3] 泄漏在有冗余下的影响\n")
D_drop <- get_val("D_有冗余_有泄漏", "mean_shap_drop")
cat("  有冗余 + 有泄漏 shap_drop:", round(D_drop, 3), "\n")

cat("\n[对比 4] 交互项（冗余 × 泄漏）\n")
A_drop <- get_val("A_无冗余_无泄漏", "mean_shap_drop")
B_drop <- get_val("B_有冗余_无泄漏", "mean_shap_drop")
interaction <- (D_drop - B_drop) - (C_drop - A_drop)
cat("  交互项 =", round(interaction, 3), "\n")
cat("  → 越大说明冗余和泄漏协同作用越强\n")

cat("\n[对比 5] AUC 稳定性（所有条件）\n")
print(summary_factorial %>%
        dplyr::select(condition, mean_abs_delta_auc, mean_clean_auc, mean_leak_auc))

# ---------------------------- 11. 可视化 ----------------------------

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  p_bar <- ggplot(summary_factorial,
                  aes(x = condition, y = mean_shap_drop, fill = condition)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = mean_shap_drop - sd_shap_drop / sqrt(n),
                      ymax = mean_shap_drop + sd_shap_drop / sqrt(n)),
                  width = 0.2) +
    labs(title = "SHAP 一致性下降（按条件）",
         x = NULL, y = "Mean SHAP drop") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))
  
  p_box <- ggplot(factorial_df,
                  aes(x = condition, y = shap_drop, fill = condition)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
    labs(title = "SHAP 下降分布（按条件）",
         x = NULL, y = "SHAP drop") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))
  
  p_scatter <- ggplot(factorial_df,
                      aes(x = leak_strength, y = shap_drop,
                          color = condition)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(title = "泄漏强度 vs SHAP 下降",
         x = "Leak Strength", y = "SHAP drop",
         color = "Condition") +
    theme_bw()
  
  print(p_bar)
  print(p_box)
  print(p_scatter)
}

# ---------------------------- 12. 结语 ----------------------------

cat("\n=== 全部完成 ===\n")
cat("结果文件: factorial_v4.rds\n")
cat("\n判读提示:\n")
cat("  若 A≈C 且 B<A 且 D>B → 冗余是前提，泄漏是放大器\n")
cat("  若 C 也明显 > 0 → 泄漏是独立因素\n")
cat("  若所有条件 shap_drop 都接近 0 → 现象可能是噪声\n")