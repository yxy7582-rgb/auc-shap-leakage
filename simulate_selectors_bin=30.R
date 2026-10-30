# ============================================================
# Leakage × SHAP Instability — P1 + P2 + Publication Figures
# 修改版：P1 每 bin n=30 + 正式统计检验（二次项 + 峰值对比）
# ============================================================

# ---------------------------- 工作目录 ----------------------------
cat("Working directory:", getwd(), "\n")

# ---------------------------- 运行开关 ----------------------------
RUN_P1 <- TRUE    # 如果 rds 已存在且不想重跑，可改成 FALSE
RUN_P2 <- TRUE

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
spearman_cor <- function(a, b) {
  s <- suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
  if (is.na(s)) 0 else s
}

jaccard_set <- function(a, b) {
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) 0 else inter / uni
}

compute_shap_importance <- function(rf, X, nsim = 20, sample_size = 150) {
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

# ---------------------------- 2. 数据生成 ----------------------------
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
    
    group_beta <- rnorm(n_groups, mean = 1.3, sd = 0.3) *
      sample(c(-1, 1), n_groups, replace = TRUE)
    logit <- rep(0, n)
    for (g in seq_len(n_groups)) {
      logit <- logit + group_beta[g] * base_list[[g]]
    }
  } else {
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- feat_names
    
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

# ---------------------------- 3. 多选择器支持 ----------------------------
select_features <- function(X, y, top_k, seed, method = "lasso") {
  set.seed(seed)
  X_mat <- as.matrix(X)
  
  if (method == "lasso" || method == "elastic_net") {
    alpha_val <- if (method == "lasso") 1 else 0.5
    cvfit <- tryCatch(
      cv.glmnet(X_mat, y, family = "binomial", alpha = alpha_val,
                nfolds = 5, standardize = TRUE),
      error = function(e) NULL
    )
    if (is.null(cvfit)) {
      cors <- sapply(X, function(col) abs(cor(col, y, use = "complete.obs")))
      return(names(sort(cors, decreasing = TRUE))[seq_len(min(top_k, ncol(X)))])
    }
    coefs <- as.numeric(coef(cvfit, s = "lambda.min"))[-1]
    names(coefs) <- colnames(X)
    nz <- names(coefs)[abs(coefs) > 1e-8]
    if (length(nz) >= top_k) {
      ord <- order(abs(coefs), decreasing = TRUE)
      return(names(coefs)[ord[seq_len(top_k)]])
    }
    cors <- sapply(X, function(col) abs(cor(col, y, use = "complete.obs")))
    backup <- setdiff(names(sort(cors, decreasing = TRUE)), nz)
    return(c(nz, backup[seq_len(top_k - length(nz))]))
    
  } else if (method == "univariate") {
    cors <- sapply(X, function(col) {
      v <- suppressWarnings(abs(cor(col, y, use = "complete.obs")))
      if (is.na(v)) 0 else v
    })
    return(names(sort(cors, decreasing = TRUE))[seq_len(min(top_k, ncol(X)))])
    
  } else {
    stop("Unknown selector: ", method)
  }
}

# ---------------------------- 4. 两种管道 ----------------------------
run_clean_iteration <- function(train, test, p, top_k, seed,
                                selector_method = "lasso") {
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
  sel <- select_features(tr[, seq_len(p)], tr$y, top_k, seed, selector_method)
  list(train = tr[, c(sel, "y")],
       test  = te[, c(sel, "y")],
       selected = sel)
}

run_leak_iteration <- function(train, test, p, top_k, leak_strength, seed,
                               selector_method = "lasso") {
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
  
  sel <- select_features(fp[, seq_len(p)], fp$y, top_k, seed, selector_method)
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
run_simulation <- function(sim_params, n_iter = 10, shap_nsim = 20,
                           selector_method = "lasso") {
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
    clean_res[[i]] <- run_clean_iteration(train, test, p, top_k, s,
                                          selector_method)
    leak_res[[i]]  <- run_leak_iteration(train, test, p, top_k,
                                         sim_params$leak_strength, s,
                                         selector_method)
  }
  
  clean_summary <- summarize_pipeline(clean_res, n_iter, shap_nsim)
  leak_summary  <- summarize_pipeline(leak_res,  n_iter, shap_nsim)
  
  tibble(
    leak_strength     = sim_params$leak_strength,
    rho_within        = sim_params$rho_within,
    selector          = selector_method,
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

# ---------------------------- 7. 优先级 1（已改为每 bin 30 次） ----------------------------
run_priority1 <- function(n_seeds_per_bin = 30,          # ← 从 5 改成 30
                          leak_bins = list(c(0.0, 0.2), c(0.2, 0.4),
                                           c(0.4, 0.6), c(0.6, 0.8),
                                           c(0.8, 1.0)),
                          n_iter = 10, shap_nsim = 20,
                          seed_start = 20000) {
  set.seed(seed_start)
  all_results <- list()
  idx <- 1
  for (rho in c(0, 0.85)) {
    for (bin in leak_bins) {
      for (k in seq_len(n_seeds_per_bin)) {
        s <- seed_start + idx; idx <- idx + 1
        leak_i <- runif(1, bin[1], bin[2])
        n_i <- sample(c(400, 600), 1)
        top_k_i <- sample(c(3, 4), 1)
        sp <- list(n = n_i, n_groups = 4, features_per_group = 5,
                   rho_within = rho, missing_rate = 0.1, label_noise = 0.1,
                   leak_strength = leak_i, top_k = top_k_i, seed = s)
        res <- tryCatch(
          run_simulation(sp, n_iter = n_iter, shap_nsim = shap_nsim,
                         selector_method = "lasso"),
          error = function(e) NULL
        )
        if (!is.null(res)) {
          res$rho_label <- ifelse(rho == 0, "No redundancy", "With redundancy")
          res$leak_bin <- sprintf("%.1f-%.1f", bin[1], bin[2])
          all_results[[length(all_results) + 1]] <- res
        }
      }
      cat("  P1: rho =", rho, ", bin =", bin[1], "-", bin[2], "done\n")
    }
  }
  bind_rows(all_results)
}

# ---------------------------- 8. 优先级 2 ----------------------------
run_priority2 <- function(n_seeds = 10, n_iter = 10, shap_nsim = 20,
                          seed_start = 30000) {
  set.seed(seed_start)
  all_results <- list()
  idx <- 1
  for (sel in c("lasso", "univariate", "elastic_net")) {
    for (k in seq_len(n_seeds)) {
      s <- seed_start + idx; idx <- idx + 1
      leak_i <- runif(1, 0.3, 0.7)
      n_i <- sample(c(400, 600), 1)
      top_k_i <- sample(c(3, 4), 1)
      sp <- list(n = n_i, n_groups = 4, features_per_group = 5,
                 rho_within = 0.85, missing_rate = 0.1, label_noise = 0.1,
                 leak_strength = leak_i, top_k = top_k_i, seed = s)
      res <- tryCatch(
        run_simulation(sp, n_iter = n_iter, shap_nsim = shap_nsim,
                       selector_method = sel),
        error = function(e) NULL
      )
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
    cat("  P2: selector =", sel, "done\n")
  }
  bind_rows(all_results)
}

# ============================================================
# 执行实验
# ============================================================
if (RUN_P1) {
  cat("\n========== P1: SHAP drop vs leakage strength (n=30 per bin) ==========\n")
  p1_df <- run_priority1(n_seeds_per_bin = 30)   # ← 关键修改
  saveRDS(p1_df, "priority1_results_n30.rds")
} else if (file.exists("priority1_results_n15.rds")) {
  p1_df <- readRDS("priority1_results_n15.rds")
  cat("P1 loaded from priority1_results_n15.rds\n")
} else {
  stop("RUN_P1 = FALSE but priority1_results_n15.rds not found.")
}

if (RUN_P2) {
  cat("\n========== P2: Selector comparison ==========\n")
  p2_df <- run_priority2(n_seeds = 10)
  saveRDS(p2_df, "priority2_results.rds")
} else if (file.exists("priority2_results.rds")) {
  p2_df <- readRDS("priority2_results.rds")
  cat("P2 loaded from priority2_results.rds\n")
} else {
  stop("RUN_P2 = FALSE but priority2_results.rds not found.")
}

# ============================================================
# 正式统计检验（新增）
# ============================================================
cat("\n========== Formal statistical tests for inverted-U ==========\n")

# 只分析有冗余的组
with_red <- p1_df %>% filter(rho_label == "With redundancy")

# 1. 把 leak_bin 转成数值中点，方便做回归
with_red <- with_red %>%
  mutate(
    leak_mid = case_when(
      leak_bin == "0.0-0.2" ~ 0.1,
      leak_bin == "0.2-0.4" ~ 0.3,
      leak_bin == "0.4-0.6" ~ 0.5,
      leak_bin == "0.6-0.8" ~ 0.7,
      leak_bin == "0.8-1.0" ~ 0.9
    )
  )

# 2. 二次多项式回归：shap_drop ~ leak_mid + I(leak_mid^2)
quad_model <- lm(shap_drop ~ leak_mid + I(leak_mid^2), data = with_red)
cat("\n----- Quadratic polynomial regression -----\n")
print(summary(quad_model))

# 提取二次项 p 值
quad_p <- summary(quad_model)$coefficients["I(leak_mid^2)", "Pr(>|t|)"]
cat("\nQuadratic term p-value:", format.pval(quad_p, digits = 4), "\n")

# 3. 峰值 bin (0.4-0.6) 与两侧的直接对比
peak   <- with_red %>% filter(leak_bin == "0.4-0.6") %>% pull(shap_drop)
left   <- with_red %>% filter(leak_bin == "0.2-0.4") %>% pull(shap_drop)
right  <- with_red %>% filter(leak_bin == "0.6-0.8") %>% pull(shap_drop)

# t 检验（也可改成 wilcox.test）
t_left  <- t.test(peak, left)
t_right <- t.test(peak, right)

cat("\n----- Peak (0.4-0.6) vs adjacent bins -----\n")
cat("Peak vs 0.2-0.4:\n")
cat("  mean difference =", round(mean(peak) - mean(left), 4), "\n")
cat("  t =", round(t_left$statistic, 3), ", p =", format.pval(t_left$p.value, digits = 4), "\n")
cat("  Cohen's d ≈", round((mean(peak) - mean(left)) / sqrt((var(peak)+var(left))/2), 3), "\n")

cat("\nPeak vs 0.6-0.8:\n")
cat("  mean difference =", round(mean(peak) - mean(right), 4), "\n")
cat("  t =", round(t_right$statistic, 3), ", p =", format.pval(t_right$p.value, digits = 4), "\n")
cat("  Cohen's d ≈", round((mean(peak) - mean(right)) / sqrt((var(peak)+var(right))/2), 3), "\n")

# 4. 简单汇总表（方便写进论文）
p1_summary_stat <- with_red %>%
  group_by(leak_bin) %>%
  summarise(
    n = n(),
    mean = mean(shap_drop, na.rm = TRUE),
    se = sd(shap_drop, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
cat("\n----- Summary by bin (With redundancy) -----\n")
print(p1_summary_stat)

# ============================================================
# Publication-Quality Figures（保持原样，仅更新 subtitle）
# ============================================================
library(tidyverse)
library(ggplot2)

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title       = element_text(size = base_size + 1),
      axis.text        = element_text(size = base_size, colour = "black"),
      axis.line        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = element_line(colour = "black", linewidth = 0.4),
      legend.title     = element_text(size = base_size),
      legend.text      = element_text(size = base_size - 1),
      legend.key.size  = unit(0.4, "cm"),
      legend.position  = "top",
      plot.title       = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle    = element_text(size = base_size, colour = "grey30"),
      plot.margin      = margin(6, 8, 6, 6)
    )
}

pal_2 <- c("No redundancy" = "#D55E00", "With redundancy" = "#0072B2")
pal_3 <- c("LASSO" = "#0072B2", "Elastic Net" = "#009E73",
           "Univariate" = "#D55E00")

p1_df <- p1_df %>%
  mutate(rho_label = dplyr::recode(as.character(rho_label),
                                   "无冗余" = "No redundancy",
                                   "有冗余" = "With redundancy",
                                   .default = as.character(rho_label)))

p2_df <- p2_df %>%
  mutate(selector = dplyr::recode(as.character(selector),
                                  "lasso"        = "LASSO",
                                  "elastic_net"  = "Elastic Net",
                                  "univariate"   = "Univariate",
                                  .default = as.character(selector)))

# Figure 1
p1_summary <- p1_df %>%
  group_by(rho_label, leak_bin) %>%
  summarise(
    n     = n(),
    mean_drop = mean(shap_drop, na.rm = TRUE),
    se        = sd(shap_drop, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(leak_bin = factor(leak_bin,
                           levels = c("0.0-0.2", "0.2-0.4",
                                      "0.4-0.6", "0.6-0.8", "0.8-1.0")))

fig1 <- ggplot(p1_summary,
               aes(x = leak_bin, y = mean_drop,
                   colour = rho_label, group = rho_label)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean_drop - 1.96 * se,
                    ymax = mean_drop + 1.96 * se),
                width = 0.10, linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = pal_2, name = NULL) +
  scale_y_continuous(limits = c(-0.15, 0.65),
                     breaks = seq(-0.1, 0.6, by = 0.1)) +
  labs(
    x        = "Leakage strength (fraction of test set)",
    y        = expression("SHAP stability drop "*Delta*rho[S]),
    title    = "(A) Redundancy gates leakage-induced SHAP instability",
    subtitle = "Error bars: 95% CI; n = 30 per bin"   # ← 更新
  ) +
  theme_pub()

# Figure 2
p2_df <- p2_df %>%
  mutate(selector = factor(selector,
                           levels = c("LASSO", "Elastic Net", "Univariate")))

fig2 <- ggplot(p2_df, aes(x = selector, y = shap_drop, fill = selector)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_boxplot(width = 0.55, outlier.shape = NA,
               colour = "black", linewidth = 0.4, alpha = 0.7) +
  geom_jitter(width = 0.10, size = 1.4, alpha = 0.75,
              colour = "grey15") +
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 2.6, fill = "white",
               colour = "black", stroke = 0.5) +
  scale_fill_manual(values = pal_3, guide = "none") +
  scale_y_continuous(limits = c(-0.25, 0.85),
                     breaks = seq(-0.2, 0.8, by = 0.2)) +
  labs(
    x        = "Feature selector",
    y        = expression("SHAP stability drop "*Delta*rho[S]),
    title    = "(B) Instability across selectors",
    subtitle = "Diamonds: means; n = 10 per selector"
  ) +
  theme_pub()

# 保存
save_pdf <- function(plot, filename, w, h) {
  ok <- tryCatch({
    ggsave(filename, plot, width = w, height = h, device = cairo_pdf)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    ggsave(filename, plot, width = w, height = h)
  }
}

save_pdf(fig1, "fig1_shap_drop_vs_leakage_n30.pdf", 6.5, 4.0)
ggsave("fig1_shap_drop_vs_leakage_n30.png", fig1,
       width = 6.5, height = 4.0, dpi = 600)

save_pdf(fig2, "fig2_selector_comparison.pdf", 5.5, 4.0)
ggsave("fig2_selector_comparison.png", fig2,
       width = 5.5, height = 4.0, dpi = 600)

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined <- fig1 / fig2 + plot_layout(heights = c(1, 1))
  save_pdf(combined, "fig_combined_n30.pdf", 9, 8)
  ggsave("fig_combined_30.png", combined,
         width = 9, height = 8, dpi = 600)
}

# ---------------------------- 最终汇总 ----------------------------
cat("\n============================================\n")
cat("DONE. Output directory:\n")
cat(" ", getwd(), "\n")
cat("\nFiles:\n")
cat("  priority1_results_n30.rds\n")
cat("  priority2_results.rds\n")
cat("  fig1_shap_drop_vs_leakage_n30.pdf / .png\n")
cat("  fig2_selector_comparison.pdf / .png\n")
cat("============================================\n")

cat("\n========== P1 Summary (n=15) ==========\n")
p1_df %>%
  group_by(rho_label, leak_bin) %>%
  summarise(
    n = n(),
    mean_drop = mean(shap_drop, na.rm = TRUE),
    se        = sd(shap_drop, na.rm = TRUE) / sqrt(n()),
    mean_abs_dauc = mean(abs(delta_auc), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = 20)

cat("\n========== P2 Summary ==========\n")
p2_df %>%
  group_by(selector) %>%
  summarise(
    n = n(),
    n_valid = sum(is.finite(shap_drop)),
    mean_drop = mean(shap_drop, na.rm = TRUE),
    sd_drop   = sd(shap_drop, na.rm = TRUE),
    mean_abs_dauc = mean(abs(delta_auc), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()