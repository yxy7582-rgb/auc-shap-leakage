# ============================================================
# Real-world validation: AUC-blind + SHAP-collapse
# Fixed 2x2 factorial design (proper A/B/C/D conditions)
# ============================================================

# ---------- 0. Packages ----------
pkgs <- c(
  "readxl", "dplyr", "tidyr", "tibble", "stringr", "purrr",
  "caret", "glmnet", "ranger", "fastshap", "ggplot2",
  "patchwork", "janitor", "pROC"
)
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
invisible(lapply(pkgs, library, character.only = TRUE))
set.seed(20240924)

# ---------- 0b. Publication theme ----------
theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      text             = element_text(family = "sans", colour = "black"),
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      axis.title       = element_text(size = base_size, colour = "black"),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      axis.line        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = element_line(colour = "black", linewidth = 0.4),
      panel.border     = element_rect(colour = "black", linewidth = 0.5, fill = NA),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.title     = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      legend.key       = element_blank(),
      legend.background= element_blank(),
      legend.position  = "right"
    )
}
save_fig <- function(p, filename, width_mm = 90, height_mm = 70) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(filename, ".pdf"), p, width = w, height = h,
         units = "in", device = cairo_pdf)
  ggsave(paste0(filename, ".tiff"), p, width = w, height = h,
         units = "in", dpi = 600, compression = "lzw")
  ggsave(paste0(filename, ".png"), p, width = w, height = h,
         units = "in", dpi = 300)
}

# ---------- 1. Load data ----------
file_path <- "your file"
raw <- read_excel(file_path, sheet = 1, skip = 1, .name_repair = "minimal")

colnames(raw) <- str_trim(as.character(colnames(raw)))
colnames(raw)[colnames(raw) == "" | is.na(colnames(raw))] <-
  paste0("V", which(colnames(raw) == "" | is.na(colnames(raw))))
colnames(raw) <- make.unique(colnames(raw), sep = "_")

cat("=== 前25个列名 ===\n")
print(head(colnames(raw), 25))

ild_col   <- "合并间质性肺病（需要预测的结果）"
malig_col <- "合并恶性肿瘤（需要预测的结果）"

if (!ild_col %in% colnames(raw)) {
  ild_candidates <- grep("间质", colnames(raw), value = TRUE)
  if (length(ild_candidates) == 0) stop("找不到 ILD 结局列，请检查列名")
  ild_col <- ild_candidates[1]
  cat("使用模糊匹配的 ILD 列：", ild_col, "\n")
}
if (!malig_col %in% colnames(raw)) {
  malig_candidates <- grep("恶性", colnames(raw), value = TRUE)
  malig_candidates <- malig_candidates[!grepl("类别", malig_candidates)]
  malig_col <- if (length(malig_candidates) > 0) malig_candidates[1] else NA_character_
}

cat("ILD 列：", ild_col, "\n")
cat("Malignancy 列：", malig_col, "\n")

id_leak_cols <- c(
  "患者主索引", "就诊号", "姓名",
  "恶性肿瘤类别（大类）", "心脑血管疾病类别"
)
id_leak_cols <- intersect(id_leak_cols, colnames(raw))
df <- raw %>% select(-any_of(id_leak_cols))

clean_val <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% c("", "NA", "N/A", "NULL", "建议CKMB质量检测", "None")] <- NA
  x <- str_replace_all(x, "^[<>≤≥]", "")
  x
}
df <- df %>% mutate(across(everything(), clean_val))

df[[ild_col]] <- suppressWarnings(as.numeric(df[[ild_col]]))
if (!is.na(malig_col) && malig_col %in% colnames(df)) {
  df[[malig_col]] <- suppressWarnings(as.numeric(df[[malig_col]]))
}

# ---------- 2. Build ILD dataset ----------
df_ild <- df %>% filter(!is.na(.data[[ild_col]]))
if (!is.na(malig_col) && malig_col %in% colnames(df_ild)) {
  df_ild <- df_ild %>% select(-all_of(malig_col))
}

df_ild <- df_ild %>% mutate(across(where(is.character),
                                   ~ suppressWarnings(as.numeric(.))))
df_ild <- df_ild %>% select(where(~ !all(is.na(.))))

miss_rate <- sapply(df_ild, function(x) mean(is.na(x)))
df_ild <- df_ild[, miss_rate <= 0.50, drop = FALSE]

nzv <- nearZeroVar(df_ild, saveMetrics = TRUE)
df_ild <- df_ild[, !(nzv$nzv | nzv$zeroVar), drop = FALSE]

y <- df_ild[[ild_col]]
X <- df_ild %>% select(-all_of(ild_col))
X <- X[, sapply(X, is.numeric), drop = FALSE]
X <- X[, apply(X, 2, function(x) {
  s <- suppressWarnings(sd(x, na.rm = TRUE))
  !is.na(s) && s > 0
}), drop = FALSE]

cat("\nILD 最终数据集：n =", nrow(X), "  p =", ncol(X),
    "  阳性数 =", sum(y == 1, na.rm = TRUE), "\n")

if (nrow(X) < 50 || ncol(X) < 5) stop("数据量太小，请检查过滤条件")

# ---------- 3. Core functions ----------
preprocess_fit <- function(train_df) {
  medians <- sapply(train_df, function(x) median(x, na.rm = TRUE))
  means   <- sapply(train_df, function(x) mean(x, na.rm = TRUE))
  sds     <- sapply(train_df, function(x) sd(x, na.rm = TRUE))
  sds[sds == 0 | is.na(sds)] <- 1
  medians[is.na(medians)] <- 0
  means[is.na(means)] <- 0
  list(medians = medians, means = means, sds = sds)
}

preprocess_apply <- function(df, params) {
  df_imp <- df
  common <- intersect(names(df_imp), names(params$medians))
  for (j in common) {
    df_imp[[j]][is.na(df_imp[[j]])] <- params$medians[j]
  }
  as.data.frame(scale(df_imp,
                      center = params$means[names(df_imp)],
                      scale  = params$sds[names(df_imp)]))
}

select_features <- function(train_x, train_y, selector = "lasso",
                            top_k = 5, seed = 1) {
  set.seed(seed)
  all_feats <- colnames(train_x)
  top_k <- max(1, min(top_k, length(all_feats)))
  if (length(unique(train_y)) < 2) return(all_feats[1:top_k])
  
  if (selector %in% c("lasso", "elastic_net")) {
    alpha <- ifelse(selector == "lasso", 1, 0.5)
    cvfit <- tryCatch(
      cv.glmnet(as.matrix(train_x), train_y, family = "binomial",
                alpha = alpha, nfolds = 5, parallel = FALSE),
      error = function(e) NULL)
    if (is.null(cvfit)) return(all_feats[1:top_k])
    coefs <- as.matrix(coef(cvfit, s = "lambda.min"))
    selected <- setdiff(rownames(coefs)[coefs[, 1] != 0], "(Intercept)")
    if (length(selected) < top_k) {
      cors <- abs(sapply(train_x, function(x)
        suppressWarnings(cor(x, train_y, use = "complete.obs"))))
      cors[is.na(cors)] <- 0
      backfill <- setdiff(names(sort(cors, decreasing = TRUE)), selected)
      selected <- c(selected, head(backfill, top_k - length(selected)))
    }
    selected <- unique(selected)[1:min(top_k, length(selected))]
  } else {
    cors <- abs(sapply(train_x, function(x)
      suppressWarnings(cor(x, train_y, use = "complete.obs"))))
    cors[is.na(cors)] <- 0
    selected <- names(sort(cors, decreasing = TRUE))[1:top_k]
  }
  selected[!is.na(selected) & selected %in% all_feats]
}

fit_rf_shap <- function(train_x, train_y, test_x, test_y,
                        selected, seed = 1, nsim = 15) {
  set.seed(seed)
  if (length(selected) < 2) {
    extra <- setdiff(colnames(train_x), selected)
    if (length(extra) > 0) selected <- c(selected, extra[1])
  }
  train_sel <- train_x[, selected, drop = FALSE]
  test_sel  <- test_x[, selected, drop = FALSE]
  
  rf <- ranger(x = train_sel, y = as.factor(train_y),
               num.trees = 150, probability = TRUE,
               seed = seed, num.threads = 1)
  
  auc <- tryCatch({
    pred <- predict(rf, data = test_sel)$predictions[, "1"]
    as.numeric(pROC::auc(roc(test_y, pred, quiet = TRUE)))
  }, error = function(e) NA_real_)
  
  shap_vals <- tryCatch(
    fastshap::explain(
      object = rf, X = train_sel, nsim = nsim,
      pred_wrapper = function(object, newdata) {
        predict(object, data = newdata)$predictions[, "1"]
      }, seed = seed),
    error = function(e) NULL)
  
  imp_full <- setNames(rep(0, ncol(train_x)), colnames(train_x))
  if (!is.null(shap_vals)) {
    imp <- colMeans(abs(as.matrix(shap_vals)))
    imp_full[names(imp)] <- imp
  }
  list(auc = auc, imp = imp_full, selected = selected)
}

run_one_pipeline <- function(X, y, train_idx, test_idx,
                             leak_strength = 0, selector = "lasso",
                             top_k = 5, seed = 1, leaky = FALSE) {
  set.seed(seed)
  train_x <- X[train_idx, , drop = FALSE]
  train_y <- y[train_idx]
  test_x  <- X[test_idx, , drop = FALSE]
  test_y  <- y[test_idx]
  
  if (leaky && leak_strength > 0) {
    n_leak <- max(1, round(leak_strength * length(test_idx)))
    n_leak <- min(n_leak, length(test_idx) - 5)
    if (n_leak < 1) n_leak <- 1
    leak_idx <- sample(test_idx, n_leak)
    merge_x <- rbind(train_x, X[leak_idx, , drop = FALSE])
    merge_y <- c(train_y, y[leak_idx])
    remain  <- setdiff(test_idx, leak_idx)
    test_x  <- X[remain, , drop = FALSE]
    test_y  <- y[remain]
  } else {
    merge_x <- train_x
    merge_y <- train_y
  }
  
  if (length(unique(merge_y)) < 2 || length(unique(test_y)) < 2) {
    return(list(auc = NA_real_,
                imp = setNames(rep(0, ncol(X)), colnames(X)),
                selected = character(0)))
  }
  
  params    <- preprocess_fit(merge_x)
  merge_std <- preprocess_apply(merge_x, params)
  test_std  <- preprocess_apply(test_x, params)
  selected  <- select_features(merge_std, merge_y, selector, top_k, seed)
  
  fit_rf_shap(merge_std, merge_y, test_std, test_y, selected, seed = seed)
}

shap_stability <- function(imp_matrix) {
  if (is.null(imp_matrix) || nrow(imp_matrix) < 2) return(NA_real_)
  cors <- c()
  for (i in 1:(nrow(imp_matrix) - 1)) {
    for (j in (i + 1):nrow(imp_matrix)) {
      cval <- suppressWarnings(
        cor(imp_matrix[i, ], imp_matrix[j, ], method = "spearman"))
      if (!is.na(cval)) cors <- c(cors, cval)
    }
  }
  if (length(cors) == 0) return(NA_real_)
  mean(cors)
}

run_2x2 <- function(X_use, y, train_idx, test_idx, leak_strength = 0.5,
                    selector = "lasso", top_k = 5, K = 8, seed_base = 100) {
  run_k <- function(leaky, ls) {
    res_list <- lapply(1:K, function(k) {
      run_one_pipeline(X_use, y, train_idx, test_idx,
                       leak_strength = ls, selector = selector,
                       top_k = top_k, seed = seed_base + k * 17,
                       leaky = leaky)
    })
    imp_mat <- do.call(rbind, lapply(res_list, `[[`, "imp"))
    aucs    <- sapply(res_list, `[[`, "auc")
    list(stab = shap_stability(imp_mat),
         auc  = mean(aucs, na.rm = TRUE))
  }
  clean <- run_k(FALSE, 0)
  leaky <- run_k(TRUE,  leak_strength)
  data.frame(
    clean_auc  = clean$auc,
    leaky_auc  = leaky$auc,
    delta_auc  = abs(leaky$auc - clean$auc),
    clean_stab = clean$stab,
    leaky_stab = leaky$stab,
    shap_drop  = clean$stab - leaky$stab
  )
}

# ---------- 4. Redundancy audit + reduced X ----------
set.seed(42)
train_idx <- createDataPartition(y, p = 0.70, list = FALSE)[, 1]
test_idx  <- setdiff(seq_along(y), train_idx)

cor_mat <- cor(X, method = "spearman", use = "pairwise.complete.obs")
cor_mat[is.na(cor_mat)] <- 0
hc <- hclust(as.dist(1 - abs(cor_mat)), method = "average")
groups <- cutree(hc, h = 0.30)   # |r| > 0.7
max_group <- max(table(groups))
cat("最大冗余组大小：", max_group, "\n")

# ---- KEY FIX: top_k < max_group to trigger within-group competition ----
top_k <- max(2, min(max_group - 1, 5))
cat("使用 top_k =", top_k, "（< max_group =", max_group, "）\n")

# ---- KEY FIX: build truly non-redundant X_reduced ----
# Keep only one representative feature per redundancy group
keep_one_per_group <- function(groups, X) {
  grp <- as.integer(groups)
  names(grp) <- names(groups)
  keep <- tapply(names(grp), grp, function(x) x[1])
  keep <- unlist(keep, use.names = FALSE)
  X[, intersect(keep, colnames(X)), drop = FALSE]
}
X_reduced <- keep_one_per_group(groups, X)
cat("X_reduced 特征数：", ncol(X_reduced),
    "（原 X：", ncol(X), "）\n")

# Sanity check: redundancy in X_reduced
cor_red <- cor(X_reduced, method = "spearman", use = "pairwise.complete.obs")
cor_red[is.na(cor_red)] <- 0
groups_red <- cutree(hclust(as.dist(1 - abs(cor_red)), method = "average"),
                     h = 0.30)
cat("X_reduced 最大组大小：", max(table(groups_red)), "\n")

# ---------- 5. Proper 2x2 factorial ----------
cat("\n开始运行 2×2 实验...\n")

# A: 无冗余 + 泄漏
res_A <- run_2x2(X_reduced, y, train_idx, test_idx,
                 leak_strength = 0.5, top_k = top_k, K = 8, seed_base = 100)
# B: 无冗余 + 无泄漏
res_B <- run_2x2(X_reduced, y, train_idx, test_idx,
                 leak_strength = 0,   top_k = top_k, K = 8, seed_base = 200)
# C: 有冗余 + 泄漏
res_C <- run_2x2(X,         y, train_idx, test_idx,
                 leak_strength = 0.5, top_k = top_k, K = 8, seed_base = 300)
# D: 有冗余 + 无泄漏
res_D <- run_2x2(X,         y, train_idx, test_idx,
                 leak_strength = 0,   top_k = top_k, K = 8, seed_base = 400)

df_2x2 <- rbind(
  cbind(condition = "A: no redundancy + leakage",    res_A),
  cbind(condition = "B: no redundancy + no leakage", res_B),
  cbind(condition = "C: redundancy + leakage",       res_C),
  cbind(condition = "D: redundancy + no leakage",    res_D)
)

cat("\n===== 2×2 结果（ILD） =====\n")
print(df_2x2)

# ---- Interaction term ----
interaction_term <- (df_2x2$shap_drop[df_2x2$condition == "C: redundancy + leakage"] -
                       df_2x2$shap_drop[df_2x2$condition == "D: redundancy + no leakage"]) -
  (df_2x2$shap_drop[df_2x2$condition == "A: no redundancy + leakage"] -
     df_2x2$shap_drop[df_2x2$condition == "B: no redundancy + no leakage"])
cat("\nInteraction term (SHAP drop) =", round(interaction_term, 4), "\n")

write.csv(df_2x2, "real_data_2x2_ild_fixed.csv", row.names = FALSE)

# ---- Figure: 2x2 factorial ----
df_2x2_plot <- df_2x2 %>%
  mutate(condition = factor(
    condition,
    levels = c("B: no redundancy + no leakage",
               "A: no redundancy + leakage",
               "D: redundancy + no leakage",
               "C: redundancy + leakage"),
    labels = c("No redundancy,\nno leakage",
               "No redundancy,\nleakage",
               "Redundancy,\nno leakage",
               "Redundancy,\nleakage")))

fig2a <- ggplot(df_2x2_plot, aes(x = condition, y = shap_drop)) +
  geom_col(fill = "grey70", colour = "black", linewidth = 0.4, width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(x = NULL, y = "SHAP stability drop  (\u0394\u03C1)") +
  theme_pub(base_size = 9) +
  theme(axis.text.x = element_text(size = 8))

fig2b <- ggplot(df_2x2_plot, aes(x = condition, y = delta_auc)) +
  geom_col(fill = "grey70", colour = "black", linewidth = 0.4, width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(x = NULL, y = "AUC change  (|\u0394AUC|)") +
  theme_pub(base_size = 9) +
  theme(axis.text.x = element_text(size = 8))

fig2 <- fig2a | fig2b
print(fig2)
save_fig(fig2, "Fig2_2x2_factorial", width_mm = 180, height_mm = 70)

# ---------- 6. Leakage-strength sweep ----------
cat("\n开始泄漏强度扫描...\n")
leak_levels <- c(0, 0.2, 0.4, 0.6, 0.8)
n_rep  <- 4
K_impl <- 5
results_sweep <- list()

for (rep in 1:n_rep) {
  set.seed(1000 + rep)
  tr <- createDataPartition(y, p = 0.70, list = FALSE)[, 1]
  te <- setdiff(seq_along(y), tr)
  
  for (ls in leak_levels) {
    clean_res <- lapply(1:K_impl, function(k) {
      run_one_pipeline(X, y, tr, te, leak_strength = 0, top_k = top_k,
                       seed = rep * 100 + k, leaky = FALSE)
    })
    leaky_res <- lapply(1:K_impl, function(k) {
      run_one_pipeline(X, y, tr, te, leak_strength = ls, top_k = top_k,
                       seed = rep * 100 + k, leaky = TRUE)
    })
    
    clean_imp <- do.call(rbind, lapply(clean_res, `[[`, "imp"))
    leaky_imp <- do.call(rbind, lapply(leaky_res, `[[`, "imp"))
    clean_auc <- mean(sapply(clean_res, `[[`, "auc"), na.rm = TRUE)
    leaky_auc <- mean(sapply(leaky_res, `[[`, "auc"), na.rm = TRUE)
    
    results_sweep[[length(results_sweep) + 1]] <- data.frame(
      rep = rep, leak_strength = ls,
      clean_auc = clean_auc, leaky_auc = leaky_auc,
      delta_auc = abs(leaky_auc - clean_auc),
      clean_stab = shap_stability(clean_imp),
      leaky_stab = shap_stability(leaky_imp),
      shap_drop  = shap_stability(clean_imp) - shap_stability(leaky_imp)
    )
  }
  cat("完成 repeat", rep, "/", n_rep, "\n")
}

df_sweep <- bind_rows(results_sweep)
summary_sweep <- df_sweep %>%
  group_by(leak_strength) %>%
  summarise(
    mean_delta_auc = mean(delta_auc, na.rm = TRUE),
    se_delta_auc   = sd(delta_auc, na.rm = TRUE) / sqrt(n()),
    mean_shap_drop = mean(shap_drop, na.rm = TRUE),
    se_shap_drop   = sd(shap_drop, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

cat("\n===== 泄漏强度扫描汇总 =====\n")
print(summary_sweep)

# Ideal separation
df_sweep <- df_sweep %>%
  mutate(ideal_strict = (delta_auc < 0.03) & (shap_drop > 0.15),
         ideal_loose  = (delta_auc < 0.05) & (shap_drop > 0.15))
cat("Ideal separation (strict, |dAUC|<0.03):",
    round(mean(df_sweep$ideal_strict, na.rm = TRUE) * 100, 1), "%\n")
cat("Ideal separation (loose,  |dAUC|<0.05):",
    round(mean(df_sweep$ideal_loose,  na.rm = TRUE) * 100, 1), "%\n")

# ---- Figure: leakage sweep ----
plot_sweep_long <- summary_sweep %>%
  select(leak_strength, mean_shap_drop, se_shap_drop,
         mean_delta_auc, se_delta_auc) %>%
  pivot_longer(cols = c(mean_shap_drop, mean_delta_auc),
               names_to = "metric", values_to = "mean") %>%
  mutate(
    se = ifelse(metric == "mean_shap_drop", se_shap_drop, se_delta_auc),
    metric_label = factor(
      metric,
      levels = c("mean_shap_drop", "mean_delta_auc"),
      labels = c("SHAP stability drop  (\u0394\u03C1)",
                 "AUC change  (|\u0394AUC|)")))

fig1 <- ggplot(plot_sweep_long,
               aes(x = leak_strength, y = mean,
                   colour = metric_label, shape = metric_label,
                   linetype = metric_label)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.03, linewidth = 0.5) +
  scale_colour_manual(values = c(
    "SHAP stability drop  (\u0394\u03C1)" = "#1F78B4",
    "AUC change  (|\u0394AUC|)"           = "#E66101")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
  labs(x = "Leakage strength (fraction of test set)",
       y = "Metric value",
       colour = NULL, shape = NULL, linetype = NULL) +
  theme_pub(base_size = 10) +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = "white", colour = NA))

print(fig1)
save_fig(fig1, "Fig1_leakage_sweep", width_mm = 90, height_mm = 65)

# ---- Figure: ideal separation scatter ----
df_sweep_plot <- df_sweep %>%
  mutate(leak_f = factor(leak_strength,
                         levels = c(0, 0.2, 0.4, 0.6, 0.8),
                         labels = c("0.0", "0.2", "0.4", "0.6", "0.8")))

fig4 <- ggplot(df_sweep_plot,
               aes(x = delta_auc, y = shap_drop, colour = leak_f)) +
  geom_point(size = 1.8, alpha = 0.85) +
  geom_hline(yintercept = 0.15, linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = 0.03, linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  scale_colour_viridis_d(name = "Leakage\nstrength",
                         option = "D", end = 0.9) +
  labs(x = "AUC change  (|\u0394AUC|)",
       y = "SHAP stability drop  (\u0394\u03C1)") +
  annotate("text", x = 0.005, y = 0.35,
           label = "Ideal separation",
           hjust = 0, vjust = 1, size = 3, colour = "black") +
  theme_pub(base_size = 10)

print(fig4)
save_fig(fig4, "Fig4_ideal_separation", width_mm = 90, height_mm = 75)

write.csv(df_sweep,      "real_data_leak_sweep.csv",         row.names = FALSE)
write.csv(summary_sweep, "real_data_leak_sweep_summary.csv", row.names = FALSE)

cat("\n全部完成！\n")
cat("关键检查：\n")
cat("  1) A vs C 的 shap_drop 差异 —— 冗余是否是前提\n")
cat("  2) A 和 B 的 shap_drop 是否接近 0 —— 无冗余时泄漏是否无效\n")
cat("  3) interaction term 是否接近 C 的 shap_drop —— 是否纯交互\n")