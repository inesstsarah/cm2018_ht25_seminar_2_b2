# ================================================================
# TASK 1 !!
# ================================================================

# IMPORTANT: Running the code automatically saves the plots into your computer

# ================================================================
# Histograms + binary counts + summary statistics
# ================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(scales)
library(rlang)   
theme_set(theme_minimal(base_size = 12))

data <- read_csv("/Users/claraazzano/Desktop/Data_T1.csv", show_col_types = FALSE)

# ================================================================
# Helpers
# ================================================================

# Freedman–Diaconis binwidth with fallback
fd_binwidth <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  bw <- 2 * IQR(x) / (length(x)^(1/3))
  if (!is.finite(bw) || bw <= 0) {
    rng <- diff(range(x))
    bw <- if (is.finite(rng) && rng > 0) rng / 30 else 1
  }
  bw
}

# Proper skew flag: third standardized moment
is_skewed <- function(x, thresh = 1) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || sd(x) == 0) return(FALSE)
  z  <- (x - mean(x)) / sd(x)
  m3 <- mean(z^3)
  abs(m3) > thresh
}

# Safe log10 (adds small offset if <=0 present)
safe_log10 <- function(x, eps = 1e-6) {
  x <- as.numeric(x)
  needs_offset <- any(x <= 0, na.rm = TRUE)
  off <- if (needs_offset) max(eps, abs(min(x, na.rm = TRUE)) + eps) else 0
  y <- log10(x + off)
  attr(y, "offset") <- off
  y
}

# Axis label with optional units and (log10) note
axis_label <- function(name, log = FALSE, units = NULL) {
  paste0(
    name,
    if (!is.null(units) && nzchar(units)) paste0(" [", units, "]") else "",
    if (log) " (log10)" else ""
  )
}

# ================================================================
# Variables and units
# ================================================================

cont_vars <- c(
  "GE","Age","Height","BW","BMI",
  "GlucoseFasting","InsulinFasting","HbA1c","MatsudaIdx","HOMAB",
  "Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY"
)

# candidates that are typically better on log-scale for EDA
log_suggest <- c("InsulinFasting","MatsudaIdx","Ghrelin","Amylin","GLP1")

binary_vars <- c("Sex","Metformin","DiabetesComplications")

# Units map (keys MUST match your column names)
units_map <- c(
  GE = "min", Age = "y", Height = "cm", BW = "kg", BMI = "kg/m^2",
  GlucoseFasting = "mmol/L", InsulinFasting = "mU/L", HbA1c = "%",
  Gastrin = "pg/mL", CCK = "pmol/L", Ghrelin = "pg/mL", Amylin = "pmol/L",
  Glucagon = "pmol/L", GLP1 = "pmol/L", PYY = "ng/L",
  MatsudaIdx = "", HOMAB = ""   # leave blank if unitless / index
)

# Safe unit lookup (prevents subscript-out-of-bounds)
get_units <- function(var) {
  nm <- as.character(var)
  if (!is.null(names(units_map)) && nm %in% names(units_map)) units_map[[nm]] else NULL
}

# ================================================================
# Plot histograms
# ================================================================

plot_hist <- function(df, var, bins = NULL, log_x = NULL,
                      units = NULL, title = NULL, annotate_n = TRUE) {
  v <- enquo(var)
  vname <- as_name(v)
  x <- df %>% pull(!!v)
  
  # auto-decide log if not specified
  if (is.null(log_x)) log_x <- all(x > 0, na.rm = TRUE) && is_skewed(x)
  
  x_plot <- if (log_x) safe_log10(x) else x
  off <- if (log_x) attr(x_plot, "offset") else 0
  
  df_plot <- tibble(x_plot = as.numeric(x_plot))
  
  # choose binning
  if (is.null(bins)) {
    bw <- fd_binwidth(df_plot$x_plot)
    p <- ggplot(df_plot, aes(x = x_plot)) +
      geom_histogram(binwidth = bw, fill = "#6A9AE2", color = "white", na.rm = TRUE)
  } else {
    p <- ggplot(df_plot, aes(x = x_plot)) +
      geom_histogram(bins = bins, fill = "#6A9AE2", color = "white", na.rm = TRUE)
  }
  
  # labels
  p <- p +
    labs(
      title = if (is.null(title)) vname else title,
      x = axis_label(vname, log = isTRUE(log_x), units = if (is.null(units)) get_units(vname) else units),
      y = "Count"
    ) +
    theme(panel.grid.minor = element_blank())
  
  # annotate n + log offset (if used)
  if (annotate_n) {
    n_use <- sum(is.finite(df_plot$x_plot))
    miss  <- sum(!is.finite(df_plot$x_plot))
    ann <- paste0("n=", n_use, if (miss > 0) paste0(" (NA=", miss, ")"))
    if (isTRUE(log_x) && off > 0) ann <- paste0(ann, "\nlog10(x + ", signif(off, 3), ")")
    p <- p + annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5, label = ann)
  }
  
  p
}

# ================================================================
# Draw histograms to screen
# ================================================================

invisible(lapply(cont_vars, function(v) {
  print(
    plot_hist(
      data, !!sym(v),
      log_x = if (v %in% log_suggest) TRUE else NULL,
      units = get_units(v)
    )
  )
}))

# ================================================================
# SAVE histograms
# ================================================================

outdir <- "task1_histograms"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

invisible(lapply(cont_vars, function(v) {
  g <- plot_hist(
    data, !!sym(v),
    log_x = if (v %in% log_suggest) TRUE else NULL,
    units = get_units(v)
  )
  ggsave(filename = file.path(outdir, paste0("hist_", v, ".png")),
         plot = g, width = 8, height = 5, dpi = 300)
}))

# ================================================================
# Binary summary TABLE (counts + percentages)
# ================================================================

binary_summary <- function(df, vars, labels = NULL) {
  out <- lapply(vars, function(v) {
    x <- df[[v]]
    tab <- tibble(Level = c(0,1)) %>%
      left_join(
        tibble(Level = x) |>
          filter(!is.na(Level)) |>
          count(Level, name = "n"),
        by = "Level"
      ) |>
      mutate(
        n   = replace_na(n, 0),
        pct = n / sum(n),
        Variable = v,
        Label = if (!is.null(labels) && !is.null(labels[[v]])) {
          labels[[v]][match(Level, c(0,1))]
        } else as.character(Level)
      ) |>
      select(Variable, Level, Label, n, pct)
  })
  bind_rows(out) |>
    mutate(pct_label = percent(pct, accuracy = 0.1))
}

label_map <- list(
  Sex = c("Female","Male"),
  Metformin = c("No","Yes"),
  DiabetesComplications = c("No","Yes")
)

bin_tbl <- binary_summary(data, binary_vars, labels = label_map)

# console view
bin_tbl %>%
  arrange(Variable, Level) %>%
  mutate(`n (%)` = paste0(n, " (", pct_label, ")")) %>%
  select(Variable, Label, `n (%)`) %>%
  print(n = Inf)

# wide CSV 
bin_tbl_wide <- bin_tbl %>%
  mutate(`n (%)` = paste0(n, " (", pct_label, ")")) %>%
  select(Variable, Label, `n (%)`) %>%
  pivot_wider(names_from = Label, values_from = `n (%)`)

write_csv(bin_tbl_wide, "binary_summary_task1.csv")

# ================================================================
# Descriptive stats for continuous + missingness
# ================================================================

desc <- data %>%
  summarise(across(all_of(cont_vars),
                   list(
                     n      = ~sum(!is.na(.)),
                     mean   = ~mean(., na.rm = TRUE),
                     sd     = ~sd(., na.rm = TRUE),
                     median = ~median(., na.rm = TRUE),
                     IQR    = ~IQR(., na.rm = TRUE),
                     min    = ~min(., na.rm = TRUE),
                     max    = ~max(., na.rm = TRUE)
                   ),
                   .names = "{.col}_{.fn}"
  )) %>%
  pivot_longer(everything(),
               names_to = c("variable",".value"),
               names_sep = "_") %>%
  mutate(across(c(mean, sd, median, IQR, min, max), ~round(., 2)))

write_csv(desc, "summary_task1.csv")

missing_tbl <- data %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to="variable", values_to="n_missing") %>%
  arrange(desc(n_missing))

print(missing_tbl, n = Inf)


# ======================================================================
# Correlations 
# ======================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(corrplot)
library(rlang)
library(patchwork)  

data_path <- "/Users/claraazzano/Desktop/Data_T1.csv"
out_dir   <- "/Users/claraazzano/correlations_task1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data <- read_csv(data_path, show_col_types = FALSE)

# ---- variable sets ----
vars_hormones <- c("Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY")
vars_demo     <- c("Age","Height","BW","BMI")
vars_gly      <- c("GlucoseFasting","HbA1c")
vars_insulin  <- c("InsulinFasting","MatsudaIdx","HOMAB")
vars_outcome  <- "GE"

# typical log-scale variables for correlations
log_skewed    <- c("InsulinFasting","MatsudaIdx","Ghrelin","Amylin","GLP1")

# ---- helpers ----
safe_log10 <- function(x, eps = 1e-6) {
  x <- as.numeric(x)
  needs_offset <- any(x <= 0, na.rm = TRUE)
  off <- if (needs_offset) max(eps, abs(min(x, na.rm = TRUE)) + eps) else 0
  y <- log10(x + off)
  attr(y, "offset") <- off
  y
}

make_pairs_df <- function(df, keep, log10_vars = character()) {
  keep <- intersect(keep, names(df))
  d <- df[keep]
  # keep numeric only; drop zero-variance columns
  d <- d[vapply(d, is.numeric, logical(1L))]
  for (v in intersect(log10_vars, names(d))) d[[v]] <- safe_log10(d[[v]])
  if (ncol(d)) {
    ok <- vapply(d, function(z) sd(z, na.rm = TRUE) > 0, logical(1L))
    d <- d[, ok, drop = FALSE]
  }
  d
}

# Spearman correlations (pairwise), BH-adjusted
spearman_table <- function(df, vars, log10_vars = character(),
                           sort_by = c("abs_r","p","r")) {
  d <- make_pairs_df(df, keep = vars, log10_vars = log10_vars)
  cn <- colnames(d); k <- length(cn)
  if (k < 2) return(tibble())
  out <- vector("list", k*(k-1)/2); z <- 0L
  for (i in seq_len(k-1)) for (j in (i+1):k) {
    x <- d[[i]]; y <- d[[j]]
    ok <- is.finite(x) & is.finite(y); n <- sum(ok)
    if (n >= 3) {
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      z <- z + 1L
      out[[z]] <- tibble(var1 = cn[i], var2 = cn[j],
                         rho = unname(ct$estimate), p = unname(ct$p.value), n = n)
    }
  }
  tab <- bind_rows(out)
  if (!nrow(tab)) return(tab)
  tab <- tab %>% mutate(padj = p.adjust(p, method = "BH"))
  sort_by <- match.arg(sort_by)
  ord <- switch(sort_by,
                abs_r = order(-abs(tab$rho), tab$padj, tab$p),
                p     = order(tab$p, -abs(tab$rho)),
                r     = order(-tab$rho))
  tab[ord, ]
}

# Base-pairs
pairs_plot <- function(df, vars, log10_vars = character(),
                       group = NULL, main = "") {
  d <- make_pairs_df(df, keep = vars, log10_vars = log10_vars)
  # custom panels
  .panel.density <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    d <- density(x[is.finite(x)])
    if (length(d$x)) lines(d$x, d$y / max(d$y), lwd = 1.1)
  }
  .panel.points <- function(x, y, col = "black", pch = 19, ...) {
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) points(x[ok], y[ok], pch = pch,
                        col = adjustcolor(col, alpha.f = 0.35), cex = 0.55)
  }
  .panel.shade <- function(x, y, palette = c("#B71C1C","white","#1B5E20"), ...) {
    ok <- is.finite(x) & is.finite(y)
    r  <- if (sum(ok) >= 3) suppressWarnings(cor(x[ok], y[ok], method = "spearman")) else NA_real_
    col_bg <- colorRampPalette(palette)(200)
    idx <- if (is.finite(r)) max(1, min(200, round((r + 1)/2 * 199 + 1))) else 100
    usr <- par("usr")
    rect(usr[1], usr[3], usr[2], usr[4], col = col_bg[idx], border = NA)
  }
  if (!is.null(group)) {
    g <- as.factor(group)
    cols <- adjustcolor(c("#6A9AE2", "#E28A6A", "#7BC47F", "#AA87FF"), alpha.f = 0.8)
    col_vec <- cols[as.integer(g)]
  } else col_vec <- "black"
  
  pairs(d,
        lower.panel = function(x, y, ...) .panel.points(x, y, col = col_vec),
        upper.panel = .panel.shade,
        diag.panel  = .panel.density,
        gap = 0.6, cex.labels = 1.0, main = main)
}

save_pairs <- function(filename, expr, width = 3200, height = 2400, res = 300) {
  png(filename, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

# ---- GLOBAL heatmap (all numeric) ----
num <- data %>%
  select(where(is.numeric)) %>%
  select(where(~ sd(., na.rm = TRUE) > 0))
M <- cor(num, method = "spearman", use = "pairwise.complete.obs")

png(file.path(out_dir, "correlation_heatmap_spearman.png"),
    width = 2600, height = 2000, res = 300)
corrplot(M, method = "color", type = "upper", diag = FALSE,
         addCoef.col = "black", number.digits = 2, number.cex = 0.5,
         tl.cex = 0.7, order = "hclust")
dev.off()

write_csv(as.data.frame(M), file.path(out_dir, "correlation_matrix_spearman.csv"))

# ---- GE vs. blocks: hormones / glycemia / insulin / demographics ----
run_block <- function(df, left, right, log_vars, prefix, group = NULL) {
  # 1) Pairplot
  save_pairs(file.path(out_dir, paste0("pairs_", prefix, ".png")),
             pairs_plot(df, vars = c(left, right),
                        log10_vars = log_vars,
                        group = group,
                        main = gsub("_", " ", toupper(prefix))))
  # 2) Table
  tab <- spearman_table(df, vars = c(left, right), log10_vars = log_vars, sort_by = "abs_r")
  write_csv(tab, file.path(out_dir, paste0("spearman_", prefix, ".csv")))
  invisible(tab)
}

tab_GE_horm <- run_block(data, vars_outcome, vars_hormones, log_skewed, "GE_hormones")
tab_GE_gly  <- run_block(data, vars_outcome, vars_gly,      log_skewed, "GE_glycemia")
tab_GE_ins  <- run_block(data, vars_outcome, vars_insulin,  log_skewed, "GE_insulin")
tab_GE_demo <- run_block(data, vars_outcome, vars_demo,     log_skewed, "GE_demographics")

# ---- GE-only summary: dotplot of correlations with GE ----
ge_dot <- (function() {
  vars_all <- unique(c(vars_hormones, vars_gly, vars_insulin, vars_demo))
  d <- make_pairs_df(data, keep = c(vars_outcome, vars_all), log10_vars = log_skewed)
  vn <- setdiff(colnames(d), vars_outcome)
  if (!length(vn)) return(invisible(NULL))
  res <- lapply(vn, function(v) {
    x <- d[[vars_outcome]]; y <- d[[v]]
    ok <- is.finite(x) & is.finite(y); n <- sum(ok)
    if (n < 3) return(NULL)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    tibble(var = v, rho = unname(ct$estimate), p = unname(ct$p.value), n = n)
  }) %>% bind_rows() %>%
    mutate(padj = p.adjust(p, method = "BH"),
           block = case_when(
             var %in% vars_hormones ~ "Hormones",
             var %in% vars_gly      ~ "Glycemia",
             var %in% vars_insulin  ~ "Insulin",
             var %in% vars_demo     ~ "Demographics",
             TRUE                   ~ "Other"
           )) %>%
    arrange(desc(abs(rho)))
  write_csv(res, file.path(out_dir, "GE_correlations_summary.csv"))
  
  # dotplot
  res$var <- factor(res$var, levels = res$var)
  p <- ggplot(res, aes(x = var, y = rho, shape = padj < 0.05)) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_point(size = 2.5) +
    coord_flip() +
    facet_grid(. ~ block, scales = "free_x", space = "free_x") +
    scale_shape_manual(values = c(16, 17), labels = c("FDR ≥ 0.05", "FDR < 0.05"), name = NULL) +
    labs(title = "Spearman correlation with GE",
         x = NULL, y = "rho (Spearman)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(out_dir, "GE_correlations_dotplot.png"),
         p, width = 9, height = 6, dpi = 300)
})()


summ <- read_csv(file.path(out_dir, "GE_correlations_summary.csv"), show_col_types = FALSE)
top4 <- summ %>% arrange(desc(abs(rho))) %>% slice(1:4) %>% pull(var)

plot_scatter <- function(df, xvar, yvar = "GE") {
  x <- df[[xvar]]; y <- df[[yvar]]
  ok <- is.finite(x) & is.finite(y)
  r  <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  ggplot(df, aes_string(x = xvar, y = yvar)) +
    geom_point(alpha = 0.35, size = 1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 0.6) +
    labs(x = xvar, y = yvar,
         subtitle = paste0("ρ = ", sprintf("%.2f", r))) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

p <- plot_scatter(data, top4[1]) + plot_scatter(data, top4[2]) +
  plot_scatter(data, top4[3]) + plot_scatter(data, top4[4])
ggsave(file.path(out_dir, "GE_top4_scatter.png"), p, width = 10, height = 8, dpi = 300)


# ============================================================
# GE modelling!!!
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(broom)
library(scales)
library(sandwich)
library(lmtest)
library(car)    

HAS_GLMNET   <- requireNamespace("glmnet",   quietly = TRUE)
HAS_PATCHWORK<- requireNamespace("patchwork",quietly = TRUE)
HAS_VISREG   <- requireNamespace("visreg",   quietly = TRUE)


# ------------------ paths ------------------
data_path <- "/Users/claraazzano/Desktop/Data_T1.csv"
out_dir   <- "/Users/claraazzano/task1_models?"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------ load/transform ------------------
dat0 <- read_csv(data_path, show_col_types = FALSE)

# pick ONE rep for each correlated block
gly_marker <- "HbA1c"          # or "GlucoseFasting"
ins_marker <- "InsulinFasting" # or "MatsudaIdx"

# log-transform skewed hormones (natural log)
dat <- dat0 %>%
  mutate(
    GLP1_ln    = ifelse(GLP1    > 0, log(GLP1),    NA_real_),
    Ghrelin_ln = ifelse(Ghrelin > 0, log(Ghrelin), NA_real_),
    Amylin_ln  = ifelse(Amylin  > 0, log(Amylin),  NA_real_),
    # treat binaries as factors if present
    Sex               = if ("Sex"               %in% names(.)) factor(Sex)               else NULL,
    Metformin         = if ("Metformin"         %in% names(.)) factor(Metformin)         else NULL,
    DiabetesComplications = if ("DiabetesComplications" %in% names(.)) factor(DiabetesComplications) else NULL
  )

base_covars <- c("Age","Sex","BMI","Metformin")               # forced-in covariates
extra_core  <- c(gly_marker, ins_marker)
hormones_ln <- c("GLP1_ln","Ghrelin_ln","Amylin_ln","Gastrin","CCK","Glucagon","PYY")

all_vars <- c("GE", base_covars, extra_core, hormones_ln)
df <- dat %>% select(any_of(all_vars)) %>% tidyr::drop_na()

# ------------------ helpers ------------------
robust_table <- function(fit) {
  ct <- coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
  tb <- data.frame(term = rownames(ct), estimate = ct[,1], std.error = ct[,2],
                   statistic = ct[,3], p.value = ct[,4], row.names = NULL)
  tb$conf.low  <- tb$estimate - 1.96*tb$std.error
  tb$conf.high <- tb$estimate + 1.96*tb$std.error
  tb
}

# ------------------ models: M0/M1 ------------------
form_M0 <- as.formula(paste("GE ~", paste(base_covars, collapse = " + ")))
form_M1 <- as.formula(paste("GE ~", paste(c(base_covars, extra_core, hormones_ln), collapse = " + ")))

M0 <- lm(form_M0, data = df)
M1 <- lm(form_M1, data = df)

write_csv(glance(M0), file.path(out_dir, "M0_baseline_glance.csv"))
write_csv(glance(M1), file.path(out_dir, "M1_extended_glance.csv"))
VIF_M1 <- tryCatch(car::vif(M1), error = function(e) NA)
if (!all(is.na(VIF_M1))) write_csv(data.frame(term=names(VIF_M1), VIF=as.numeric(VIF_M1)),
                                   file.path(out_dir, "M1_VIF.csv"))
write_csv(robust_table(M0), file.path(out_dir, "M0_coef_robust.csv"))
write_csv(robust_table(M1), file.path(out_dir, "M1_coef_robust.csv"))

# ------------------ selection: LASSO (or stepwise fallback) ------------------

if (HAS_GLMNET) {
  library(glmnet)
  
  # Build model.matrix and use its 'assign' attribute to map columns -> terms
  all_preds <- attr(terms(M1), "term.labels")
  Terms     <- terms(reformulate(all_preds))          # terms object
  MM        <- model.matrix(Terms, data = df)         # includes intercept
  assignVec <- attr(MM, "assign")                     # column -> term index
  termLabs  <- attr(Terms, "term.labels")
  
  # Drop intercept column
  X       <- MM[, -1, drop = FALSE]
  assignX <- assignVec[-1]
  
  # Unpenalize columns that belong to base covariates (handles factor dummies)
  pf <- rep(1, ncol(X)); names(pf) <- colnames(X)
  base_idx <- which(termLabs %in% base_covars)
  pf[assignX %in% base_idx] <- 0
  
  # CV LASSO
  set.seed(123)
  cvfit <- cv.glmnet(x = X, y = df$GE, alpha = 1, family = "gaussian",
                     nfolds = 10, penalty.factor = pf)
  
  png(file.path(out_dir, "lasso_cv_curve.png"), width = 1200, height = 900, res = 180)
  plot(cvfit); dev.off()
  
  # 1-SE rule, map selected columns back to their originating term names
  lam1se   <- cvfit$lambda.1se
  beta     <- as.numeric(coef(cvfit, s = lam1se))[-1]     # drop intercept
  keep_cols <- which(beta != 0)
  sel_terms <- unique(termLabs[ assignX[keep_cols] ])
  write_csv(data.frame(selected = sel_terms),
            file.path(out_dir, "lasso_selected_lambda1se.csv"))
  
  # Final model = base covariates + selected terms
  final_terms <- unique(c(base_covars, sel_terms))
  form_final  <- reformulate(final_terms, response = "GE")
  Mfinal      <- lm(form_final, data = df)
  
} else {
  # Fallback: stepwise with base covariates forced in
  library(MASS)
  scope  <- list(lower = reformulate(base_covars, response = "GE"),
                 upper  = form_M1)
  Mfinal <- stepAIC(M1, scope = scope, direction = "both", trace = FALSE)
}

# save final model summaries
write_csv(glance(Mfinal),              file.path(out_dir, "Mfinal_glance.csv"))
write_csv(robust_table(Mfinal),        file.path(out_dir, "Mfinal_coef_robust.csv"))

# diagnostics
png(file.path(out_dir, "Mfinal_diagnostics.png"), width=1600, height=1600, res=180)
par(mfrow=c(2,2)); plot(Mfinal); dev.off()

# ============================================================
# VISUALISATIONS (integrated)
# ============================================================

# 1) Coefficient forest (robust CIs)
coef_final <- robust_table(Mfinal) %>% 
  filter(term != "(Intercept)") %>%
  mutate(term = forcats::fct_reorder(term, estimate))
p_coef <- ggplot(coef_final, aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.3) +
  coord_flip() +
  labs(title = "Final model coefficients (HC3 robust 95% CI)",
       x = NULL, y = "β") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "Mfinal_coef_forest.png"), p_coef, width = 8, height = 6, dpi = 300)

# helpers for CV predictions and metrics
cv_predict_lm <- function(formula, data, k = 10, seed = 123) {
  set.seed(seed)
  n <- nrow(data); fold <- sample(rep(1:k, length.out = n))
  pred <- rep(NA_real_, n)
  for (i in seq_len(k)) {
    fit <- lm(formula, data = data[fold != i, , drop = FALSE])
    pred[fold == i] <- predict(fit, newdata = data[fold == i, , drop = FALSE])
  }
  pred
}
perf_metrics <- function(y, yhat) {
  rmse <- sqrt(mean((y - yhat)^2))
  r2   <- 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
  list(rmse = rmse, r2 = r2)
}
plot_obs_pred <- function(obs, pred, title, subtitle) {
  m <- perf_metrics(obs, pred)
  ggplot(data.frame(obs, pred), aes(pred, obs)) +
    geom_abline(linetype = 2) +
    geom_point(alpha = .45, size = 1) +
    labs(title = title,
         subtitle = sprintf("%s  •  RMSE=%.2f,  R²=%.3f", subtitle, m$rmse, m$r2),
         x = "Predicted", y = "Observed GE") +
    theme_minimal(base_size = 12)
}

# 2) Observed vs Predicted (in-sample & 10-fold CV) for M0 & Mfinal
df$.pred_M0_in     <- fitted(M0)
df$.pred_M0_cv     <- cv_predict_lm(form_M0, df)
df$.pred_Mfinal_in <- fitted(Mfinal)
df$.pred_Mfinal_cv <- cv_predict_lm(formula(Mfinal), df)

p1 <- plot_obs_pred(df$GE, df$.pred_M0_in,     "Baseline (M0)",        "In-sample")
p2 <- plot_obs_pred(df$GE, df$.pred_M0_cv,     "Baseline (M0)",        "10-fold CV")
p3 <- plot_obs_pred(df$GE, df$.pred_Mfinal_in, "Final model (Mfinal)", "In-sample")
p4 <- plot_obs_pred(df$GE, df$.pred_Mfinal_cv, "Final model (Mfinal)", "10-fold CV")

if (HAS_PATCHWORK) {
  library(patchwork)
  g <- (p1 | p2) / (p3 | p4)
  ggsave(file.path(out_dir, "obs_vs_pred_M0_Mfinal.png"), g, width = 10, height = 8, dpi = 300)
} else {
  ggsave(file.path(out_dir, "obs_pred_M0_in.png"),     p1, width=6, height=5, dpi=300)
  ggsave(file.path(out_dir, "obs_pred_M0_cv.png"),     p2, width=6, height=5, dpi=300)
  ggsave(file.path(out_dir, "obs_pred_Mfinal_in.png"), p3, width=6, height=5, dpi=300)
  ggsave(file.path(out_dir, "obs_pred_Mfinal_cv.png"), p4, width=6, height=5, dpi=300)
}

# 3) Calibration (deciles) for Mfinal (CV preds)
cal_df <- tibble(obs = df$GE, pred = df$.pred_Mfinal_cv) |> tidyr::drop_na() |>
  mutate(bin = ntile(pred, 10)) |>
  group_by(bin) |> summarise(pred = mean(pred), obs = mean(obs), .groups = "drop")
p_cal <- ggplot(cal_df, aes(pred, obs)) +
  geom_abline(linetype = 2) +
  geom_point(size = 2) + geom_line() +
  labs(title = "Calibration of final model (10-fold CV)",
       x = "Mean predicted GE (by decile)", y = "Mean observed GE") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "Mfinal_calibration.png"), p_cal, width = 7, height = 5, dpi = 300)

# 4) Partial effects for most influential numeric predictors
# choose top 4 by standardized beta (numeric terms only)
std_df <- df
num_vars <- names(std_df)[sapply(std_df, is.numeric)]
num_vars <- setdiff(num_vars, c("GE"))
std_df[num_vars] <- lapply(std_df[num_vars], scale)

M_std <- lm(formula(Mfinal), data = std_df)
std_tab <- tidy(M_std) %>% 
  filter(term != "(Intercept)") %>%
  arrange(desc(abs(estimate)))
top4 <- head(std_tab$term[!grepl("^Sex|^Metformin", std_tab$term)], 4)

if (HAS_VISREG) {
  library(visreg)
  ps <- lapply(top4, function(v) {
    visreg(Mfinal, v, gg = TRUE, overlay = FALSE, rug = TRUE) +
      labs(title = paste("Adjusted effect:", v), x = v, y = "GE (partial)") +
      theme_minimal(base_size = 12)
  })
  if (HAS_PATCHWORK) {
    library(patchwork)
    pg <- wrap_plots(ps, ncol = 2)
    ggsave(file.path(out_dir, "Mfinal_partial_effects.png"), pg, width = 10, height = 8, dpi = 300)
  } else {
    # save individually
    for (i in seq_along(ps)) {
      ggsave(file.path(out_dir, paste0("partial_", gsub("[^A-Za-z0-9]+","_", top4[i]), ".png")),
             ps[[i]], width=6, height=5, dpi=300)
    }
  }
} else {
  # fallback: added-variable plots
  png(file.path(out_dir, "Mfinal_added_variable_top4.png"), width = 1600, height = 1200, res = 180)
  par(mfrow = c(2,2))
  car::avPlots(Mfinal, terms = top4)
  dev.off()
}

message("Done. Outputs in: ", normalizePath(out_dir, mustWork = FALSE))

## ===============================
## Reporting add-ons for the paper
## ===============================

library(dplyr)
library(ggplot2)
library(broom)
library(readr)
library(sandwich)
library(lmtest)
theme_set(theme_minimal(base_size = 13))

#--- where to save
out_dir <- "model_report"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "added_variable_plots"), showWarnings = FALSE)

#--- helper: RMSE and simple K-fold CV RMSE for lm()
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))

cv_rmse_lm <- function(formula, data, K = 10, seed = 1) {
  set.seed(seed)
  n <- nrow(data); idx <- sample.int(n)
  folds <- split(idx, cut(seq_along(idx), K, labels = FALSE))
  errs <- numeric(length(folds))
  for (i in seq_along(folds)) {
    test_id  <- folds[[i]]
    train_id <- setdiff(idx, test_id)
    fit <- lm(formula, data = data[train_id, , drop = FALSE])
    pred <- predict(fit, newdata = data[test_id, , drop = FALSE])
    errs[i] <- rmse(data$GE[test_id], pred)
  }
  mean(errs)
}

#--- 1) Model comparison: M0 vs Mfinal
stopifnot(exists("M0"), exists("Mfinal"))
df_final <- model.frame(Mfinal)   # exact data used by final model

tab_models <- tibble(
  model     = c("M0_baseline", "Mfinal"),
  adj_r2    = c(glance(M0)$adj.r.squared, glance(Mfinal)$adj.r.squared),
  AIC       = c(AIC(M0), AIC(Mfinal)),
  BIC       = c(BIC(M0), BIC(Mfinal)),
  RMSE_in   = c(rmse(df_final$GE, fitted(M0)), rmse(df_final$GE, fitted(Mfinal))),
  RMSE_cv10 = c(cv_rmse_lm(formula(M0), df_final, K = 10, seed = 42),
                cv_rmse_lm(formula(Mfinal), df_final, K = 10, seed = 42))
)

print(tab_models)
write_csv(tab_models, file.path(out_dir, "M0_vs_Mfinal_metrics.csv"))

#--- 2) Added-variable (partial regression) plots for selected terms
avp <- function(fit, term, data, title_prefix = NULL) {
  all_terms <- attr(terms(fit), "term.labels")
  stopifnot(term %in% all_terms)
  others <- setdiff(all_terms, term)
  
  # residuals of GE ~ others
  f_y <- reformulate(others, response = all.vars(formula(fit))[1])
  ry  <- resid(lm(f_y, data = data))
  
  # residuals of term ~ others (term must be on LHS of reformulate)
  f_x <- reformulate(others, response = term)
  rx  <- resid(lm(f_x, data = data))
  
  d <- tibble(rx = rx, ry = ry)
  ggplot(d, aes(rx, ry)) +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_vline(xintercept = 0, linewidth = 0.2) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = paste(term, "| others"),
         y = paste(all.vars(formula(fit))[1], "| others"),
         title = paste0(title_prefix %||% "Added-variable: ", term))
}

sel_terms <- attr(terms(Mfinal), "term.labels")
# focus on hormones that survived (any *_ln or hormone names present)
candidate_horms <- c("GLP1_ln","Ghrelin_ln","Amylin_ln","Gastrin","CCK","PYY","Glucagon")
av_terms <- intersect(sel_terms, candidate_horms)

for (trm in av_terms) {
  p <- avp(Mfinal, trm, df_final, title_prefix = "Added-variable: ")
  ggsave(file.path(out_dir, "added_variable_plots", paste0("AVP_", trm, ".png")),
         p, width = 7, height = 4.5, dpi = 300)
}

#--- 3) Robust coefficient table (HC3) for Mfinal
rob <- coeftest(Mfinal, vcov = sandwich::vcovHC(Mfinal, type = "HC3"))
tbl_coef <- tibble(
  term = rownames(rob),
  estimate = rob[,1],
  std.error = rob[,2],
  statistic = rob[,3],
  p.value = rob[,4]
) %>%
  mutate(
    conf.low  = estimate - 1.96*std.error,
    conf.high = estimate + 1.96*std.error
  )

print(tbl_coef)
write_csv(tbl_coef, file.path(out_dir, "Mfinal_coef_robust.csv"))

#--- 4) Plain-language effect-size translations

LN_10pct <- log(1.1)   # = 0.09531
LN_2x    <- log(2)     # = 0.69315

describe_effect <- function(term, beta, data) {
  out <- tibble(term = term,
                beta = beta,
                delta_10pct = NA_real_,
                delta_2x    = NA_real_,
                delta_per1  = NA_real_,
                delta_perIQR = NA_real_,
                note = NA_character_)
  if (grepl("_ln$", term)) {
    out$delta_10pct <- beta * LN_10pct
    out$delta_2x    <- beta * LN_2x
    out$note <- "natural log: ΔGE for 10% and 2× increase"
    return(out)
  }
  if (grepl("^Sex", term)) {
    out$note <- "difference vs Sex reference level"
    return(out)
  }
  if (grepl("^Metformin", term)) {
    out$note <- "difference vs Metformin reference"
    return(out)
  }
  v <- term
  if (v %in% names(data) && is.numeric(data[[v]])) {
    out$delta_per1   <- beta
    out$delta_perIQR <- beta * IQR(data[[v]], na.rm = TRUE)
    out$note <- "continuous: per-unit and per-IQR change"
  } else {
    out$note <- "contrast or transformed term"
  }
  out
}

terms_in_model <- setdiff(tbl_coef$term, "(Intercept)")
fx_table <- purrr::map_dfr(
  terms_in_model,
  ~describe_effect(.x, beta = tbl_coef$estimate[tbl_coef$term == .x], data = df_final)
) %>%
  left_join(tbl_coef %>% select(term, conf.low, conf.high, p.value), by = "term")

write_csv(fx_table, file.path(out_dir, "Mfinal_effect_size_translations.csv"))

# Optional quick dot plot
plot_df <- fx_table %>%
  filter(!is.na(delta_perIQR) | !is.na(delta_2x)) %>%
  mutate(group = dplyr::case_when(
    grepl("_ln$", term)            ~ "Hormones (ln)",
    grepl("^Sex|^Metformin", term) ~ "Binary contrasts",
    TRUE                           ~ "Continuous"
  ))

p_dot <- ggplot(plot_df,
                aes(x = ifelse(is.na(delta_perIQR), delta_2x, delta_perIQR),
                    y = term)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Effect on GE (per-IQR or 2× / 10%)", y = NULL,
       title = "Effect-size translations (final model)")
ggsave(file.path(out_dir, "Mfinal_effect_sizes_dotplot.png"),
       p_dot, width = 9, height = 6, dpi = 300)

#--- 5) (optional) quick dot plot of per-IQR effects for continuous terms
plot_df <- fx_table %>%
  filter(!is.na(delta_perIQR) | !is.na(delta_2x)) %>%
  mutate(group = dplyr::case_when(
    grepl("_ln$", term)              ~ "Hormones (log10)",
    grepl("^Sex|^Metformin", term)   ~ "Binary contrasts",
    TRUE                             ~ "Continuous"
  ))

p_dot <- ggplot(plot_df,
                aes(x = ifelse(is.na(delta_perIQR), delta_2x, delta_perIQR),
                    y = term)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Effect on GE (per-IQR or 2× / 10% as available)", y = NULL,
       title = "Effect-size translations (final model)") +
  facet_wrap(~group, scales = "free_x")
ggsave(file.path(out_dir, "Mfinal_effect_sizes_dotplot.png"),
       p_dot, width = 9, height = 6, dpi = 300)

cat("\nWrote:",
    "\n  -", file.path(out_dir, "M0_vs_Mfinal_metrics.csv"),
    "\n  -", file.path(out_dir, "Mfinal_coef_robust.csv"),
    "\n  -", file.path(out_dir, "Mfinal_effect_size_translations.csv"),
    "\n  - added-variable plots in:", file.path(out_dir, "added_variable_plots"),
    "\n  -", file.path(out_dir, "Mfinal_effect_sizes_dotplot.png"), "\n")

