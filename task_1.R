# ================================================================
# TASK 1 - VISUALISATIONS
# ================================================================
# Goal of this part:
#   1) Show the empirical distribution of all continuous variables with histograms.
#   2) Summarise binary variables as clean "n (%)" tables (Console print).
#   3) Provide descriptive stats and a missingness table (Console print).
# Design choices:
#   - Light auto-detection of skew → suggest log10 scale for certain variables.
#   - All plots print to the RStudio Plots pane.
#   - Nothing is written to disk unless SAVE_OUTPUTS <- TRUE.

# Toggle: set TRUE only if you also want to save PNG/CSV outputs.
SAVE_OUTPUTS <- FALSE

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(scales)
library(rlang)
theme_set(theme_minimal(base_size = 12))


data <- read_csv("Data_T1.csv", show_col_types = FALSE)


# Helpers ==========================================================================
# Purpose: small utilities used by multiple plots/tables, documented inline.

# Compute Freedman–Diaconis bin width (fallback to 30 bins if needed).
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

# Quick skew flag via 3rd standardized moment (|skew| > thresh → TRUE).
is_skewed <- function(x, thresh = 1) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || sd(x) == 0) return(FALSE)
  z  <- (x - mean(x)) / sd(x)
  m3 <- mean(z^3)
  abs(m3) > thresh
}

# Log10 transform that safely handles non-positive values by adding a small offset.
# The chosen offset is stored as an attribute so we can annotate the plot.
safe_log10 <- function(x, eps = 1e-6) {
  x <- as.numeric(x)
  needs_offset <- any(x <= 0, na.rm = TRUE)
  off <- if (needs_offset) max(eps, abs(min(x, na.rm = TRUE)) + eps) else 0
  y <- log10(x + off)
  attr(y, "offset") <- off
  y
}

# Build an x-axis label with optional units and a "(log10)" suffix when used.
axis_label <- function(name, log = FALSE, units = NULL) {
  paste0(
    name,
    if (!is.null(units) && nzchar(units)) paste0(" [", units, "]") else "",
    if (log) " (log10)" else ""
  )
}


# Variables and units ==============================================================
# What we plot/summarise and which units to show in axis labels.

cont_vars <- c(
  "GE","Age","Height","BW","BMI",
  "GlucoseFasting","InsulinFasting","HbA1c","MatsudaIdx","HOMAB",
  "Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY"
)

# Variables that are typically very skewed → default to log scale in histograms.
log_suggest <- c("InsulinFasting","MatsudaIdx","Ghrelin","Amylin","GLP1")

# 0/1 indicators we summarise as counts & percentages.
binary_vars <- c("Sex","Metformin","DiabetesComplications")

# Units used in histogram x-axis labels (only cosmetic; analysis uses raw units).
units_map <- c(
  GE = "min", Age = "y", Height = "cm", BW = "kg", BMI = "kg/m^2",
  GlucoseFasting = "mmol/L", InsulinFasting = "mU/L", HbA1c = "%",
  Gastrin = "pg/mL", CCK = "pmol/L", Ghrelin = "pg/mL", Amylin = "pmol/L",
  Glucagon = "pmol/L", GLP1 = "pmol/L", PYY = "ng/L",
  MatsudaIdx = "", HOMAB = ""
)

# Safe lookup for units (avoids indexing errors if a key is missing).
get_units <- function(var) {
  nm <- as.character(var)
  if (!is.null(names(units_map)) && nm %in% names(units_map)) units_map[[nm]] else NULL
}


# Plot histograms ==================================================================
# Single function that decides on log scale (if requested/auto), computes FD bin width,
# prints n / NA counts, and annotates the added offset if log10(x + offset) was used.

plot_hist <- function(df, var, bins = NULL, log_x = NULL,
                      units = NULL, title = NULL, annotate_n = TRUE) {
  v <- enquo(var)
  vname <- as_name(v)
  x <- df %>% pull(!!v)
  
  # If log_x not specified, auto-enable when data are positive and skewed.
  if (is.null(log_x)) log_x <- all(x > 0, na.rm = TRUE) && is_skewed(x)
  
  # Compute plotting vector (possibly with safe log10 + offset).
  x_plot <- if (log_x) safe_log10(x) else x
  off <- if (log_x) attr(x_plot, "offset") else 0
  df_plot <- tibble(x_plot = as.numeric(x_plot))
  
  # FD binning by default; explicit bins overrides it.
  if (is.null(bins)) {
    bw <- fd_binwidth(df_plot$x_plot)
    p <- ggplot(df_plot, aes(x = x_plot)) +
      geom_histogram(binwidth = bw, fill = "#6A9AE2", color = "white", na.rm = TRUE)
  } else {
    p <- ggplot(df_plot, aes(x = x_plot)) +
      geom_histogram(bins = bins, fill = "#6A9AE2", color = "white", na.rm = TRUE)
  }
  
  # Labels + lighter minor grid.
  p <- p +
    labs(
      title = if (is.null(title)) vname else title,
      x = axis_label(vname, log = isTRUE(log_x), units = if (is.null(units)) get_units(vname) else units),
      y = "Count"
    ) +
    theme(panel.grid.minor = element_blank())
  
  # Optional annotation with n and NA; add offset info if safe_log10 was used.
  if (annotate_n) {
    n_use <- sum(is.finite(df_plot$x_plot))
    miss  <- sum(!is.finite(df_plot$x_plot))
    ann <- paste0("n=", n_use, if (miss > 0) paste0(" (NA=", miss, ")"))
    if (isTRUE(log_x) && off > 0) ann <- paste0(ann, "\nlog10(x + ", signif(off, 3), ")")
    p <- p + annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5, label = ann)
  }
  
  p
}

# ---- print each histogram to Plots pane ------------------------------------------
invisible(lapply(cont_vars, function(v) {
  print(
    plot_hist(
      data, !!sym(v),
      log_x = if (v %in% log_suggest) TRUE else NULL,
      units = get_units(v)
    )
  )
}))

# ---- optionally save the same histograms to ./task1_histograms -------------------
if (SAVE_OUTPUTS) {
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
}


# Binary summary (prints tidy tables to Console) ===================================
# Clean "n (%)" by level for Sex, Metformin, and DiabetesComplications.

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

# Human-readable labels for 0/1.
label_map <- list(
  Sex = c("Female","Male"),
  Metformin = c("No","Yes"),
  DiabetesComplications = c("No","Yes")
)

# Print long and wide tables to Console; optionally save wide.
bin_tbl <- binary_summary(data, binary_vars, labels = label_map)

cat("\n--- Binary variables: counts and percentages ---\n")
bin_tbl %>%
  arrange(Variable, Level) %>%
  mutate(`n (%)` = paste0(n, " (", pct_label, ")")) %>%
  select(Variable, Label, `n (%)`) %>%
  print(n = Inf)

# Wide "quick glance" view (columns per level).
bin_tbl_wide <- bin_tbl %>%
  mutate(`n (%)` = paste0(n, " (", pct_label, ")")) %>%
  select(Variable, Label, `n (%)`) %>%
  pivot_wider(names_from = Label, values_from = `n (%)`)

cat("\n--- Binary variables (wide view) ---\n")
print(bin_tbl_wide, n = Inf)

if (SAVE_OUTPUTS) write_csv(bin_tbl_wide, "binary_summary_task1.csv")


# Descriptive stats for continuous + missingness (Console prints) ==================
# Summary: n, mean, sd, median, IQR, min, max for each continuous variable.
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

cat("\n--- Descriptive stats (continuous) ---\n")
print(desc, n = Inf)

if (SAVE_OUTPUTS) write_csv(desc, "summary_task1.csv")

# Per-variable NA counts (descending).
missing_tbl <- data %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to="variable", values_to="n_missing") %>%
  arrange(desc(n_missing))

cat("\n--- Missingness by variable ---\n")
print(missing_tbl, n = Inf)



# ================================================================
# TASK 1 - CORRELATIONS
# ================================================================
# Goal of this part:
#   1) Global Spearman correlation heatmap (numeric variables only).
#   2) Pairwise plots + tables for GE vs. each variable block (hormones, glycemia, insulin, demographics).
#   3) A dotplot that ranks correlations with GE across all blocks.
# Printing only; optional saving guarded by SAVE_OUTPUTS.

SAVE_OUTPUTS <- FALSE   # <- set TRUE to also save PNG/CSV
HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(corrplot)
library(rlang)
theme_set(theme_minimal(base_size = 12))
if (HAS_PATCHWORK) library(patchwork)

data <- read_csv("Data_T1.csv", show_col_types = FALSE)


#  Variable sets ================================================================
# Outcome (GE) and the blocks of predictors we will correlate against it.
vars_hormones <- c("Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY")
vars_demo     <- c("Age","Height","BW","BMI")
vars_gly      <- c("GlucoseFasting","HbA1c")
vars_insulin  <- c("InsulinFasting","MatsudaIdx","HOMAB")
vars_outcome  <- "GE"

# Variables we’ll log10 for correlation EDA because of skew.
log_skewed    <- c("InsulinFasting","MatsudaIdx","Ghrelin","Amylin","GLP1")


# Helpers ========================================================================
# safe_log10: same idea as before—adds an offset if non-positive values present.
safe_log10 <- function(x, eps = 1e-6) {
  x <- as.numeric(x)
  needs_offset <- any(x <= 0, na.rm = TRUE)
  off <- if (needs_offset) max(eps, abs(min(x, na.rm = TRUE)) + eps) else 0
  y <- log10(x + off)
  attr(y, "offset") <- off
  y
}

# make_pairs_df: subset + numeric-only + drop zero-variance columns + optional log10.
make_pairs_df <- function(df, keep, log10_vars = character()) {
  keep <- intersect(keep, names(df))
  d <- df[keep]
  d <- d[vapply(d, is.numeric, logical(1L))]                 # numeric only
  for (v in intersect(log10_vars, names(d))) d[[v]] <- safe_log10(d[[v]])
  if (ncol(d)) {
    ok <- vapply(d, function(z) sd(z, na.rm = TRUE) > 0, logical(1L))
    d <- d[, ok, drop = FALSE]
  }
  d
}


# Spearman correlations (pairwise), BH-adjusted ===================================
# Returns a sorted tibble of (var1, var2, rho, p, n, padj).
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


# Base-pairs ======================================================================
# A compact base::pairs() with:
#   - lower: point cloud with alpha,
#   - upper: shaded background by |rho| (Spearman),
#   - diag:  kernel density (normalized to [0,1.5] y-range for readability).
pairs_plot <- function(df, vars, log10_vars = character(),
                       group = NULL, main = "") {
  d <- make_pairs_df(df, keep = vars, log10_vars = log10_vars)
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
  invisible(NULL)
}


# Optional saver (used only if SAVE_OUTPUTS) --------------------------------------
save_pairs <- function(filename, expr, width = 3200, height = 2400, res = 300) {
  png(filename, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}


# GLOBAL heatmap ===================================================================
# Spearman correlation across all numeric variables (pairwise NAs).
num <- data %>%
  select(where(is.numeric)) %>%
  select(where(~ sd(., na.rm = TRUE) > 0))
M <- cor(num, method = "spearman", use = "pairwise.complete.obs")

cat("\n--- Global Spearman correlation heatmap (printed) ---\n")
corrplot(M, method = "color", type = "upper", diag = FALSE,
         addCoef.col = "black", number.digits = 2, number.cex = 0.5,
         tl.cex = 0.7, order = "hclust")

if (SAVE_OUTPUTS) {
  out_dir <- "correlations_task1"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  png(file.path(out_dir, "correlation_heatmap_spearman.png"),
      width = 2600, height = 2000, res = 300)
  corrplot(M, method = "color", type = "upper", diag = FALSE,
           addCoef.col = "black", number.digits = 2, number.cex = 0.5,
           tl.cex = 0.7, order = "hclust")
  dev.off()
  write_csv(as.data.frame(M), file.path(out_dir, "correlation_matrix_spearman.csv"))
}


# GE vs. blocks: hormones / glycemia / insulin / demographics =====================
# For each block:
#   1) Print the pairplot (GE + block variables).
#   2) Print the top correlations sorted by |rho|.
#   3) Optionally save the image and the full table.

run_block <- function(df, left, right, log_vars, prefix, group = NULL, save = FALSE) {
  title <- gsub("_", " ", toupper(prefix))
  
  # 1) Pairplot (printed)
  cat("\n--- Pairplot:", title, "---\n")
  pairs_plot(df, vars = c(left, right), log10_vars = log_vars, group = group, main = title)
  
  # 2) Table (printed)
  tab <- spearman_table(df, vars = c(left, right), log10_vars = log_vars, sort_by = "abs_r")
  cat("\nTop correlations (|rho|):", title, "\n")
  print(head(tab, 20))
  
  # 3) Optional save
  if (save) {
    out_dir <- "correlations_task1"
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    save_pairs(file.path(out_dir, paste0("pairs_", prefix, ".png")),
               pairs_plot(df, vars = c(left, right), log10_vars = log_vars, group = group, main = title))
    write_csv(tab, file.path(out_dir, paste0("spearman_", prefix, ".csv")))
  }
  
  invisible(tab)
}

tab_GE_horm <- run_block(data, vars_outcome, vars_hormones, log_skewed, "GE_hormones", save = SAVE_OUTPUTS)
tab_GE_gly  <- run_block(data, vars_outcome, vars_gly,      log_skewed, "GE_glycemia", save = SAVE_OUTPUTS)
tab_GE_ins  <- run_block(data, vars_outcome, vars_insulin,  log_skewed, "GE_insulin",  save = SAVE_OUTPUTS)
tab_GE_demo <- run_block(data, vars_outcome, vars_demo,     log_skewed, "GE_demographics", save = SAVE_OUTPUTS)


# GE-only summary: dotplot of correlations with GE =================================
# Ranks all variables by |rho| (Spearman), BH-adjusts p, and facets by block.
cat("\n--- GE correlations summary (printed) ---\n")
vars_all <- unique(c(vars_hormones, vars_gly, vars_insulin, vars_demo))
d_ge <- make_pairs_df(data, keep = c(vars_outcome, vars_all), log10_vars = log_skewed)
vn <- setdiff(colnames(d_ge), vars_outcome)

res <- lapply(vn, function(v) {
  x <- d_ge[[vars_outcome]]; y <- d_ge[[v]]
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

print(head(res, 20))  # show top 20 in Console

# Dotplot (printed)
res$var <- factor(res$var, levels = res$var)
p_dot <- ggplot(res, aes(x = var, y = rho, shape = padj < 0.05)) +
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
print(p_dot)

if (SAVE_OUTPUTS) {
  out_dir <- "correlations_task1"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out_dir, "GE_correlations_dotplot.png"), p_dot, width = 9, height = 6, dpi = 300)
  write_csv(res, file.path(out_dir, "GE_correlations_summary.csv"))
}


# Quick scatter checks for top-4 absolute correlations with GE =====================
# Convenience check: LOESS trend and annotated Spearman rho.
top4 <- res %>% arrange(desc(abs(rho))) %>% slice(1:4) %>% pull(var)

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

if (length(top4) >= 4) {
  p1 <- plot_scatter(data, top4[1])
  p2 <- plot_scatter(data, top4[2])
  p3 <- plot_scatter(data, top4[3])
  p4 <- plot_scatter(data, top4[4])
  
  if (HAS_PATCHWORK) {
    print((p1 | p2) / (p3 | p4))
    if (SAVE_OUTPUTS) {
      out_dir <- "correlations_task1"
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      ggsave(file.path(out_dir, "GE_top4_scatter.png"),
             (p1 | p2) / (p3 | p4), width = 10, height = 8, dpi = 300)
    }
  } else {
    print(p1); print(p2); print(p3); print(p4)
  }
}



# ============================================================
# TASK 1 - GE modelling
# ============================================================
# Goal of this part:
#   1) Fit baseline model M0 (Age, Sex, BMI, Metformin).
#   2) Fit extended model M1 (M0 + glycemia + insulin + hormones).
#   3) Use LASSO (keeping baseline unpenalized) to select a sparse final model Mfinal.
#   4) Print diagnostics, coefficient forest, obs-vs-pred, calibration, and partial effects.
#   5) Print model comparison metrics (adj R², AIC/BIC, RMSE, CV-RMSE).
# Printing only; optional saving via SAVE_OUTPUTS.

SAVE_OUTPUTS <- FALSE  # <- set TRUE to also save PNG/CSV
HAS_GLMNET    <- requireNamespace("glmnet",    quietly = TRUE)
HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)
HAS_VISREG    <- requireNamespace("visreg",    quietly = TRUE)

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(broom)
library(scales)
library(sandwich)
library(lmtest)
library(car)
theme_set(theme_minimal(base_size = 12))

data <- read_csv("Data_T1.csv", show_col_types = FALSE)

# pick ONE representative for each correlated block (avoid multicollinearity).
gly_marker <- "HbA1c"          # or "GlucoseFasting"
ins_marker <- "InsulinFasting" # or "MatsudaIdx"

# ------------------ transform -----------------------------------------------------
# Natural-log transform selected hormones to linearize their relationships.
# Also coerce binaries to factors so contrasts are explicit in the model.
dat <- dat0 %>%
  mutate(
    GLP1_ln    = ifelse(GLP1    > 0, log(GLP1),    NA_real_),
    Ghrelin_ln = ifelse(Ghrelin > 0, log(Ghrelin), NA_real_),
    Amylin_ln  = ifelse(Amylin  > 0, log(Amylin),  NA_real_),
    Sex                   = if ("Sex" %in% names(.)) factor(Sex) else NULL,
    Metformin             = if ("Metformin" %in% names(.)) factor(Metformin) else NULL,
    DiabetesComplications = if ("DiabetesComplications" %in% names(.)) factor(DiabetesComplications) else NULL
  )

# Variables included in models.
base_covars <- c("Age","Sex","BMI","Metformin")               # always adjust for these
extra_core  <- c(gly_marker, ins_marker)
hormones_ln <- c("GLP1_ln","Ghrelin_ln","Amylin_ln","Gastrin","CCK","Glucagon","PYY")

# Assemble complete modelling frame and drop rows with any NA in selected variables.
all_vars <- c("GE", base_covars, extra_core, hormones_ln)
df <- dat %>% select(any_of(all_vars)) %>% tidyr::drop_na()

# ------------------ helpers -------------------------------------------------------
# HC3-robust coefficient table (β, SE, z, p, 95% CI).
robust_table <- function(fit) {
  ct <- coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
  tb <- data.frame(term = rownames(ct), estimate = ct[,1], std.error = ct[,2],
                   statistic = ct[,3], p.value = ct[,4], row.names = NULL)
  tb$conf.low  <- tb$estimate - 1.96*tb$std.error
  tb$conf.high <- tb$estimate + 1.96*tb$std.error
  tb
}

# ------------------ models: M0 (baseline) & M1 (extended) ------------------------
form_M0 <- as.formula(paste("GE ~", paste(base_covars, collapse = " + ")))
form_M1 <- as.formula(paste("GE ~", paste(c(base_covars, extra_core, hormones_ln), collapse = " + ")))

M0 <- lm(form_M0, data = df)
M1 <- lm(form_M1, data = df)

cat("\n=== Baseline model (M0) ===\n"); print(glance(M0)); cat("\n"); print(robust_table(M0))
cat("\n=== Extended model (M1) ===\n"); print(glance(M1)); cat("\n")
VIF_M1 <- tryCatch(car::vif(M1), error = function(e) NA)
if (!all(is.na(VIF_M1))) { cat("\nVIF (M1):\n"); print(VIF_M1) }

# ------------------ selection: LASSO (or stepwise fallback) ----------------------
# LASSO with baseline covariates *unpenalized* so they are always kept.
if (HAS_GLMNET) {
  library(glmnet)
  
  # Build a model matrix and use 'assign' to map columns back to term names.
  all_preds <- attr(terms(M1), "term.labels")
  Terms     <- terms(reformulate(all_preds))
  MM        <- model.matrix(Terms, data = df)   # includes intercept
  assignVec <- attr(MM, "assign")
  termLabs  <- attr(Terms, "term.labels")
  
  X       <- MM[, -1, drop = FALSE]            # drop intercept column
  assignX <- assignVec[-1]
  
  # Unpenalize any column that belongs to a baseline covariate (handles factor dummies).
  pf <- rep(1, ncol(X)); names(pf) <- colnames(X)
  base_idx <- which(termLabs %in% base_covars)
  pf[assignX %in% base_idx] <- 0
  
  set.seed(123)
  cvfit <- cv.glmnet(x = X, y = df$GE, alpha = 1, family = "gaussian",
                     nfolds = 10, penalty.factor = pf)
  
  cat("\n=== LASSO CV curve (printed) ===\n")
  plot(cvfit)  # Plots pane
  if (SAVE_OUTPUTS) {
    out_dir <- "task1_models"
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    png(file.path(out_dir, "lasso_cv_curve.png"), width = 1200, height = 900, res = 180)
    plot(cvfit); dev.off()
  }
  
  # 1-SE rule → selected terms; map back from columns to term labels.
  lam1se <- cvfit$lambda.1se
  beta   <- as.numeric(coef(cvfit, s = lam1se))[-1]
  keep_cols <- which(beta != 0)
  sel_terms <- unique(termLabs[ assignX[keep_cols] ])
  
  cat("\nSelected (λ_1SE):\n"); print(sel_terms)
  if (SAVE_OUTPUTS) {
    out_dir <- "task1_models"
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(data.frame(selected = sel_terms),
                     file.path(out_dir, "lasso_selected_lambda1se.csv"))
  }
  
  final_terms <- unique(c(base_covars, sel_terms))
  form_final  <- reformulate(final_terms, response = "GE")
  Mfinal      <- lm(form_final, data = df)
  
} else {
  # Fallback: stepwise selection (AIC) with baseline terms forced in.
  library(MASS)
  scope  <- list(lower = reformulate(base_covars, response = "GE"), upper = form_M1)
  Mfinal <- stepAIC(M1, scope = scope, direction = "both", trace = FALSE)
}

# ------------------ final model outputs ------------------------------------------
cat("\n=== Final model (Mfinal) — glance ===\n"); print(glance(Mfinal))
cat("\n=== Final model (Mfinal) — robust coefficients ===\n"); print(robust_table(Mfinal))

# Diagnostics: residuals, QQ, scale-location, Cook's distance (prints to device).
cat("\n=== Mfinal diagnostics (printed) ===\n")
op <- par(mfrow = c(2,2)); on.exit(par(op), add = TRUE); plot(Mfinal)

if (SAVE_OUTPUTS) {
  out_dir <- "task1_models"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  png(file.path(out_dir, "Mfinal_diagnostics.png"), width = 1600, height = 1600, res = 180)
  par(mfrow = c(2,2)); plot(Mfinal); dev.off()
}

# ============================================================
# VISUALISATIONS 
# ============================================================
# 1) Coefficient forest (robust CI), 2) Obs vs Pred (in-sample & CV),
# 3) Calibration curve, 4) Partial effects (visreg or added-variable plots).

# 1) Coefficient forest (HC3 robust CI)
coef_final <- robust_table(Mfinal) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = forcats::fct_reorder(term, estimate))

p_coef <- ggplot(coef_final, aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.3) +
  coord_flip() +
  labs(title = "Final model coefficients (HC3 robust 95% CI)", x = NULL, y = "β")
print(p_coef)
if (SAVE_OUTPUTS) {
  out_dir <- "task1_models"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, "Mfinal_coef_forest.png"), p_coef, width = 8, height = 6, dpi = 300)
}

# --- helpers for CV predictions & plots ------------------------------------------
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
         x = "Predicted", y = "Observed GE")
}

# 2) Observed vs Predicted (in-sample & 10-fold CV) for M0 & Mfinal ---------------
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
  print((p1 | p2) / (p3 | p4))
  if (SAVE_OUTPUTS) {
    out_dir <- "task1_models"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(out_dir, "obs_vs_pred_M0_Mfinal.png"), (p1 | p2) / (p3 | p4),
           width = 10, height = 8, dpi = 300)
  }
} else {
  print(p1); print(p2); print(p3); print(p4)
  if (SAVE_OUTPUTS) {
    out_dir <- "task1_models"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(out_dir, "obs_pred_M0_in.png"),     p1, width=6, height=5, dpi=300)
    ggsave(file.path(out_dir, "obs_pred_M0_cv.png"),     p2, width=6, height=5, dpi=300)
    ggsave(file.path(out_dir, "obs_pred_Mfinal_in.png"), p3, width=6, height=5, dpi=300)
    ggsave(file.path(out_dir, "obs_pred_Mfinal_cv.png"), p4, width=6, height=5, dpi=300)
  }
}

# 3) Calibration (deciles) for Mfinal (CV preds) ----------------------------------
cal_df <- tibble(obs = df$GE, pred = df$.pred_Mfinal_cv) %>%
  tidyr::drop_na() %>%
  mutate(bin = ntile(pred, 10)) %>%
  group_by(bin) %>% summarise(pred = mean(pred), obs = mean(obs), .groups = "drop")

p_cal <- ggplot(cal_df, aes(pred, obs)) +
  geom_abline(linetype = 2) +
  geom_point(size = 2) + geom_line() +
  labs(title = "Calibration of final model (10-fold CV)",
       x = "Mean predicted GE (by decile)", y = "Mean observed GE")
print(p_cal)
if (SAVE_OUTPUTS) {
  out_dir <- "task1_models"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, "Mfinal_calibration.png"), p_cal, width = 7, height = 5, dpi = 300)
}

# 4) Partial effects for most influential numeric predictors ----------------------
# We standardize numeric predictors to compare absolute beta magnitudes,
# pick top 4 (excluding binary dummies), and display their adjusted partial effects.
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
      labs(title = paste("Adjusted effect:", v), x = v, y = "GE (partial)")
  })
  if (HAS_PATCHWORK) {
    library(patchwork)
    print(wrap_plots(ps, ncol = 2))
    if (SAVE_OUTPUTS) {
      out_dir <- "task1_models"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      ggsave(file.path(out_dir, "Mfinal_partial_effects.png"),
             wrap_plots(ps, ncol = 2), width = 10, height = 8, dpi = 300)
    }
  } else {
    for (p in ps) print(p)
  }
} else {
  cat("\n(visreg not installed) Showing added-variable plots via car::avPlots\n")
  car::avPlots(Mfinal, terms = top4)  # prints to device
  if (SAVE_OUTPUTS) {
    out_dir <- "task1_models"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    png(file.path(out_dir, "Mfinal_added_variable_top4.png"), width = 1600, height = 1200, res = 180)
    par(mfrow = c(2,2)); car::avPlots(Mfinal, terms = top4); dev.off()
  }
}

# ===============================
# Reporting add-ons
# ===============================
# Compare M0 vs Mfinal using adjusted R², AIC, BIC, RMSE (in-sample) and 10-fold CV RMSE.

# simple RMSE helpers --------------------------------------------------------------
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

stopifnot(exists("M0"), exists("Mfinal"))
df_final <- model.frame(Mfinal)

tab_models <- tibble::tibble(
  model     = c("M0_baseline", "Mfinal"),
  adj_r2    = c(glance(M0)$adj.r.squared, glance(Mfinal)$adj.r.squared),
  AIC       = c(AIC(M0), AIC(Mfinal)),
  BIC       = c(BIC(M0), BIC(Mfinal)),
  RMSE_in   = c(rmse(df_final$GE, fitted(M0)), rmse(df_final$GE, fitted(Mfinal))),
  RMSE_cv10 = c(cv_rmse_lm(formula(M0), df_final, K = 10, seed = 42),
                cv_rmse_lm(formula(Mfinal), df_final, K = 10, seed = 42))
)
cat("\n=== Model comparison (M0 vs Mfinal) ===\n"); print(tab_models)

# Added-variable plots for selected hormones (printed) -----------------------------
# Partial regression plots constructed manually (rx, ry residuals).

avp <- function(fit, term, data, title_prefix = NULL) {
  all_terms <- attr(terms(fit), "term.labels"); stopifnot(term %in% all_terms)
  others <- setdiff(all_terms, term)
  f_y <- reformulate(others, response = all.vars(formula(fit))[1]); ry <- resid(lm(f_y, data = data))
  f_x <- reformulate(others, response = term);                          rx <- resid(lm(f_x, data = data))
  d <- tibble::tibble(rx = rx, ry = ry)
  ggplot(d, aes(rx, ry)) +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_vline(xintercept = 0, linewidth = 0.2) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = paste(term, "| others"), y = paste(all.vars(formula(fit))[1], "| others"),
         title = paste0(title_prefix %||% "Added-variable: ", term))
}

sel_terms <- attr(terms(Mfinal), "term.labels")
candidate_horms <- c("GLP1_ln","Ghrelin_ln","Amylin_ln","Gastrin","CCK","PYY","Glucagon")
av_terms <- intersect(sel_terms, candidate_horms)
if (length(av_terms)) {
  aps <- lapply(av_terms, function(trm) avp(Mfinal, trm, df_final, title_prefix = "Added-variable: "))
  if (HAS_PATCHWORK) { library(patchwork); print(wrap_plots(aps, ncol = 2)) } else { for (p in aps) print(p) }
}

# Robust coefficient table (printed; optional save) --------------------------------
rob <- coeftest(Mfinal, vcov = sandwich::vcovHC(Mfinal, type = "HC3"))
tbl_coef <- tibble::tibble(
  term = rownames(rob),
  estimate = rob[,1],
  std.error = rob[,2],
  statistic = rob[,3],
  p.value = rob[,4]
) %>% mutate(conf.low = estimate - 1.96*std.error,
             conf.high = estimate + 1.96*std.error)
cat("\n=== Robust coefficients (HC3) ===\n"); print(tbl_coef)

# Plain-language effect-size translations (natural log) ----------------------------
LN_10pct <- log(1.1); LN_2x <- log(2)
describe_effect <- function(term, beta, data) {
  out <- tibble::tibble(term = term, beta = beta,
                        delta_10pct = NA_real_, delta_2x = NA_real_,
                        delta_per1 = NA_real_, delta_perIQR = NA_real_,
                        note = NA_character_)

}
