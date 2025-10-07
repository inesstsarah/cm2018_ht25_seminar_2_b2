# ================================================================
# TASK 1 — SETUP
# ================================================================
# Goal: Explore and model gastric emptying half-life (GE) in T2DM.

# --- Libraries used across EDA, correlations, and modelling -----
library(ggplot2)    
library(dplyr)       
library(readr)       
library(tidyr)     
library(forcats)     
library(scales)     
library(rlang)   
library(broom)    
library(corrplot)    
library(sandwich)    
library(lmtest)   
library(car)     
library(glmnet)      
library(patchwork)  
library(visreg)    
library(purrr)

theme_set(theme_minimal(base_size = 12))

data <- read_csv("/Users/claraazzano/Desktop/Data_T1.csv", show_col_types = FALSE)

# ================================================================
# VARIABLES & LABELS
# ================================================================

cont_vars <- c(
  "GE","Age","Height","BW","BMI",
  "GlucoseFasting","InsulinFasting","HbA1c","MatsudaIdx","HOMAB",
  "Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY"
)

log_suggest <- c("InsulinFasting","MatsudaIdx","Ghrelin","Amylin","GLP1")
binary_vars <- c("Sex","Metformin","DiabetesComplications")

units_map <- c(
  GE = "min", Age = "y", Height = "cm", BW = "kg", BMI = "kg/m^2",
  GlucoseFasting = "mmol/L", InsulinFasting = "mU/L", HbA1c = "%",
  Gastrin = "pg/mL", CCK = "pmol/L", Ghrelin = "pg/mL", Amylin = "pmol/L",
  Glucagon = "pmol/L", GLP1 = "pmol/L", PYY = "ng/L", MatsudaIdx = "", HOMAB = ""
)

label_map <- list(
  Sex = c("Female","Male"),
  Metformin = c("No","Yes"),
  DiabetesComplications = c("No","Yes")
)

get_units <- function(var) {
  nm <- as.character(var)
  if (!is.null(names(units_map)) && nm %in% names(units_map)) units_map[[nm]] else NULL
}

data <- data %>%
  mutate(across(all_of(cont_vars), ~ suppressWarnings(as.numeric(.))))

# ================================================================
# SHARED HELPERS
# ================================================================

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

is_skewed <- function(x, thresh = 1) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || sd(x) == 0) return(FALSE)
  z <- (x - mean(x)) / sd(x)
  abs(mean(z^3)) > thresh
}

safe_log10 <- function(x, eps = 1e-6) {
  x <- as.numeric(x)
  needs_offset <- any(x <= 0, na.rm = TRUE)
  off <- if (needs_offset) max(eps, abs(min(x, na.rm = TRUE)) + eps) else 0
  y <- log10(x + off)
  attr(y, "offset") <- off
  y
}

axis_label <- function(name, log = FALSE, units = NULL) {
  paste0(name,
         if (!is.null(units) && nzchar(units)) paste0(" [", units, "]") else "",
         if (log) " (log10)" else "")
}

as_pretty_factor <- function(x, var, labels = label_map) {
  xn <- x[!is.na(x)]
  if (is.numeric(x) && dplyr::n_distinct(xn) == 2) {
    lv <- sort(unique(xn))
    f  <- factor(x, levels = lv)
    if (!is.null(labels[[var]]) && length(labels[[var]]) == 2) {
      levels(f) <- labels[[var]]
    } else if (all(lv %in% c(0,1))) {
      levels(f) <- c("No","Yes")
    }
    return(f)
  }
  if (is.character(x) || is.factor(x)) {
    f <- forcats::fct_drop(factor(x))
    if (nlevels(f) == 2) return(f)
  }
  NULL
}

cohen_d <- function(y, g) {
  g <- droplevels(g); if (nlevels(g) != 2) return(NA_real_)
  y1 <- y[g == levels(g)[1]]; y2 <- y[g == levels(g)[2]]
  n1 <- sum(is.finite(y1)); n2 <- sum(is.finite(y2))
  s1 <- var(y1, na.rm = TRUE); s2 <- var(y2, na.rm = TRUE)
  sp <- sqrt(((n1 - 1)*s1 + (n2 - 1)*s2) / (n1 + n2 - 2))
  (mean(y2, na.rm = TRUE) - mean(y1, na.rm = TRUE)) / sp
}

plot_group_boxes <- function(df, group_var, y = "GE",
                             labels = label_map,
                             save_dir = NULL, width = 6.5, height = 4.5, dpi = 300) {
  stopifnot(group_var %in% names(df), y %in% names(df))
  g <- as_pretty_factor(df[[group_var]], group_var, labels)
  if (is.null(g)) { message("Skipping '", group_var, "': not 2-level."); return(NULL) }
  
  d <- df %>% dplyr::transmute(group = g, y = !!rlang::sym(y)) %>% tidyr::drop_na()
  if (nrow(d) == 0 || nlevels(d$group) < 2) {
    message("Skipping '", group_var, "': insufficient data after drop_na."); return(NULL)
  }
  
  n_tab <- d %>% dplyr::count(group, name = "n")
  subtitle_txt <- paste(paste0(n_tab$group, ": n=", n_tab$n), collapse = "   |   ")
  
  p <- ggplot(d, aes(group, y)) +
    geom_boxplot(outlier.alpha = 0.35, width = 0.55) +
    geom_jitter(width = 0.08, alpha = 0.35, size = 1) +
    labs(title = paste("GE by", group_var), subtitle = subtitle_txt,
         x = group_var, y = "GE [min]") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  print(p)
  
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(file.path(save_dir, paste0("box_GE_by_", group_var, ".png")),
           plot = p, width = width, height = height, dpi = dpi)
  }
  
  t_out <- t.test(y ~ group, data = d, var.equal = FALSE)
  w_out <- wilcox.test(y ~ group, data = d, exact = FALSE)
  d_out <- cohen_d(d$y, d$group)
  
  tests <- dplyr::bind_rows(
    broom::tidy(t_out) %>% dplyr::mutate(test = "Welch t-test", effect = d_out),
    broom::tidy(w_out) %>% dplyr::mutate(test = "Wilcoxon rank-sum", effect = NA_real_)
  ) %>%
    dplyr::mutate(group_var = group_var) %>%
    dplyr::select(group_var, test, estimate, estimate1, estimate2,
                  statistic, parameter, p.value, conf.low, conf.high, effect, method)
  
  list(plot = p, tests = tests)
}

# ================================================================
# PART A — EXPLORATORY VISUALISATIONS
# ================================================================

plot_hist <- function(df, var, bins = NULL, log_x = NULL,
                      units = NULL, title = NULL, annotate_n = TRUE) {
  v <- enquo(var); vname <- as_name(v)
  x <- df %>% pull(!!v)
  
  # Auto-log if not specified: only if all positive and skewed
  if (is.null(log_x)) log_x <- all(x > 0, na.rm = TRUE) && is_skewed(x)
  
  # Prepare plotting vector (keeps chosen offset so we can annotate it)
  x_plot <- if (log_x) safe_log10(x) else x
  off <- if (log_x) attr(x_plot, "offset") else 0
  df_plot <- tibble(x_plot = as.numeric(x_plot))
  
  # Bin choice: FD by default, or user-specified bin count
  if (is.null(bins)) {
    bw <- fd_binwidth(df_plot$x_plot)
    p <- ggplot(df_plot, aes(x = x_plot)) +
      geom_histogram(binwidth = bw, fill = "#6A9AE2", color = "white", na.rm = TRUE)
  } else {
    p <- ggplot(df_plot, aes(x = x_plot)) +
      geom_histogram(bins = bins, fill = "#6A9AE2", color = "white", na.rm = TRUE)
  }
  
  # Labels (append units and “(log10)” as needed)
  p <- p + labs(
    title = if (is.null(title)) vname else title,
    x = axis_label(vname, log = isTRUE(log_x),
                   units = if (is.null(units)) get_units(vname) else units),
    y = "Count"
  ) + theme(panel.grid.minor = element_blank())
  
  # Annotate with n and missingness; if offset was used, display it explicitly.
  if (annotate_n) {
    n_use <- sum(is.finite(df_plot$x_plot))
    miss  <- sum(!is.finite(df_plot$x_plot))
    ann <- paste0("n=", n_use, if (miss > 0) paste0(" (NA=", miss, ")"))
    if (isTRUE(log_x) && off > 0) ann <- paste0(ann, "\nlog10(x + ", signif(off, 3), ")")
    p <- p + annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5, label = ann)
  }
  p
}

# Print histograms for all continuous variables (to Plots)
invisible(lapply(cont_vars, function(v) {
  print(plot_hist(data, !!sym(v),
                  log_x = if (v %in% log_suggest) TRUE else NULL,
                  units = get_units(v)))
}))

# ---- Binary summaries: neat "n (%)" tables ---------------------
binary_summary <- function(df, vars, labels = NULL) {
  out <- lapply(vars, function(v) {
    x <- df[[v]]
    tab <- tibble(Level = c(0,1)) %>%
      left_join(
        tibble(Level = x) %>% filter(!is.na(Level)) %>% count(Level, name = "n"),
        by = "Level"
      ) %>%
      mutate(
        n   = replace_na(n, 0),
        pct = n / sum(n),
        Variable = v,
        Label = if (!is.null(labels) && !is.null(labels[[v]])) {
          labels[[v]][match(Level, c(0,1))]
        } else as.character(Level)
      ) %>%
      select(Variable, Level, Label, n, pct)
  })
  bind_rows(out) %>% mutate(pct_label = percent(pct, accuracy = 0.1))
}

# Presentable labels for 0/1 indicators
label_map <- list(Sex = c("Female","Male"),
                  Metformin = c("No","Yes"),
                  DiabetesComplications = c("No","Yes"))

# Print long and wide forms
bin_tbl <- binary_summary(data, binary_vars, labels = label_map)
cat("\n--- Binary variables: counts and percentages ---\n")
bin_tbl %>%
  arrange(Variable, Level) %>%
  mutate(`n (%)` = paste0(n, " (", pct_label, ")")) %>%
  select(Variable, Label, `n (%)`) %>% print(n = Inf)

cat("\n--- Binary variables (wide view) ---\n")
bin_tbl %>%
  mutate(`n (%)` = paste0(n, " (", pct_label, ")")) %>%
  select(Variable, Label, `n (%)`) %>%
  pivot_wider(names_from = Label, values_from = `n (%)`) %>%
  print(n = Inf)

# ---- Descriptives and missingness ------------------------------
# Summarise: n, mean, sd, median, IQR, min, max (rounded for readability)
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
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(everything(),
               names_to = c("variable",".value"),
               names_sep = "_") %>%
  mutate(across(c(mean, sd, median, IQR, min, max), ~round(., 2)))
cat("\n--- Descriptive stats (continuous) ---\n"); print(desc, n = Inf)

# Missingness: how many NAs per variable (descending)
missing_tbl <- data %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to="variable", values_to="n_missing") %>%
  arrange(desc(n_missing))
cat("\n--- Missingness by variable ---\n"); print(missing_tbl, n = Inf)

# ================================================================
# PART B — CORRELATIONS
# ================================================================
# Variable blocks for structured correlation checks vs GE
vars_hormones <- c("Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY")
vars_demo     <- c("Age","Height","BW","BMI")
vars_gly      <- c("GlucoseFasting","HbA1c")
vars_insulin  <- c("InsulinFasting","MatsudaIdx","HOMAB")
vars_outcome  <- "GE"

# EDA-only log10 transforms to tame skew (keeps interpretation on a monotone scale)
log_skewed <- c("InsulinFasting","MatsudaIdx","Ghrelin","Amylin","GLP1")

# Prepare a numeric-only, zero-variance-free, optionally log-transformed data frame
make_pairs_df <- function(df, keep, log10_vars = character()) {
  keep <- intersect(keep, names(df))
  d <- df[keep]
  d <- d[vapply(d, is.numeric, logical(1L))]         # numeric only
  for (v in intersect(log10_vars, names(d))) d[[v]] <- safe_log10(d[[v]])
  if (ncol(d)) {
    ok <- vapply(d, function(z) sd(z, na.rm = TRUE) > 0, logical(1L))
    d <- d[, ok, drop = FALSE]                       # drop constants
  }
  d
}

# Pairwise Spearman with BH FDR; returns tidy tibble sorted by |rho|
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

# Custom pairs-plot (lower: points; upper: shaded by rho; diag: densities)
pairs_plot <- function(df, vars, log10_vars = character(),
                       group = NULL, main = "") {
  d <- make_pairs_df(df, keep = vars, log10_vars = log10_vars)
  if (ncol(d) < 2L) {
    message("Not enough numeric variables to draw pairs(): ", paste(vars, collapse = ", "))
    return(invisible(NULL))
  }
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

# --- Global Spearman correlation heatmap on all numeric, non-constant ----------
num <- data %>% select(where(is.numeric)) %>% select(where(~ sd(., na.rm = TRUE) > 0))
M <- cor(num, method = "spearman", use = "pairwise.complete.obs")
cat("\n--- Global Spearman correlation heatmap (printed) ---\n")
corrplot(M, method = "color", type = "upper", diag = FALSE,
         addCoef.col = "black", number.digits = 2, number.cex = 0.5,
         tl.cex = 0.7, order = "hclust")

# --- GE vs blocks: pairplots + top correlations ----------------
run_block <- function(df, left, right, log_vars, prefix, group = NULL) {
  title <- gsub("_", " ", toupper(prefix))
  cat("\n--- Pairplot:", title, "---\n")
  pairs_plot(df, vars = c(left, right), log10_vars = log_vars, group = group, main = title)
  tab <- spearman_table(df, vars = c(left, right), log10_vars = log_vars, sort_by = "abs_r")
  cat("\nTop correlations (|rho|):", title, "\n"); print(head(tab, 20))
  invisible(tab)
}
tab_GE_horm <- run_block(data, vars_outcome, vars_hormones, log_skewed, "GE_hormones")
tab_GE_gly  <- run_block(data, vars_outcome, vars_gly,      log_skewed, "GE_glycemia")
tab_GE_ins  <- run_block(data, vars_outcome, vars_insulin,  log_skewed, "GE_insulin")
tab_GE_demo <- run_block(data, vars_outcome, vars_demo,     log_skewed, "GE_demographics")

# --- GE correlation dot-plot summary: rho + FDR per block -------
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
print(head(res, 20))

# Reorder within each facet by |rho| for a clean top-to-bottom reading
res <- res %>%
  dplyr::group_by(block) %>%
  dplyr::mutate(var = forcats::fct_reorder(var, abs(rho))) %>%
  dplyr::ungroup()

p_dot <- ggplot(res, aes(x = var, y = rho, shape = padj < 0.05)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_point(size = 2.5) +
  coord_flip() +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_shape_manual(values = c(16, 17),
                     labels = c("FDR ≥ 0.05", "FDR < 0.05"), name = NULL) +
  labs(title = "Spearman correlation with GE",
       x = NULL, y = "rho (Spearman)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
print(p_dot)

# --- Top-4 scatter checks with LOESS trend and rho annotation ---
top4 <- res %>%
  arrange(desc(abs(rho))) %>%
  filter(var %in% names(data)) %>%
  slice(1:4) %>%
  pull(var) %>% as.character()

plot_scatter <- function(df, xvar, yvar = "GE") {
  x <- suppressWarnings(as.numeric(df[[xvar]]))
  y <- suppressWarnings(as.numeric(df[[yvar]]))
  ok <- is.finite(x) & is.finite(y)
  r  <- if (sum(ok) >= 3L) suppressWarnings(cor(x[ok], y[ok], method = "spearman")) else NA_real_
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_point(alpha = 0.35, size = 1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 0.6) +
    labs(x = xvar, y = yvar,
         subtitle = paste0("rho = ", if (is.finite(r)) sprintf("%.2f", r) else "NA")) +
    theme(panel.grid.minor = element_blank())
}

if (length(top4) >= 4) {
  p1 <- plot_scatter(data, top4[1]); p2 <- plot_scatter(data, top4[2])
  p3 <- plot_scatter(data, top4[3]); p4 <- plot_scatter(data, top4[4])
  print((p1 | p2) / (p3 | p4))   # use patchwork grid
} else {
  for (v in top4) print(plot_scatter(data, v))
}



# ================================================================
# PART C — GROUP DIFFERENCES (Boxplots: GE ~ Sex/Metformin/Complications)
# ================================================================

cat("\n--- Boxplots & group comparisons ---\n")

group_vars <- c("Sex", "Metformin", "DiabetesComplications")

tests_all <- purrr::map_dfr(group_vars, function(gv) {
  if (!gv %in% names(data)) { message("'", gv, "' not found."); return(NULL) }
  out <- plot_group_boxes(data, gv, y = "GE", labels = label_map)
  if (is.null(out)) return(NULL)
  cat("\n-- Tests for", gv, "--\n"); print(out$tests); cat("\n")
  out$tests
})

if (nrow(tests_all) > 0) {
  tests_all <- tests_all %>%
    group_by(test) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    ungroup()
  cat("\n=== Summary across splits (BH-FDR) ===\n"); print(tests_all, n = Inf)
} else {
  message("No group tests available.")
}

# ================================================================
# PART D — MODELLING GE
# ================================================================
# Strategy:
#   1) Baseline clinical covariates (M0).
#   2) Add glycaemia/insulin + hormones (M1).
#   3) Use LASSO (with baseline unpenalized) to pick a parsimonious set → refit OLS (Mfinal).
#   4) Robust inference (HC3), diagnostics, CV performance, calibration, partial effects.

# Choose representatives for correlated blocks (helps interpretability)
gly_marker <- "HbA1c"          # representative of glycaemia (long-term)
ins_marker <- "InsulinFasting" # representative of insulin axis (level)

# Transform skewed hormones for modelling (natural log) and cast binaries to factors
dat <- data %>%
  mutate(
    GLP1_ln    = ifelse(GLP1    > 0, log(GLP1),    NA_real_),
    Ghrelin_ln = ifelse(Ghrelin > 0, log(Ghrelin), NA_real_),
    Amylin_ln  = ifelse(Amylin  > 0, log(Amylin),  NA_real_),
    Sex                   = if ("Sex" %in% names(.)) factor(Sex) else NULL,
    Metformin             = if ("Metformin" %in% names(.)) factor(Metformin) else NULL,
    DiabetesComplications = if ("DiabetesComplications" %in% names(.)) factor(DiabetesComplications) else NULL
  )

# Define modelling sets
base_covars <- c("Age","Sex","BMI","Metformin")
extra_core  <- c(gly_marker, ins_marker)
hormones_ln <- c("GLP1_ln","Ghrelin_ln","Amylin_ln","Gastrin","CCK","Glucagon","PYY")

# Build modelling frame and drop rows missing any included predictor (complete-case)
# NOTE: If missingness is not MCAR, consider multiple imputation in future work.
all_vars <- c("GE", base_covars, extra_core, hormones_ln)
df <- dat %>% select(any_of(all_vars)) %>% tidyr::drop_na()

# Robust coefficient table helper (HC3 SEs + 95% CI)
robust_table <- function(fit) {
  ct <- coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
  tb <- data.frame(term = rownames(ct), estimate = ct[,1], std.error = ct[,2],
                   statistic = ct[,3], p.value = ct[,4], row.names = NULL)
  tb$conf.low  <- tb$estimate - 1.96*tb$std.error
  tb$conf.high <- tb$estimate + 1.96*tb$std.error
  tb
}

# --- M0 (baseline) and M1 (extended) ----------------------------
form_M0 <- as.formula(paste("GE ~", paste(base_covars, collapse = " + ")))
form_M1 <- as.formula(paste("GE ~", paste(c(base_covars, extra_core, hormones_ln), collapse = " + ")))

M0 <- lm(form_M0, data = df)
M1 <- lm(form_M1, data = df)

cat("\n=== Baseline model (M0) ===\n"); print(glance(M0)); cat("\n"); print(robust_table(M0))
cat("\n=== Extended model (M1) ===\n"); print(glance(M1))

# Multicollinearity check on M1 (VIF > ~5–10 indicates potential concern)
VIF_M1 <- tryCatch(car::vif(M1), error = function(e) NA)
if (!all(is.na(VIF_M1))) { cat("\nVIF (M1):\n"); print(VIF_M1) }

# --- LASSO selection with baseline unpenalized ------------------
# Build model matrix for all predictors, keep mapping to original terms
all_preds <- attr(terms(M1), "term.labels")
Terms     <- terms(reformulate(all_preds))
MM        <- model.matrix(Terms, data = df)      # includes intercept
assignVec <- attr(MM, "assign")                  # each column's term index
termLabs  <- attr(Terms, "term.labels")

X       <- MM[, -1, drop = FALSE]                # drop intercept column
assignX <- assignVec[-1]

# Penalty factors: 0 = never penalize (force keep), 1 = normal penalty
pf <- rep(1, ncol(X)); names(pf) <- colnames(X)
base_idx <- which(termLabs %in% base_covars)
pf[assignX %in% base_idx] <- 0                   # keep baseline terms regardless of lambda

set.seed(123)
cvfit <- cv.glmnet(x = X, y = df$GE, alpha = 1, family = "gaussian",
                   nfolds = 10, penalty.factor = pf)
cat("\n=== LASSO CV curve (printed) ===\n"); plot(cvfit)

# Use 1-SE rule (simpler model with comparable CV error)
lam1se    <- cvfit$lambda.1se
beta      <- as.numeric(coef(cvfit, s = lam1se))[-1]
keep_cols <- which(beta != 0)
sel_terms <- unique(termLabs[ assignX[keep_cols] ])
cat("\nSelected (lambda_1se):\n"); print(sel_terms)

# Final OLS refit for unbiased coefficients + robust inference
final_terms <- unique(c(base_covars, sel_terms))
form_final  <- reformulate(final_terms, response = "GE")
Mfinal      <- lm(form_final, data = df)

# --- Final model outputs ----------------------------------------
cat("\n=== Final model (Mfinal) — glance ===\n"); print(glance(Mfinal))
cat("\n=== Final model (Mfinal) — robust coefficients ===\n"); print(robust_table(Mfinal))

# Diagnostics: residual patterns, QQ, scale-location, leverage/Cook's
cat("\n=== Mfinal diagnostics (printed) ===\n")
op <- par(mfrow = c(2,2)); on.exit(par(op), add = TRUE); plot(Mfinal)

# --- Coefficient forest (robust 95% CI) -------------------------
coef_final <- robust_table(Mfinal) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = fct_reorder(term, estimate))
p_coef <- ggplot(coef_final, aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.3) +
  coord_flip() +
  labs(title = "Final model coefficients (HC3 robust 95% CI)", x = NULL, y = "beta")
print(p_coef)

# ================================================================
# PERFORMANCE & INTERPRETATION
# ================================================================

# K-fold CV helpers and observed-vs-predicted plotting
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
         subtitle = sprintf("%s  •  RMSE=%.2f,  R^2=%.3f", subtitle, m$rmse, m$r2),
         x = "Predicted", y = "Observed GE")
}

# Compare in-sample vs 10-fold CV for M0 and Mfinal
df$.pred_M0_in     <- fitted(M0)
df$.pred_M0_cv     <- cv_predict_lm(form_M0, df)
df$.pred_Mfinal_in <- fitted(Mfinal)
df$.pred_Mfinal_cv <- cv_predict_lm(formula(Mfinal), df)

p1 <- plot_obs_pred(df$GE, df$.pred_M0_in,     "Baseline (M0)",        "In-sample")
p2 <- plot_obs_pred(df$GE, df$.pred_M0_cv,     "Baseline (M0)",        "10-fold CV")
p3 <- plot_obs_pred(df$GE, df$.pred_Mfinal_in, "Final model (Mfinal)", "In-sample")
p4 <- plot_obs_pred(df$GE, df$.pred_Mfinal_cv, "Final model (Mfinal)", "10-fold CV")
print((p1 | p2) / (p3 | p4))

# Calibration curve: mean observed vs mean predicted within deciles
cal_df <- tibble(obs = df$GE, pred = df$.pred_Mfinal_cv) %>%
  tidyr::drop_na() %>%
  mutate(bin = ntile(pred, 10)) %>%
  group_by(bin) %>% summarise(pred = mean(pred), obs = mean(obs), .groups = "drop")
p_cal <- ggplot(cal_df, aes(pred, obs)) +
  geom_abline(linetype = 2) + geom_point(size = 2) + geom_line() +
  labs(title = "Calibration of final model (10-fold CV)",
       x = "Mean predicted GE (by decile)", y = "Mean observed GE")
print(p_cal)

# ---- Adjusted (partial) effects for top standardized predictors ---------------
# Standardize numerics for comparability; pick top 4 by |beta| (excluding dummies)
std_df <- df
num_vars <- names(std_df)[sapply(std_df, is.numeric)]
num_vars <- setdiff(num_vars, "GE")
std_df[num_vars] <- lapply(std_df[num_vars], scale)

M_std <- lm(formula(Mfinal), data = std_df)
std_tab <- tidy(M_std) %>% filter(term != "(Intercept)") %>% arrange(desc(abs(estimate)))
top4 <- head(std_tab$term[!grepl("^Sex|^Metformin", std_tab$term)], 4)

ps <- lapply(top4, function(v) {
  visreg(Mfinal, v, gg = TRUE, overlay = FALSE, rug = TRUE) +
    labs(title = paste("Adjusted effect:", v), x = v, y = "GE (partial)")
})
print(wrap_plots(ps, ncol = 2))

# ---- Compact model comparison (report-ready table) --------------
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

# ---- Added-variable (partial regression) plots for selected hormones ----------
# Visual cross-check for unique contribution of hormone terms in Mfinal.
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
    labs(x = paste(term, "| others"),
         y = paste(all.vars(formula(fit))[1], "| others"),
         title = paste0(title_prefix %||% "Added-variable: ", term))
}
sel_terms <- attr(terms(Mfinal), "term.labels")
candidate_horms <- c("GLP1_ln","Ghrelin_ln","Amylin_ln","Gastrin","CCK","PYY","Glucagon")
av_terms <- intersect(sel_terms, candidate_horms)
if (length(av_terms)) {
  aps <- lapply(av_terms, function(trm) avp(Mfinal, trm, df_final, title_prefix = "Added-variable: "))
  print(wrap_plots(aps, ncol = 2))
}

# ============================== END TASK 1 =============================




