
# THIS IS FOR VIZUALIZATIONS !!!!!!!!! NOT FINAL CODE - NEEDS CELANING AND BETTER COMMENTING

data <- read.csv("/Users/claraazzano/Desktop/Data_T1.csv")

data <- data %>% select(-X)   # drop the redundant index column

head(data)
summary(data)

library(ggplot2)
library(rlang)
library(dplyr)
library(GGally)

#-----DATA VIZUALISATION------------------------------------------------


#MAYBE???

## ======================================================================
## EDA: distributions (continuous + binary) for Data_T1.csv
## ======================================================================

## ---- packages ----
# install.packages(c("readr","dplyr","ggplot2","rlang","forcats","skimr","scales"))
library(readr)
library(dplyr)
library(ggplot2)
library(rlang)
library(forcats)
library(skimr)
library(scales)

theme_set(theme_minimal(base_size = 12))

## ---- load data (keep original column names) ----
# change the path if needed
data <- read_csv("/Users/claraazzano/Desktop/Data_T1.csv", show_col_types = FALSE)

# quick structure + skim
str(data)
skimr::skim(data)

## ======================================================================
## Helpers
## ======================================================================

## -- sensible binwidth (Freedman–Diaconis) --
fd_binwidth <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  bw <- 2 * IQR(x) / (length(x)^(1/3))
  if (!is.finite(bw) || bw <= 0) bw <- diff(range(x)) / 30
  bw
}

## -- quick skew detector (for auto-log suggestion) --
is_skewed <- function(x, thresh = 1) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(FALSE)
  m3 <- mean((x - mean(x))/sd(x))^3
  isTRUE(abs(m3) > thresh)
}

## -- main histogram helper (density scaled to counts) --
#   log_x: TRUE/FALSE/NULL(=auto if skewed & >0)
#   bins:  integer #bins (overrides fd binwidth). leave NULL for FD.
plot_hist_density <- function(df, var, bins = NULL, log_x = NULL, title = NULL,
                              annotate_n = TRUE, density_adjust = 1) {
  v_sym <- ensym(var); v_str <- as_string(v_sym)
  x <- df[[v_str]]
  
  # decide log if not specified
  if (is.null(log_x)) log_x <- all(x > 0, na.rm = TRUE) && is_skewed(x)
  
  x_plot <- if (log_x) log10(x) else x
  df_plot <- data.frame(x_plot = x_plot)
  
  # choose histogram geometry: binwidth (FD) or bins
  if (is.null(bins)) {
    bw <- fd_binwidth(df_plot$x_plot)
    p <- ggplot(df_plot, aes(x = x_plot)) +
      geom_histogram(binwidth = bw, fill = "#6A9AE2", color = "white", na.rm = TRUE)
    scale_const <- sum(!is.na(df_plot$x_plot)) * bw
  } else {
    rng <- range(df_plot$x_plot, na.rm = TRUE)
    bw <- diff(rng) / bins
    p <- ggplot(df_plot, aes(x = x_plot)) +
      geom_histogram(bins = bins, fill = "#6A9AE2", color = "white", na.rm = TRUE)
    scale_const <- sum(!is.na(df_plot$x_plot)) * bw
  }
  
  p <- p +
    geom_density(aes(y = after_stat(density) * scale_const),
                 adjust = density_adjust, linewidth = 1, color = "black", na.rm = TRUE) +
    labs(
      title = title %||% v_str,
      x = paste0(v_str, if (log_x) " (log10)" else ""),
      y = "Count"
    ) +
    theme(panel.grid.minor = element_blank())
  
  if (annotate_n) {
    n  <- sum(is.finite(x_plot))
    miss <- sum(!is.finite(x_plot))
    p <- p + annotate("text", x = Inf, y = Inf,
                      label = paste0("n=", n, if (miss > 0) paste0(" (NA=", miss, ")")),
                      hjust = 1.1, vjust = 1.5, size = 3.5)
  }
  
  p
}

## -- optional: same histogram faceted by a binary/categorical group --
plot_hist_by_group <- function(df, var, group, log_x = NULL, density_adjust = 1) {
  v <- ensym(var); g <- ensym(group)
  v_str <- as_string(v)
  x <- df[[v_str]]
  
  if (is.null(log_x)) log_x <- all(x > 0, na.rm = TRUE) && is_skewed(x)
  x_plot <- if (log_x) log10(x) else x
  
  d <- df %>% mutate(x_plot = x_plot, .grp = as.factor(!!g))
  bw <- fd_binwidth(d$x_plot); n <- sum(!is.na(d$x_plot)); scale_const <- n * bw
  
  ggplot(d, aes(x = x_plot)) +
    geom_histogram(binwidth = bw, fill = "#6A9AE2", color = "white") +
    geom_density(aes(y = after_stat(density) * scale_const), adjust = density_adjust) +
    facet_wrap(vars(.grp), nrow = 1) +
    labs(title = v_str, x = paste0(v_str, if (log_x) " (log10)" else ""), y = "Count") +
    theme(panel.grid.minor = element_blank())
}

## -- binary pie with % labels + NA option --
plot_binary_pie <- function(df, var,
                            labels = c("0","1"),
                            title = NULL,
                            palette = c("#6A9AE2", "#9BC3FF", "#9E9E9E"),
                            include_na = FALSE,
                            na_label = "Missing",
                            min_pct_label = 3) {
  v_sym <- ensym(var); v_str <- as_string(v_sym)
  d <- df %>% mutate(.g = !!v_sym)
  
  if (include_na) {
    d <- d %>% mutate(.g = forcats::fct_explicit_na(factor(.g, levels = c(0,1), labels = labels),
                                                    na_level = na_label))
  } else {
    d <- d %>% filter(!is.na(.g)) %>%
      mutate(.g = factor(.g, levels = c(0,1), labels = labels))
  }
  
  tab <- d %>%
    count(.g, name = "n") %>%
    mutate(pct = 100 * n / sum(n),
           lbl = paste0(sprintf("%.1f", pct), "% (", format(n, big.mark=","), ")"),
           show = pct >= min_pct_label)
  
  ggplot(tab, aes(x = "", y = n, fill = .g)) +
    geom_col(color = "white", width = 1) +
    coord_polar(theta = "y") +
    geom_text(data = subset(tab, show), aes(label = lbl),
              position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = palette[seq_len(nrow(tab))], name = NULL) +
    labs(title = title %||% v_str, x = NULL, y = NULL,
         caption = if (any(!tab$show))
           paste0("Labels suppressed for slices < ", min_pct_label, "%") else NULL) +
    theme_void(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
}

## -- optional: binary bar (% with counts) --
plot_binary_bar <- function(df, var, labels = c("0","1"), title = NULL) {
  v <- ensym(var); v_str <- as_string(v)
  d <- df %>% filter(!is.na(!!v)) %>%
    mutate(.g = factor(!!v, levels = c(0,1), labels = labels)) %>%
    count(.g) %>%
    mutate(pct = n/sum(n), lbl = percent(pct, accuracy = 0.1))
  
  ggplot(d, aes(x = .g, y = pct, fill = .g)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0(lbl, " (", n, ")")), vjust = -0.3) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1.05)) +
    scale_fill_manual(values = c("#6A9AE2", "#9BC3FF"), guide = "none") +
    labs(title = title %||% v_str, x = NULL, y = "Percentage")
}

## ======================================================================
## Which variables to plot
## ======================================================================

cont_vars <- c(
  "GE","Age","Height","BW","BMI",
  "GlucoseFasting","InsulinFasting","HbA1c","MatsudaIdx","HOMAB",
  "Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY"
)

# variables you *know* are better on log scale for EDA
log_suggest <- c("InsulinFasting","MatsudaIdx","Ghrelin","Amylin","GLP1")

binary_vars <- c(
  "Sex",                   # 0/1 (Female/Male) in your dataset
  "Metformin",             # 0/1
  "DiabetesComplications"  # 0/1
)

## ======================================================================
## Draw histograms (to screen)
## ======================================================================

invisible(lapply(cont_vars, function(v) {
  print(
    plot_hist_density(
      data, !!sym(v),
      log_x = if (v %in% log_suggest) TRUE else NULL  # TRUE for suggested; auto for others
    )
  )
}))

## (optional) a couple faceted examples
# print(plot_hist_by_group(data, GLP1, Metformin))
# print(plot_hist_by_group(data, GE, Sex))

## ======================================================================
## Binary plots (pies by assignment, bar alternative)
## ======================================================================

# pies (drop NA from percentages)
print(plot_binary_pie(data, Sex, labels = c("Female","Male"), title = "Sex"))
print(plot_binary_pie(data, Metformin, labels = c("No","Yes"), title = "Metformin"))
print(plot_binary_pie(data, DiabetesComplications, labels = c("No","Yes"),
                      title = "DiabetesComplications"))

# include NA slice version (if you want one example)
# print(plot_binary_pie(data, DiabetesComplications, labels = c("No","Yes"),
#                       include_na = TRUE, na_label = "NA", title = "DiabetesComplications"))

# bar alternatives (often clearer than pies)
# print(plot_binary_bar(data, Sex, labels = c("Female","Male")))
# print(plot_binary_bar(data, Metformin, labels = c("No","Yes")))
# print(plot_binary_bar(data, DiabetesComplications, labels = c("No","Yes")))

## ======================================================================
## Save all continuous histograms to files (appendix-ready)
## ======================================================================

outdir <- "fig_histograms"
if (!dir.exists(outdir)) dir.create(outdir)

invisible(lapply(cont_vars, function(v) {
  g <- plot_hist_density(
    data, !!sym(v),
    log_x = if (v %in% log_suggest) TRUE else NULL
  )
  ggsave(filename = file.path(outdir, paste0("hist_", v, ".png")),
         plot = g, width = 8, height = 5, dpi = 300)
}))

## save pies, too (optional)
outdir_bin <- "fig_binary"
if (!dir.exists(outdir_bin)) dir.create(outdir_bin)

ggsave(file.path(outdir_bin, "pie_Sex.png"),
       plot_binary_pie(data, Sex, labels = c("Female","Male"), title = "Sex"),
       width = 6, height = 6, dpi = 300)
ggsave(file.path(outdir_bin, "pie_Metformin.png"),
       plot_binary_pie(data, Metformin, labels = c("No","Yes"), title = "Metformin"),
       width = 6, height = 6, dpi = 300)
ggsave(file.path(outdir_bin, "pie_DiabetesComplications.png"),
       plot_binary_pie(data, DiabetesComplications, labels = c("No","Yes"), title = "DiabetesComplications"),
       width = 6, height = 6, dpi = 300)

## ======================================================================
## Compact descriptive stats for continuous variables
## ======================================================================

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
  tidyr::pivot_longer(everything(),
                      names_to = c("variable",".value"),
                      names_sep = "_") %>%
  mutate(across(c(mean, sd, median, IQR, min, max), ~round(., 2)))

print(desc, n = Inf)
# write_csv(desc, "descriptives_continuous.csv")

#browseURL(file.path(getwd(), "fig_histograms"))
#browseURL(file.path(getwd(), "fig_binary"))

missing_tbl <- data %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  tidyr::pivot_longer(everything(), names_to="variable", values_to="n_missing") %>%
  arrange(desc(n_missing))
print(missing_tbl, n = Inf)

# CORRELATIONS PLOTS

#=========================
  #==========================================

  
data <- read.csv("/Users/claraazzano/Desktop/Data_T1.csv")

library(ggplot2)
library(rlang)
library(dplyr)
library(GGally)


# Scatter plot correlation--------------------------------

## ---------- panels (no r/p text) ----------

# diagonal density
.panel.density <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  d <- density(x[is.finite(x)])
  lines(d$x, d$y / max(d$y), lwd = 1.2)
}

# lower: scatter with light alpha
.panel.points <- function(x, y, col = "black", pch = 19, ...) {
  ok <- is.finite(x) & is.finite(y)
  points(x[ok], y[ok], pch = pch,
         col = adjustcolor(col, alpha.f = 0.35), cex = 0.6)
}

# upper: just a background tint by |rho| (no text)
.panel.shade <- function(x, y, palette = c("#B71C1C","white","#1B5E20"), ...) {
  ok <- is.finite(x) & is.finite(y)
  r  <- if (sum(ok) >= 3) suppressWarnings(cor(x[ok], y[ok], method = "spearman")) else NA_real_
  col_bg <- colorRampPalette(palette)(200)
  idx <- if (is.finite(r)) max(1, min(200, round((r + 1)/2 * 199 + 1))) else 100
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = col_bg[idx], border = NA)
}

# build plotting data.frame with optional log10 transform
make_pairs_df <- function(df, keep, log10_vars = character()) {
  keep <- intersect(keep, names(df))
  d <- df[keep]
  d <- d[vapply(d, is.numeric, logical(1L))]
  for (v in intersect(log10_vars, names(d))) d[[v]] <- log10(d[[v]])
  if (ncol(d)) {
    sd_ok <- vapply(d, function(z) sd(z, na.rm = TRUE) > 0, logical(1L))
    d <- d[, sd_ok, drop = FALSE]
  }
  d
}

# wrapper: pairs plot with our panels (no numbers)
pairs_chunk <- function(df, vars, log10_vars = character(),
                        group = NULL, group_colors = c("#6A9AE2", "#9BC3FF"),
                        main = NULL) {
  d <- make_pairs_df(df, keep = vars, log10_vars = log10_vars)
  if (!is.null(group)) {
    g <- as.factor(group)
    cols <- adjustcolor(group_colors, alpha.f = 0.8)
    col_vec <- cols[as.integer(g)]
  } else col_vec <- "black"
  
  pairs(d,
        lower.panel = function(x, y, ...) .panel.points(x, y, col = col_vec),
        upper.panel = .panel.shade,
        diag.panel  = .panel.density,
        gap = 0.6,
        cex.labels = 1.0,
        main = main)
}

## ---------- tidy Spearman table (ρ, p, N) ----------

spearman_table <- function(df, vars, log10_vars = character(),
                           sort_by = c("abs_r","p","r")) {
  d <- make_pairs_df(df, keep = vars, log10_vars = log10_vars)
  cn <- colnames(d); k <- length(cn)
  out <- list()
  for (i in seq_len(k - 1)) {
    for (j in (i + 1):k) {
      x <- d[[i]]; y <- d[[j]]
      ok <- is.finite(x) & is.finite(y)
      n  <- sum(ok)
      if (n >= 3) {
        ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
        out[[length(out) + 1]] <- data.frame(
          var1 = cn[i],
          var2 = cn[j],
          rho  = unname(ct$estimate),
          p    = unname(ct$p.value),
          n    = n,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  tab <- do.call(rbind, out)
  sort_by <- match.arg(sort_by)
  ord <- switch(sort_by,
                abs_r = order(-abs(tab$rho), tab$p),
                p     = order(tab$p, -abs(tab$rho)),
                r     = order(-tab$rho))
  tab[ord, ]
}

# define your groups
vars_hormones <- c("Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY")
vars_demo     <- c("Age","Height","BW","BMI")
vars_gly      <- c("GlucoseFasting","HbA1c")
vars_insulin  <- c("InsulinFasting","MatsudaIdx","HOMAB")
vars_outcome  <- "GE"
log_skewed    <- c("InsulinFasting","MatsudaIdx","Ghrelin","Amylin","GLP1")

## A) GE + hormones
pairs_chunk(data, vars = c(vars_outcome, vars_hormones),
            log10_vars = log_skewed, main = "GE + GI hormones")
tab_GE_horm <- spearman_table(data, vars = c(vars_outcome, vars_hormones),
                              log10_vars = log_skewed, sort_by = "abs_r")
print(tab_GE_horm, row.names = FALSE)
# write.csv(tab_GE_horm, "spearman_GE_hormones.csv", row.names = FALSE)

## B) Hormones + demographics
pairs_chunk(data, vars = c(vars_demo, vars_hormones),
            log10_vars = log_skewed, main = "Hormones + Demographics")
tab_demo <- spearman_table(data, vars = c(vars_demo, vars_hormones),
                           log10_vars = log_skewed, sort_by = "abs_r")

## C) Hormones + fasting glycemia
pairs_chunk(data, vars = c(vars_gly, vars_hormones),
            log10_vars = log_skewed, main = "Hormones + Fasting glycemia")
tab_gly <- spearman_table(data, vars = c(vars_gly, vars_hormones),
                          log10_vars = log_skewed, sort_by = "abs_r")

## D) Hormones + insulin dynamics
pairs_chunk(data, vars = c(vars_insulin, vars_hormones),
            log10_vars = log_skewed, main = "Hormones + Insulin dynamics")
tab_ins <- spearman_table(data, vars = c(vars_insulin, vars_hormones),
                          log10_vars = log_skewed, sort_by = "abs_r")

#-------MAYBE??--------
# ============================ #
#   Correlation (pairs + tab)  #
# ============================ #

data <- read.csv("/Users/claraazzano/Desktop/Data_T1.csv")

library(ggplot2)
library(rlang)
library(dplyr)
library(GGally)

out_dir <- "/Users/claraazzano/results_corr"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

vars_hormones <- c("Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY")
vars_demo     <- c("Age","Height","BW","BMI")
vars_gly      <- c("GlucoseFasting","HbA1c")
vars_insulin  <- c("InsulinFasting","MatsudaIdx","HOMAB")
vars_outcome  <- "GE"
log_skewed    <- c("InsulinFasting","MatsudaIdx","Ghrelin","Amylin","GLP1")

safe_log10 <- function(x) {
  eps <- 1e-8
  y <- x
  bad <- !is.finite(x) | x <= 0
  y[bad] <- NA_real_
  log10(y + eps)
}

make_pairs_df <- function(df, keep, log10_vars = character()) {
  keep <- intersect(keep, names(df))
  d <- df[keep]
  d <- d[vapply(d, is.numeric, logical(1L))]
  for (v in intersect(log10_vars, names(d))) d[[v]] <- safe_log10(d[[v]])
  if (ncol(d)) {
    sd_ok <- vapply(d, function(z) sd(z, na.rm = TRUE) > 0, logical(1L))
    d <- d[, sd_ok, drop = FALSE]
  }
  d
}

.panel.density <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  d <- density(x[is.finite(x)])
  if (length(d$x)) lines(d$x, d$y / max(d$y), lwd = 1.2)
}
.panel.points <- function(x, y, col = "black", pch = 19, ...) {
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) points(x[ok], y[ok], pch = pch,
                      col = adjustcolor(col, alpha.f = 0.35), cex = 0.6)
}
.panel.shade <- function(x, y, palette = c("#B71C1C","white","#1B5E20"), ...) {
  ok <- is.finite(x) & is.finite(y)
  r  <- if (sum(ok) >= 3) suppressWarnings(cor(x[ok], y[ok], method = "spearman")) else NA_real_
  col_bg <- colorRampPalette(palette)(200)
  idx <- if (is.finite(r)) max(1, min(200, round((r + 1)/2 * 199 + 1))) else 100
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = col_bg[idx], border = NA)
}

pairs_chunk <- function(df, vars, log10_vars = character(),
                        group = NULL, group_colors = c("#6A9AE2", "#E28A6A"),
                        main = NULL) {
  d <- make_pairs_df(df, keep = vars, log10_vars = log10_vars)
  if (!ncol(d)) stop("No numeric columns after filtering.")
  if (!is.null(group)) {
    g <- as.factor(group)
    cols <- adjustcolor(group_colors, alpha.f = 0.8)
    col_vec <- cols[as.integer(g)]
  } else col_vec <- "black"
  
  pairs(d,
        lower.panel = function(x, y, ...) .panel.points(x, y, col = col_vec),
        upper.panel = .panel.shade,
        diag.panel  = .panel.density,
        gap = 0.6,
        cex.labels = 1.0,
        main = main)
}

save_pairs <- function(filename, expr, width = 3200, height = 2400, res = 300) {
  png(filename, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

spearman_table <- function(df, vars, log10_vars = character(),
                           sort_by = c("abs_r","p","r")) {
  d <- make_pairs_df(df, keep = vars, log10_vars = log10_vars)
  cn <- colnames(d); k <- length(cn)
  out <- vector("list", max(0, k*(k-1)/2)); z <- 0L
  for (i in seq_len(k - 1)) for (j in (i + 1):k) {
    x <- d[[i]]; y <- d[[j]]
    ok <- is.finite(x) & is.finite(y)
    n  <- sum(ok)
    if (n >= 3) {
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      z <- z + 1L
      out[[z]] <- data.frame(var1 = cn[i], var2 = cn[j],
                             rho = unname(ct$estimate), p = unname(ct$p.value),
                             n = n, stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, out[seq_len(z)])
  if (is.null(tab) || !nrow(tab)) return(tab)
  tab$padj <- p.adjust(tab$p, method = "BH")
  tab$sig  <- tab$padj < 0.05
  sort_by <- match.arg(sort_by)
  ord <- switch(sort_by,
                abs_r = order(-abs(tab$rho), tab$padj, tab$p),
                p     = order(tab$p, -abs(tab$rho)),
                r     = order(-tab$rho))
  tab[ord, c("var1","var2","rho","p","padj","n","sig")]
}

# ========== A) GE + hormones ==========
save_pairs(file.path(out_dir, "pairs_GE_hormones.png"),
           pairs_chunk(data, vars = c(vars_outcome, vars_hormones),
                       log10_vars = log_skewed,
                       # group = data$Metformin,  # <- optional color
                       main = "GE + GI hormones")
)
# print to Plots pane
pairs_chunk(data, vars = c(vars_outcome, vars_hormones),
            log10_vars = log_skewed,
            # group = data$Metformin,
            main = "GE + GI hormones")
tab_GE_horm <- spearman_table(data, vars = c(vars_outcome, vars_hormones),
                              log10_vars = log_skewed, sort_by = "abs_r")
print(tab_GE_horm, row.names = FALSE)
write.csv(tab_GE_horm, file.path(out_dir, "spearman_GE_hormones.csv"), row.names = FALSE)

# ========== B) Hormones + demographics ==========
save_pairs(file.path(out_dir, "pairs_hormones_demographics.png"),
           pairs_chunk(data, vars = c(vars_demo, vars_hormones),
                       log10_vars = log_skewed,
                       main = "Hormones + Demographics")
)
pairs_chunk(data, vars = c(vars_demo, vars_hormones),
            log10_vars = log_skewed,
            main = "Hormones + Demographics")
tab_demo <- spearman_table(data, vars = c(vars_demo, vars_hormones),
                           log10_vars = log_skewed, sort_by = "abs_r")
write.csv(tab_demo, file.path(out_dir, "spearman_demo.csv"), row.names = FALSE)

# ========== C) Hormones + fasting glycemia ==========
save_pairs(file.path(out_dir, "pairs_hormones_glycemia.png"),
           pairs_chunk(data, vars = c(vars_gly, vars_hormones),
                       log10_vars = log_skewed,
                       main = "Hormones + Fasting glycemia")
)
pairs_chunk(data, vars = c(vars_gly, vars_hormones),
            log10_vars = log_skewed,
            main = "Hormones + Fasting glycemia")
tab_gly <- spearman_table(data, vars = c(vars_gly, vars_hormones),
                          log10_vars = log_skewed, sort_by = "abs_r")
write.csv(tab_gly, file.path(out_dir, "spearman_glycemia.csv"), row.names = FALSE)

# ========== D) Hormones + insulin dynamics ==========
save_pairs(file.path(out_dir, "pairs_hormones_insulin.png"),
           pairs_chunk(data, vars = c(vars_insulin, vars_hormones),
                       log10_vars = log_skewed,
                       main = "Hormones + Insulin dynamics")
)
pairs_chunk(data, vars = c(vars_insulin, vars_hormones),
            log10_vars = log_skewed,
            main = "Hormones + Insulin dynamics")
tab_ins <- spearman_table(data, vars = c(vars_insulin, vars_hormones),
                          log10_vars = log_skewed, sort_by = "abs_r")
write.csv(tab_ins, file.path(out_dir, "spearman_insulin.csv"), row.names = FALSE)

# ========== E) GE + Demographics ==========
save_pairs(file.path(out_dir, "pairs_GE_demographics.png"),
           pairs_chunk(data, vars = c(vars_outcome, vars_demo),
                       log10_vars = log_skewed,
                       # group = data$Metformin,  # optional color
                       main = "GE + Demographics")
)
pairs_chunk(data, vars = c(vars_outcome, vars_demo),
            log10_vars = log_skewed,
            # group = data$Metformin,
            main = "GE + Demographics")
tab_GE_demo <- spearman_table(data, vars = c(vars_outcome, vars_demo),
                              log10_vars = log_skewed, sort_by = "abs_r")
print(tab_GE_demo, row.names = FALSE)
write.csv(tab_GE_demo, file.path(out_dir, "spearman_GE_demographics.csv"), row.names = FALSE)

# ========== F) GE + Fasting glycemia ==========
save_pairs(file.path(out_dir, "pairs_GE_glycemia.png"),
           pairs_chunk(data, vars = c(vars_outcome, vars_gly),
                       log10_vars = log_skewed,
                       # group = data$Metformin,
                       main = "GE + Fasting glycemia")
)
pairs_chunk(data, vars = c(vars_outcome, vars_gly),
            log10_vars = log_skewed,
            # group = data$Metformin,
            main = "GE + Fasting glycemia")
tab_GE_gly <- spearman_table(data, vars = c(vars_outcome, vars_gly),
                             log10_vars = log_skewed, sort_by = "abs_r")
print(tab_GE_gly, row.names = FALSE)
write.csv(tab_GE_gly, file.path(out_dir, "spearman_GE_glycemia.csv"), row.names = FALSE)

# ========== G) GE + Insulin dynamics ==========
save_pairs(file.path(out_dir, "pairs_GE_insulin.png"),
           pairs_chunk(data, vars = c(vars_outcome, vars_insulin),
                       log10_vars = log_skewed,
                       # group = data$Metformin,
                       main = "GE + Insulin dynamics")
)
pairs_chunk(data, vars = c(vars_outcome, vars_insulin),
            log10_vars = log_skewed,
            # group = data$Metformin,
            main = "GE + Insulin dynamics")
tab_GE_ins <- spearman_table(data, vars = c(vars_outcome, vars_insulin),
                             log10_vars = log_skewed, sort_by = "abs_r")
print(tab_GE_ins, row.names = FALSE)
write.csv(tab_GE_ins, file.path(out_dir, "spearman_GE_insulin.csv"), row.names = FALSE)

# Files saved to: /Users/claraazzano/results_corr

library(dplyr)
library(corrplot)

# keep only numeric columns and drop any with zero variance (sd = 0 -> all same value)
num <- data %>%
  select(where(is.numeric)) %>%
  select(where(~ sd(., na.rm = TRUE) > 0))

# correlation matrix (Spearman is safer for these biomarkers; use = "pairwise.complete.obs" handles NAs)
M <- cor(num, method = "spearman", use = "pairwise.complete.obs")

# basic numeric upper-triangle plot
corrplot(M,
         method = "color",
         type   = "upper",
         diag   = FALSE,
         addCoef.col = "black",   # adds numbers on tiles
         number.digits = 2,
         number.cex    = 0.5,
         tl.cex        = 0.7,
         order         = "hclust")

#MDOEL TEST!!!!----------------------------------------------------                    
#======================================
                    #=============================


data <- read.csv("/Users/claraazzano/Desktop/Data_T1.csv")


library(ggplot2)
library(rlang)
library(dplyr)
library(GGally)
library(car)
library(lmtest)


#------- MODEL TIME------------------------------------------------------------

df <- data  # your dataset

## 0) Utilities ------------------------------------------------------------
safe_log <- function(x) {
  x <- suppressWarnings(as.numeric(x))        # coerce if factor/char
  ifelse(is.finite(x) & x > 0, log(x), NA_real_)
}

# Map of common aliases -> first that exists will be used
aliases <- list(
  GE             = c("GE","GE_t12","GE_t_half","GE_t1_2"),
  Age            = c("Age"),
  Sex            = c("Sex"),
  Height         = c("Height"),
  BW             = c("BW","Weight"),
  BMI            = c("BMI"),
  HbA1c          = c("HbA1c","HbA1C","A1c"),
  Metformin      = c("Metformin"),
  InsulinFasting = c("InsulinFasting","Insulin_Fasting","Insulin","FastingInsulin"),
  MatsudaIdx     = c("MatsudaIdx","Matsuda","Matsuda_Index"),
  HOMAB          = c("HOMAB","HOMA_B","HOMA.B","HOMA_Beta","HOMA%20B"),
  
  Gastrin        = c("Gastrin"),
  CCK            = c("CCK"),
  Ghrelin        = c("Ghrelin"),
  Amylin         = c("Amylin"),
  Glucagon       = c("Glucagon"),
  GLP1           = c("GLP1","GLP-1","GLP_1"),
  PYY            = c("PYY")
)

resolve_name <- function(df, key) {
  cand <- aliases[[key]]
  if (is.null(cand)) cand <- key
  hit <- intersect(cand, names(df))
  if (length(hit)) hit[1] else NA_character_
}

## 1) Resolve all needed columns ------------------------------------------
# choose insulin metric (TRUE = MatsudaIdx, FALSE = InsulinFasting)
use_matsuda <- TRUE

hormones_keys <- c("Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY")
covar_keys    <- c("Age","Sex","BMI","HbA1c","Metformin",
                   if (use_matsuda) "MatsudaIdx" else "InsulinFasting")

GE_col        <- resolve_name(df, "GE")
horm_cols     <- setNames(vapply(hormones_keys, resolve_name, "", df = df), hormones_keys)
covar_cols    <- setNames(vapply(covar_keys,    resolve_name, "", df = df), covar_keys)

# Report what was found / missing
message("\nResolved columns:")
print(list(GE = GE_col, hormones = horm_cols, covariates = covar_cols))

missing_h <- names(horm_cols)[horm_cols == ""]
missing_c <- names(covar_cols)[covar_cols == ""]
if (length(missing_h)) message("Missing hormone columns: ", paste(missing_h, collapse = ", "))
if (length(missing_c)) message("Missing covariate columns: ", paste(missing_c, collapse = ", "))

# Keep only present hormones
horm_present <- names(horm_cols[horm_cols != ""])
if (!length(horm_present)) stop("No hormone columns found. Check names(df).")

# Build working frame with resolved names
need_cols <- unique(c(GE_col, unname(horm_cols[horm_present]), unname(covar_cols[covar_cols != ""])))
dfW <- df[, need_cols, drop = FALSE]

## 2) Create log variables --------------------------------------------------
dfW$log_GE <- safe_log(dfW[[GE_col]])

# insulin metric (optional)
if (use_matsuda && covar_cols[["MatsudaIdx"]] != "") {
  dfW$log_Matsuda <- safe_log(dfW[[ covar_cols[["MatsudaIdx"]] ]])
} else if (!use_matsuda && covar_cols[["InsulinFasting"]] != "") {
  dfW$log_Insulin <- safe_log(dfW[[ covar_cols[["InsulinFasting"]] ]])
}

# log hormones
for (h in horm_present) {
  colname <- horm_cols[[h]]
  dfW[[paste0("log_", h)]] <- safe_log(dfW[[colname]])
}

## 3) Model helper ----------------------------------------------------------
fit_one_hormone <- function(h) {
  logH <- paste0("log_", h)
  
  # build RHS covariate list using LOGGED insulin metric if present
  insulin_term <- if (use_matsuda && "log_Matsuda" %in% names(dfW)) {
    "log_Matsuda"
  } else if (!use_matsuda && "log_Insulin" %in% names(dfW)) {
    "log_Insulin"
  } else {
    NULL
  }
  
  base_covars <- c(
    if (covar_cols[["Age"]]       != "") "Age"       else NULL,
    if (covar_cols[["Sex"]]       != "") "Sex"       else NULL,
    if (covar_cols[["BMI"]]       != "") "BMI"       else NULL,
    if (covar_cols[["HbA1c"]]     != "") "HbA1c"     else NULL,
    if (covar_cols[["Metformin"]] != "") "Metformin" else NULL,
    insulin_term
  )
  
  need <- c("log_GE", logH, base_covars)
  if (!all(need %in% names(dfW))) return(NULL)
  
  dd <- dfW[stats::complete.cases(dfW[, need, drop = FALSE]), need, drop = FALSE]
  if (nrow(dd) < 30) return(NULL)
  
  fml_full <- as.formula(paste0("log_GE ~ ", paste(c(logH, base_covars), collapse = " + ")))
  fml_base <- as.formula(paste0("log_GE ~ ", paste(base_covars, collapse = " + ")))
  
  m0 <- lm(fml_base, data = dd)
  m1 <- lm(fml_full, data = dd)
  
  s1  <- summary(m1); co <- coef(s1)
  if (!(logH %in% rownames(co))) return(NULL)
  
  beta <- unname(co[logH, "Estimate"])
  se   <- unname(co[logH, "Std. Error"])
  pval <- unname(co[logH, "Pr(>|t|)"])
  ci   <- beta + c(-1, 1) * qt(0.975, df = s1$df[2]) * se
  
  ln2  <- log(2)
  pct  <- (exp(beta * ln2) - 1) * 100
  pctL <- (exp(ci[1] * ln2) - 1) * 100
  pctU <- (exp(ci[2] * ln2) - 1) * 100
  
  RSS0 <- deviance(m0); RSS1 <- deviance(m1)
  pR2  <- if (is.finite(RSS0) && RSS0 > 0) (RSS0 - RSS1) / RSS0 else NA_real_
  
  data.frame(
    hormone = h,
    n = nrow(dd),
    pct_per_doubling = pct,
    pct_ci_lo = pctL,
    pct_ci_hi = pctU,
    p_value = pval,
    partial_R2 = pR2,
    adj_R2_full = unname(s1$adj.r.squared),
    beta = beta,
    beta_ci_lo = ci[1],
    beta_ci_hi = ci[2],
    stringsAsFactors = FALSE
  )
}

## 4) Run across hormones ---------------------------------------------------
res_list <- lapply(horm_present, fit_one_hormone)
res_tab  <- do.call(rbind, res_list)
res_tab  <- res_tab[order(res_tab$p_value), ]

res_tab_disp <- within(res_tab, {
  pct_per_doubling <- round(pct_per_doubling, 1)
  pct_ci_lo        <- round(pct_ci_lo, 1)
  pct_ci_hi        <- round(pct_ci_hi, 1)
  p_value          <- signif(p_value, 3)
  partial_R2       <- round(partial_R2, 3)
  adj_R2_full      <- round(adj_R2_full, 3)
})

print(res_tab_disp, row.names = FALSE)

cat("\nInterpretation:\n",
    "- pct_per_doubling = % change in GE t½ per doubling of the hormone (adjusted).\n",
    "- partial_R2 = extra variance explained beyond covariates.\n",
    "- adj_R2_full = overall adjusted R² of the full model.\n")

#-------FUN ADD ONS FOR THE MODEL---------------------------------------------

## ============================================================
##  Task 1 add-ons: Forest plot, Sensitivity, GLP-1×Metformin
##  (depends only on ggplot2; base R otherwise)
## ============================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# --- Utilities from the robust run (name resolver + logging) --------------

safe_log <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x) & x > 0, log(x), NA_real_)
}

aliases <- list(
  GE             = c("GE","GE_t12","GE_t_half","GE_t1_2"),
  Age            = c("Age"),
  Sex            = c("Sex"),
  Height         = c("Height"),
  BW             = c("BW","Weight"),
  BMI            = c("BMI"),
  HbA1c          = c("HbA1c","HbA1C","A1c"),
  Metformin      = c("Metformin","metformin"),
  InsulinFasting = c("InsulinFasting","Insulin_Fasting","Insulin","FastingInsulin"),
  MatsudaIdx     = c("MatsudaIdx","Matsuda","Matsuda_Index"),
  HOMAB          = c("HOMAB","HOMA_B","HOMA.B","HOMA_Beta"),
  
  Gastrin        = c("Gastrin"),
  CCK            = c("CCK"),
  Ghrelin        = c("Ghrelin"),
  Amylin         = c("Amylin"),
  Glucagon       = c("Glucagon"),
  GLP1           = c("GLP1","GLP-1","GLP_1"),
  PYY            = c("PYY")
)

resolve_name <- function(df, key) {
  cand <- aliases[[key]]; if (is.null(cand)) cand <- key
  hit <- intersect(cand, names(df))
  if (length(hit)) hit[1] else NA_character_
}

# --- Build a working frame with logs + chosen insulin covariate ------------

build_work_df <- function(df, use_matsuda = TRUE) {
  GE_col <- resolve_name(df, "GE")
  stopifnot(!is.na(GE_col))
  
  hormones_keys <- c("Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY")
  horm_cols <- setNames(vapply(hormones_keys, resolve_name, "", df = df), hormones_keys)
  horm_present <- names(horm_cols[horm_cols != ""])
  if (!length(horm_present)) stop("No hormone columns found.")
  
  covar_keys <- c("Age","Sex","BMI","HbA1c","Metformin",
                  if (use_matsuda) "MatsudaIdx" else "InsulinFasting")
  covar_cols <- setNames(vapply(covar_keys, resolve_name, "", df = df), covar_keys)
  
  need_cols <- unique(c(GE_col, unname(horm_cols[horm_present]), unname(covar_cols[covar_cols != ""])))
  d <- df[, need_cols, drop = FALSE]
  
  # create log outcome and log covariate (insulin metric)
  d$log_GE <- safe_log(d[[GE_col]])
  if (use_matsuda && covar_cols[["MatsudaIdx"]] != "") {
    d$log_Matsuda <- safe_log(d[[covar_cols[["MatsudaIdx"]]]])
  } else if (!use_matsuda && covar_cols[["InsulinFasting"]] != "") {
    d$log_Insulin <- safe_log(d[[covar_cols[["InsulinFasting"]]]])
  }
  
  # log hormones
  for (h in horm_present) d[[paste0("log_", h)]] <- safe_log(d[[horm_cols[[h]]]])
  
  list(
    d = d,
    GE_col = GE_col,
    horm_present = horm_present,
    covar_cols = covar_cols,
    use_matsuda = use_matsuda
  )
}

# --- Fit one adjusted model (same spec you used) ---------------------------

fit_one_hormone <- function(d, covar_cols, use_matsuda, hormone_name) {
  logH <- paste0("log_", hormone_name)
  
  insulin_term <- if (use_matsuda && "log_Matsuda" %in% names(d)) {
    "log_Matsuda"
  } else if (!use_matsuda && "log_Insulin" %in% names(d)) {
    "log_Insulin"
  } else NULL
  
  base_covars <- c(
    if (covar_cols[["Age"]]       != "") "Age"       else NULL,
    if (covar_cols[["Sex"]]       != "") "Sex"       else NULL,
    if (covar_cols[["BMI"]]       != "") "BMI"       else NULL,
    if (covar_cols[["HbA1c"]]     != "") "HbA1c"     else NULL,
    if (covar_cols[["Metformin"]] != "") "Metformin" else NULL,
    insulin_term
  )
  
  need <- c("log_GE", logH, base_covars)
  if (!all(need %in% names(d))) return(NULL)
  
  dd <- d[stats::complete.cases(d[, need, drop = FALSE]), need, drop = FALSE]
  if (nrow(dd) < 30) return(NULL)
  
  fml_full <- as.formula(paste0("log_GE ~ ", paste(c(logH, base_covars), collapse = " + ")))
  fml_base <- as.formula(paste0("log_GE ~ ", paste(base_covars, collapse = " + ")))
  
  m0 <- lm(fml_base, data = dd)
  m1 <- lm(fml_full, data = dd)
  
  s1  <- summary(m1); co <- coef(s1)
  if (!(logH %in% rownames(co))) return(NULL)
  
  beta <- unname(co[logH, "Estimate"])
  se   <- unname(co[logH, "Std. Error"])
  pval <- unname(co[logH, "Pr(>|t|)"])
  ci   <- beta + c(-1, 1) * qt(0.975, df = s1$df[2]) * se
  
  ln2  <- log(2)
  pct  <- (exp(beta * ln2) - 1) * 100
  pctL <- (exp(ci[1] * ln2) - 1) * 100
  pctU <- (exp(ci[2] * ln2) - 1) * 100
  
  RSS0 <- deviance(m0); RSS1 <- deviance(m1)
  pR2  <- if (is.finite(RSS0) && RSS0 > 0) (RSS0 - RSS1) / RSS0 else NA_real_
  
  list(
    row = data.frame(
      hormone = hormone_name,
      n = nrow(dd),
      pct_per_doubling = pct,
      pct_ci_lo = pctL,
      pct_ci_hi = pctU,
      p_value = pval,
      partial_R2 = pR2,
      adj_R2_full = unname(s1$adj.r.squared),
      beta = beta,
      beta_ci_lo = ci[1],
      beta_ci_hi = ci[2],
      stringsAsFactors = FALSE
    ),
    model = m1
  )
}

run_models <- function(df, use_matsuda = TRUE) {
  built <- build_work_df(df, use_matsuda = use_matsuda)
  d <- built$d
  out <- lapply(built$horm_present, function(h)
    fit_one_hormone(d, built$covar_cols, built$use_matsuda, h))
  rows <- do.call(rbind, lapply(out, `[[`, "row"))
  models <- setNames(lapply(out, `[[`, "model"), rows$hormone)
  
  rows_disp <- within(rows, {
    pct_per_doubling <- round(pct_per_doubling, 1)
    pct_ci_lo        <- round(pct_ci_lo, 1)
    pct_ci_hi        <- round(pct_ci_hi, 1)
    p_value          <- signif(p_value, 3)
    partial_R2       <- round(partial_R2, 3)
    adj_R2_full      <- round(adj_R2_full, 3)
  })
  rows_disp <- rows_disp[order(rows$p_value), ]
  list(table = rows_disp, table_raw = rows, models = models, built = built)
}

# --- Forest plot (uses the display table from run_models) ------------------

forest_plot <- function(res_table, title = "Adjusted association with GE (per doubling)") {
  dfp <- res_table
  dfp$hormone <- factor(dfp$hormone, levels = dfp$hormone[order(dfp$p_value)])  # order by p
  ggplot(dfp, aes(y = hormone, x = pct_per_doubling)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = pct_ci_lo, xmax = pct_ci_hi), height = 0.2) +
    geom_point(size = 2.5) +
    labs(
      title = title,
      x = "% change in GE t½ per doubling (95% CI)", y = NULL,
      caption = "Adjusted for Age, Sex, BMI, HbA1c, Metformin, and one insulin metric"
    ) +
    theme_minimal(base_size = 12)
}

# --- Interaction: GLP1 × Metformin ---------------------------------------

glp1_metformin_interaction <- function(df, use_matsuda = TRUE) {
  b <- build_work_df(df, use_matsuda = use_matsuda)
  d <- b$d
  # need GLP1 and Metformin present
  if (!("log_GLP1" %in% names(d))) stop("GLP1 not found.")
  if (!("Metformin" %in% names(d))) stop("Metformin not found.")
  # ensure Metformin is 0/1 numeric
  d$Metformin <- as.numeric(as.character(d$Metformin))
  if (!all(na.omit(d$Metformin) %in% c(0,1))) {
    # try to coerce factor labels "No/Yes"
    if (is.factor(d$Metformin)) d$Metformin <- as.numeric(d$Metformin) - 1
  }
  
  insulin_term <- if (use_matsuda && "log_Matsuda" %in% names(d)) "log_Matsuda"
  else if (!use_matsuda && "log_Insulin" %in% names(d)) "log_Insulin"
  else NULL
  
  base_covars <- c("Age","Sex","BMI","HbA1c","Metformin", insulin_term)
  need <- c("log_GE","log_GLP1","Metformin", base_covars)
  need <- unique(need[need != ""])
  dd <- d[stats::complete.cases(d[, need, drop=FALSE]), need, drop=FALSE]
  
  fml <- as.formula(paste0("log_GE ~ log_GLP1 * Metformin + ",
                           paste(setdiff(base_covars, "Metformin"), collapse = " + ")))
  m <- lm(fml, data = dd)
  s <- summary(m)
  
  # simple slopes at Metformin = 0 and = 1
  b1   <- coef(m)["log_GLP1"]
  bInt <- coef(m)["log_GLP1:Metformin"]
  V    <- vcov(m)
  ln2  <- log(2)
  
  slope0 <- b1
  var0   <- V["log_GLP1","log_GLP1"]
  slope1 <- b1 + bInt
  var1   <- V["log_GLP1","log_GLP1"] + V["log_GLP1:Metformin","log_GLP1:Metformin"] +
    2*V["log_GLP1","log_GLP1:Metformin"]
  
  ci_fun <- function(b, v) b + c(-1,1) * qt(0.975, df = s$df[2]) * sqrt(v)
  ci0 <- ci_fun(slope0, var0)
  ci1 <- ci_fun(slope1, var1)
  
  # convert to % change per doubling
  pct  <- function(b) (exp(b * ln2) - 1) * 100
  out <- data.frame(
    group = c("Metformin = 0 (No)","Metformin = 1 (Yes)"),
    pct_per_doubling = c(pct(slope0), pct(slope1)),
    pct_ci_lo        = c(pct(ci0[1]), pct(ci1[1])),
    pct_ci_hi        = c(pct(ci0[2]), pct(ci1[2])),
    stringsAsFactors = FALSE
  )
  list(model = m, summary = s, simple_slopes = out)
}

#-----Forest plot-------------------

primary <- run_models(data, use_matsuda = TRUE)   # using MatsudaIdx as before
primary$table                              # view table
forest_plot(primary$table,
            title = "Hormones vs GE — % change per doubling (adjusted)")   # draw forest

#------Sensitivity plot---------------

sens_ins <- run_models(data, use_matsuda = FALSE)  # now adjust for InsulinFasting instead

# Side-by-side comparison table
cmp <- merge(
  primary$table[, c("hormone","pct_per_doubling","pct_ci_lo","pct_ci_hi","p_value","partial_R2")],
  sens_ins$table[, c("hormone","pct_per_doubling","pct_ci_lo","pct_ci_hi","p_value","partial_R2")],
  by = "hormone", suffixes = c("_Matsuda","_Insulin")
)
print(cmp[order(cmp$p_value_Matsuda), ], row.names = FALSE)

#--------Interaction test--------------

ix <- glp1_metformin_interaction(data, use_matsuda = TRUE)
ix$simple_slopes   # %Δ GE per doubling GLP-1, separately in Metformin No vs Yes
summary(ix$model)  # regression table (look at the log_GLP1:Metformin term)

## ============================================================
## Plot GLP-1 × Metformin interaction
## ============================================================

plot_glp1_metformin <- function(df, use_matsuda = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  library(ggplot2)
  
  # rebuild dataset
  b <- build_work_df(df, use_matsuda = use_matsuda)
  d <- b$d
  if (!("log_GLP1" %in% names(d))) stop("Need log_GLP1")
  if (!("Metformin" %in% names(d))) stop("Need Metformin")
  
  # ensure Metformin is binary numeric
  d$Metformin <- as.numeric(as.character(d$Metformin))
  if (!all(na.omit(d$Metformin) %in% c(0,1))) {
    if (is.factor(d$Metformin)) d$Metformin <- as.numeric(d$Metformin) - 1
  }
  d$Metformin <- factor(d$Metformin, labels = c("No Metformin", "Yes Metformin"))
  
  # fit model with interaction
  insulin_term <- if (use_matsuda && "log_Matsuda" %in% names(d)) "log_Matsuda"
  else if (!use_matsuda && "log_Insulin" %in% names(d)) "log_Insulin"
  else NULL
  base_covars <- c("Age","Sex","BMI","HbA1c","Metformin", insulin_term)
  need <- unique(c("log_GE","log_GLP1","Metformin", base_covars))
  dd <- d[stats::complete.cases(d[, need, drop=FALSE]), need, drop=FALSE]
  
  fml <- as.formula(paste0("log_GE ~ log_GLP1 * Metformin + ",
                           paste(setdiff(base_covars, "Metformin"), collapse = " + ")))
  m <- lm(fml, data = dd)
  
  # plot
  ggplot(dd, aes(x = log_GLP1, y = log_GE, color = Metformin)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
    labs(
      title = "GLP-1 × Metformin interaction on GE half-life",
      x = "log(GLP-1)",
      y = "log(GE half-life)"
    ) +
    scale_color_manual(values = c("#2C7BB6", "#D7191C")) +
    theme_minimal(base_size = 13)
}

plot_glp1_metformin(data, use_matsuda = TRUE)   # use Matsuda adjustment

#---------ASSESS MODEL TIME--------------------------------------------

## ============================================================
## MODEL DIAGNOSTICS — base R + car + lmtest
## ============================================================

# install once if needed
# install.packages(c("car","lmtest"))



# ---- helper: refit without using `update()` (avoids env/formula issues) ---
refit_without <- function(mod, drop_idx = integer(0)) {
  stopifnot(inherits(mod, "lm"))
  mf <- model.frame(mod)                                     # data used in 'mod'
  if (length(drop_idx)) mf <- mf[-drop_idx, , drop = FALSE]  # drop rows
  lm(formula(mod), data = mf)                                # refit from stored formula + data
}

# ---- main diagnostics function -------------------------------------------
diag_report <- function(mod, main_predictor) {
  stopifnot(inherits(mod, "lm"))
  
  cat("\n=== MODEL DIAGNOSTICS ===\n")
  cat("Formula: ", deparse(formula(mod)), "\n\n")
  
  # 1) classic 4 plots
  op <- par(mfrow = c(2,2)); on.exit(par(op), add = TRUE)
  plot(mod)  # Residuals vs Fitted, QQ plot, Scale-Location, Residuals vs Leverage
  
  # 2) numbers to report
  cat("\n-- Summary --\n"); print(summary(mod))
  cat("\n-- VIF --\n");        print(vif(mod))                 # multicollinearity
  cat("\n-- Shapiro-Wilk (residuals) --\n"); print(shapiro.test(residuals(mod)))
  cat("\n-- Breusch–Pagan (homoscedasticity) --\n"); print(bptest(mod))
  
  # 3) influence (Cook's distance)
  cooks <- cooks.distance(mod)
  thr <- 4/length(cooks)                                      # common rule-of-thumb
  par(mfrow = c(1,1))
  plot(cooks, type = "h", main = "Cook's distance"); abline(h = thr, col = "red", lty = 2)
  cat("\n-- Influence --\n")
  cat("Threshold 4/n =", signif(thr, 3), "\n")
  cat("Top 5 Cook's:\n"); print(head(sort(cooks, decreasing = TRUE), 5))
  idx_infl <- which(cooks > thr)
  cat("N influential (Cook's > 4/n):", length(idx_infl), "\n")
  
  # 4) robustness refit (drop influential rows) and compare main effect
  if (length(idx_infl)) {
    m_refit <- refit_without(mod, drop_idx = idx_infl)
    s0 <- summary(mod); s1 <- summary(m_refit)
    
    if (main_predictor %in% rownames(coef(s0))) {
      ln2 <- log(2)
      b0 <- coef(s0)[main_predictor, "Estimate"]; se0 <- coef(s0)[main_predictor, "Std. Error"]
      b1 <- coef(s1)[main_predictor, "Estimate"]; se1 <- coef(s1)[main_predictor, "Std. Error"]
      pct0 <- (exp(b0 * ln2) - 1) * 100
      pct1 <- (exp(b1 * ln2) - 1) * 100
      
      cat(sprintf("\nEffect of %s (%% per doubling): full = %.1f%%, refit = %.1f%%\n",
                  main_predictor, pct0, pct1))
      cat(sprintf("SE change: full = %.4f, refit = %.4f\n", se0, se1))
    } else {
      cat("\n(main_predictor not found in coefficient table; skipping effect comparison.)\n")
    }
  } else {
    cat("\nNo observations above the 4/n Cook’s threshold.\n")
  }
}

# ---- quick convenience: save 4-panel plot to PNG -------------------------
save_diag_png <- function(mod, filename = "diagnostics.png", width = 1600, height = 1200, res = 150) {
  png(filename, width = width, height = height, res = res)
  par(mfrow = c(2,2))
  plot(mod)
  dev.off()
  cat("Saved:", filename, "\n")
}

## ============================================================
## HOW TO RUN (examples)
## ============================================================

# assuming you already have: primary <- run_models(data, use_matsuda = TRUE)

# 1) GLP-1 model diagnostics
diag_report(primary$models$GLP1,    main_predictor = "log_GLP1")
# save_diag_png(primary$models$GLP1, "diag_glp1.png")

# 2) PYY model diagnostics
diag_report(primary$models$PYY,     main_predictor = "log_PYY")
# save_diag_png(primary$models$PYY, "diag_pyy.png")

# 3) Ghrelin model diagnostics
diag_report(primary$models$Ghrelin, main_predictor = "log_Ghrelin")
# save_diag_png(primary$models$Ghrelin, "diag_ghrelin.png")

# 4) quick fit comparators
AIC(primary$models$GLP1, primary$models$PYY, primary$models$Ghrelin)
sapply(primary$models[c("GLP1","PYY","Ghrelin")], function(m) summary(m)$adj.r.squared)


#SAVE EVERYTHING

# ============================================================
# EXPORT: results + plots for the regression models
# ============================================================

# 0) output folders ----------------------------------------------------------
out_dir <- "/Users/claraazzano/plots_model"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "diagnostics"), showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"),       showWarnings = FALSE)

# 1) save result tables ------------------------------------------------------
# main adjusted table
write.csv(primary$table,
          file.path(out_dir, "tables", "model_primary_Matsuda.csv"),
          row.names = FALSE)

# sensitivity (InsulinFasting) table
write.csv(sens_ins$table,
          file.path(out_dir, "tables", "model_sensitivity_InsulinFasting.csv"),
          row.names = FALSE)

# side-by-side comparison
cmp <- merge(
  primary$table[, c("hormone","pct_per_doubling","pct_ci_lo","pct_ci_hi","p_value","partial_R2")],
  sens_ins$table[, c("hormone","pct_per_doubling","pct_ci_lo","pct_ci_hi","p_value","partial_R2")],
  by = "hormone", suffixes = c("_Matsuda","_Insulin")
)
write.csv(cmp, file.path(out_dir, "tables", "comparison_Matsuda_vs_Insulin.csv"),
          row.names = FALSE)

# GLP-1 × Metformin simple slopes
write.csv(ix$simple_slopes,
          file.path(out_dir, "tables", "interaction_GLP1xMetformin_simple_slopes.csv"),
          row.names = FALSE)

# 2) forest plot (fix deprecation of geom_errorbarh) ------------------------
forest_plot <- function(res_table, title = "Hormones vs GE — % change per doubling (adjusted)") {
  dfp <- res_table
  dfp$hormone <- factor(dfp$hormone, levels = dfp$hormone[order(dfp$p_value)])
  ggplot(dfp, aes(y = hormone, x = pct_per_doubling)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    # use horizontal error bars via orientation = "y"
    geom_errorbar(aes(xmin = pct_ci_lo, xmax = pct_ci_hi), width = 0.2, orientation = "y") +
    geom_point(size = 2.5) +
    labs(
      title = title,
      x = "% change in GE t½ per doubling (95% CI)", y = NULL,
      caption = "Adjusted for Age, Sex, BMI, HbA1c, Metformin, and one insulin metric"
    ) +
    theme_minimal(base_size = 12)
}

p_forest <- forest_plot(primary$table)
ggsave(file.path(out_dir, "forest_primary.png"), p_forest, width = 7.5, height = 6, dpi = 300)

# 3) GLP-1 × Metformin interaction plot -------------------------------------
p_inter <- plot_glp1_metformin(data, use_matsuda = TRUE)
ggsave(file.path(out_dir, "interaction_GLP1xMetformin.png"), p_inter, width = 7.5, height = 6, dpi = 300)

# 4) diagnostics: residual panels + Cook’s distance -------------------------
# helper from your script (re-use)
save_diag_png <- function(mod, filename = "diagnostics.png", width = 1600, height = 1200, res = 150) {
  png(filename, width = width, height = height, res = res)
  par(mfrow = c(2,2)); plot(mod); dev.off()
  message("Saved: ", filename)
}
save_cooks_png <- function(mod, filename) {
  png(filename, width = 1200, height = 900, res = 150)
  cooks <- cooks.distance(mod); thr <- 4/length(cooks)
  plot(cooks, type = "h", main = "Cook's distance"); abline(h = thr, col = "red", lty = 2)
  dev.off()
  message("Saved: ", filename)
}

mods_to_save <- c("GLP1","PYY","Ghrelin")
for (nm in mods_to_save) {
  m <- primary$models[[nm]]
  stopifnot(inherits(m, "lm"))
  save_diag_png(m, file.path(out_dir, "diagnostics", paste0("diag_", nm, ".png")))
  save_cooks_png(m, file.path(out_dir, "diagnostics", paste0("cooks_", nm, ".png")))
}

# 5) also print key results to the console ----------------------------------
cat("\n--- Primary adjusted results (Matsuda) ---\n"); print(primary$table)
cat("\n--- Sensitivity (InsulinFasting) ---\n"); print(sens_ins$table)
cat("\n--- GLP-1 × Metformin simple slopes (%Δ per doubling) ---\n"); print(ix$simple_slopes)
cat("\nFiles written to: ", out_dir, "\n")

#T TESTS AND MORE STATS ANALYSIS ETC ------------------------------------------
       #===========================================
       #=================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(broom)


data <- read.csv("/Users/claraazzano/Desktop/Data_T1.csv")

# Variables
vars_hormones <- c("Gastrin","CCK","Ghrelin","Amylin","Glucagon","GLP1","PYY")
vars_outcome  <- "GE"
vars_all      <- c(vars_outcome, vars_hormones)

# Output folders
root_dir   <- "/Users/claraazzano/plots"
box_dir    <- file.path(root_dir, "boxplots")
table_dir  <- file.path(root_dir, "tables")
dirs_to_make <- c(root_dir, box_dir, table_dir,
                  file.path(box_dir, "Metformin"),
                  file.path(box_dir, "DiabetesComplications"))
invisible(lapply(dirs_to_make, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))

# Helper to keep binary 0/1 groups and drop NA/other levels
.clean_binary <- function(df, group_var, yvar) {
  df %>%
    filter(!is.na(.data[[group_var]]), !is.na(.data[[yvar]])) %>%
    mutate(!!group_var := as.integer(.data[[group_var]])) %>%
    filter(.data[[group_var]] %in% c(0L, 1L)) %>%
    mutate(!!group_var := factor(.data[[group_var]], levels = c(0,1),
                                 labels = c("0","1")))
}

# ============================================================
# 2) Welch t-tests (returns a tibble and writes CSV)
# ============================================================
run_t_tests <- function(df, group_var, variables, table_dir) {
  res <- bind_rows(lapply(variables, function(v) {
    tmp <- .clean_binary(df, group_var, v)
    if (length(unique(tmp[[group_var]])) < 2) return(NULL)
    tt  <- t.test(as.formula(paste(v, "~", group_var)), data = tmp)
    tibble(
      variable = v,
      group    = group_var,
      p_value  = tt$p.value,
      mean_0   = mean(tmp[[v]][tmp[[group_var]] == "0"], na.rm = TRUE),
      mean_1   = mean(tmp[[v]][tmp[[group_var]] == "1"], na.rm = TRUE),
      n_0      = sum(tmp[[group_var]] == "0"),
      n_1      = sum(tmp[[group_var]] == "1")
    )
  }))
  if (nrow(res)) {
    out_csv <- file.path(table_dir, paste0("ttests_", group_var, ".csv"))
    write.csv(res, out_csv, row.names = FALSE)
  }
  res
}

# ============================================================
# 3) Boxplots (save AND show)
# ============================================================
plot_box_save <- function(df, group_var, variables, box_dir, show = TRUE) {
  subdir <- file.path(box_dir, group_var)
  if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
  
  for (v in variables) {
    tmp <- .clean_binary(df, group_var, v)
    if (nrow(tmp) == 0 || length(unique(tmp[[group_var]])) < 2) next
    
    p <- ggplot(tmp, aes(x = .data[[group_var]], y = .data[[v]], fill = .data[[group_var]])) +
      geom_boxplot(alpha = 0.65, outlier.shape = NA) +
      geom_jitter(width = 0.18, alpha = 0.55, size = 1, color = "grey20") +
      stat_compare_means(method = "t.test",
                         label = "p.format", label.y.npc = "top") +
      scale_fill_manual(values = c("#F28E8E", "#7BD1D1")) +
      labs(title = paste(v, "by", group_var),
           x = group_var, y = v) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none")
    
    # Save
    ggsave(filename = file.path(subdir, paste0("box_", v, "_by_", group_var, ".png")),
           plot = p, width = 10, height = 6, dpi = 300)
    
    # Show
    if (isTRUE(show)) print(p)
  }
  invisible(NULL)
}

# ============================================================
# 4) Run for Metformin and DiabetesComplications
# ============================================================

## A) Metformin
cat("\n--- Welch t-tests (Metformin) ---\n")
tt_metformin <- run_t_tests(data, "Metformin", vars_all, table_dir)
print(tt_metformin)
plot_box_save(data, "Metformin", vars_all, box_dir, show = TRUE)

## B) DiabetesComplications
cat("\n--- Welch t-tests (DiabetesComplications) ---\n")
tt_complic <- run_t_tests(data, "DiabetesComplications", vars_all, table_dir)
print(tt_complic)
plot_box_save(data, "DiabetesComplications", vars_all, box_dir, show = TRUE)

cat("\nSaved boxplots to:\n",
    file.path(box_dir, "Metformin"), "\n",
    file.path(box_dir, "DiabetesComplications"), "\n",
    "T-test tables saved to:\n", table_dir, "\n", sep = "")


# DEMOGRAHICS PART------------------------------------

# ============================================================
# 1. Setup
# ============================================================


# Demographic predictors
vars_demo <- c("Age","Sex","Height","BW","BMI")
outcome   <- "DiabetesComplications"

# Output dirs
root_dir   <- "/Users/claraazzano/plots"
box_dir    <- file.path(root_dir, "boxplots_demo")
table_dir  <- file.path(root_dir, "tables_demo")
dirs_to_make <- c(root_dir, box_dir, table_dir)
invisible(lapply(dirs_to_make, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))

# ============================================================
# 2. Logistic regression helper
# ============================================================
log_reg_demo <- function(data, outcome, predictors) {
  form <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  fit  <- glm(form, data = data, family = binomial)
  est  <- summary(fit)$coefficients
  
  out  <- data.frame(term = rownames(est),
                     estimate = est[,1],
                     se = est[,2],
                     z = est[,3],
                     p = est[,4])
  out$OR   <- exp(out$estimate)
  out$LCL  <- exp(out$estimate - 1.96*out$se)
  out$UCL  <- exp(out$estimate + 1.96*out$se)
  rownames(out) <- NULL
  out
}

# ============================================================
# 3. Run regressions
# ============================================================

# Univariate
uni_results <- lapply(vars_demo, function(v) log_reg_demo(data, outcome, v))
names(uni_results) <- vars_demo

# Save univariate results
for (nm in names(uni_results)) {
  cat("\n--- DiabetesComplications ~", nm, "---\n")
  print(uni_results[[nm]])
  write.csv(uni_results[[nm]], file.path(table_dir, paste0("logreg_uni_", nm, ".csv")),
            row.names = FALSE)
}

# Multivariable
multi_results <- log_reg_demo(data, outcome, vars_demo)
cat("\n--- Multivariable logistic regression ---\n")
print(multi_results)
write.csv(multi_results, file.path(table_dir, "logreg_multivariable.csv"),
          row.names = FALSE)

# ============================================================
# 4. Boxplots (continuous) + barplot (categorical)
# ============================================================
plot_box_demo <- function(df, group_var, variables, out_dir, show = TRUE) {
  for (v in variables) {
    if (v == "Sex") {
      p <- ggplot(df, aes(x = factor(.data[[v]]), fill = factor(.data[[group_var]]))) +
        geom_bar(position = "fill") +
        scale_y_continuous(labels = scales::percent) +
        labs(title = paste(v, "by", group_var),
             x = v, y = "Proportion with complications") +
        theme_minimal(base_size = 13)
    } else {
      p <- ggplot(df, aes(x = factor(.data[[group_var]]),
                          y = .data[[v]],
                          fill = factor(.data[[group_var]]))) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
        stat_compare_means(method = "t.test") +
        labs(title = paste(v, "by", group_var),
             x = group_var, y = v) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "none")
    }
    
    # Save
    fname <- paste0("box_", v, "_by_", group_var, ".png")
    ggsave(file.path(out_dir, fname), plot = p,
           width = 10, height = 6, dpi = 300)
    
    # Show
    if (isTRUE(show)) print(p)
  }
}

# Run boxplots
plot_box_demo(data, outcome, vars_demo, box_dir, show = TRUE)

cat("\nBoxplots saved to:", box_dir,
    "\nTables saved to:", table_dir, "\n")


#DEMOGRAPHICS-----------------------------------

# Demographic predictors
vars_demo <- c("Age","Sex","Height","BW","BMI")
outcome   <- "DiabetesComplications"

# Output dirs
root_dir   <- "/Users/claraazzano/plots"
box_dir    <- file.path(root_dir, "boxplots_demo")
table_dir  <- file.path(root_dir, "tables_demo")
dirs_to_make <- c(root_dir, box_dir, table_dir)
invisible(lapply(dirs_to_make, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))

# ============================================================
# 2. Logistic regression helper
# ============================================================
log_reg_demo <- function(data, outcome, predictors) {
  form <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  fit  <- glm(form, data = data, family = binomial)
  est  <- summary(fit)$coefficients
  
  out  <- data.frame(term = rownames(est),
                     estimate = est[,1],
                     se = est[,2],
                     z = est[,3],
                     p = est[,4])
  out$OR   <- exp(out$estimate)
  out$LCL  <- exp(out$estimate - 1.96*out$se)
  out$UCL  <- exp(out$estimate + 1.96*out$se)
  rownames(out) <- NULL
  out
}

# ============================================================
# 3. Run regressions
# ============================================================

# Univariate
uni_results <- lapply(vars_demo, function(v) log_reg_demo(data, outcome, v))
names(uni_results) <- vars_demo

# Save univariate results
for (nm in names(uni_results)) {
  cat("\n--- DiabetesComplications ~", nm, "---\n")
  print(uni_results[[nm]])
  write.csv(uni_results[[nm]], file.path(table_dir, paste0("logreg_uni_", nm, ".csv")),
            row.names = FALSE)
}

# Multivariable
multi_results <- log_reg_demo(data, outcome, vars_demo)
cat("\n--- Multivariable logistic regression ---\n")
print(multi_results)
write.csv(multi_results, file.path(table_dir, "logreg_multivariable.csv"),
          row.names = FALSE)

# ============================================================
# 4. Boxplots (continuous) + barplot (categorical)
# ============================================================
plot_box_demo <- function(df, group_var, variables, out_dir, show = TRUE) {
  for (v in variables) {
    if (v == "Sex") {
      p <- ggplot(df, aes(x = factor(.data[[v]]), fill = factor(.data[[group_var]]))) +
        geom_bar(position = "fill") +
        scale_y_continuous(labels = scales::percent) +
        labs(title = paste(v, "by", group_var),
             x = v, y = "Proportion with complications") +
        theme_minimal(base_size = 13)
    } else {
      p <- ggplot(df, aes(x = factor(.data[[group_var]]),
                          y = .data[[v]],
                          fill = factor(.data[[group_var]]))) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
        stat_compare_means(method = "t.test") +
        labs(title = paste(v, "by", group_var),
             x = group_var, y = v) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "none")
    }
    
    # Save
    fname <- paste0("box_", v, "_by_", group_var, ".png")
    ggsave(file.path(out_dir, fname), plot = p,
           width = 10, height = 6, dpi = 300)
    
    # Show
    if (isTRUE(show)) print(p)
  }
}

# Run boxplots
plot_box_demo(data, outcome, vars_demo, box_dir, show = TRUE)

cat("\nBoxplots saved to:", box_dir,
    "\nTables saved to:", table_dir, "\n")
