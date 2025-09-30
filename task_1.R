
THIS IS FOR VIZUALIZATIONS

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


