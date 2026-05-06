
## ------------------------------------------------------------------------------------
## -------------------------------------- LIBRARIES -----------------------------------
## ------------------------------------------------------------------------------------
library(baseballr)
library(readr)
library(dplyr)
library(tidyverse)
library(depmixS4)
library(zoo)
library(forcats)
library(moments)
library(DescTools)
library(ggplot2)
library(tidyr)
library(ggtext)
library(knitr)

## ------------------------------------------------------------------------------------
## -------------------------------------- DATA -----------------------------------
## ------------------------------------------------------------------------------------


seasons = c(2022, 2023, 2024, 2025)  

fetch_season_data = function(season) {
  url = paste0(
    'https://baseballsavant.mlb.com/statcast_search/csv?all=true',
    '&hfGT=R%7C&hfSea=', season, '%7C&player_type=batter',
    '&batters_lookup%5B%5D=677594&min_pitches=0&min_results=0',
    '&group_by=name&sort_col=pitches&player_event_sort=h_launch_speed',
    '&sort_order=desc&min_abs=0&type=details'
  )
  
  message("Fetching season: ", season)
  
  tryCatch({
    # Use read_csv and define the na strings explicitly
    df = read_csv(url, na = c("", "NA", "null"), show_col_types = FALSE)
    df$season = season
    return(df)
  }, error = function(e) {
    warning("Failed to fetch season ", season, ": ", e$message)
    return(NULL)
  })
}

all_seasons_list = lapply(seasons, fetch_season_data)
data_raw = bind_rows(Filter(Negate(is.null), all_seasons_list))


jrod_data = data_raw %>% 
  dplyr::select(
    game_date,
    season,
    events,
    launch_speed,
    launch_angle,
    estimated_ba_using_speedangle,
    estimated_woba_using_speedangle,
    bb_type
  )

jrod_data = jrod_data %>% 
  # Filter out mid-at-bat pitches, but KEEP all plate appearance enders
  filter(!is.na(events)) %>% 
  mutate(
    game_date = as.Date(game_date),
    events = as.factor(events),
    bb_type = as.factor(bb_type),
    
    # Create the continuous sequence metric
    xwoba_continuous = case_when(
      # 1. If he hit the ball, use the Statcast speed/angle estimate
      !is.na(estimated_woba_using_speedangle) ~ estimated_woba_using_speedangle,
      
      # 2. Strikeouts are worth zero
      events %in% c("strikeout", "strikeout_double_play") ~ 0.00,
      
      # 3. Walks are positive outcomes
      events %in% c("walk", "intentional_walk") ~ 0.69,
      
      # 4. Hit-by-pitch is slightly more valuable than a walk in wOBA
      events == "hit_by_pitch" ~ 0.72,
      
      # 5. Catch-all for rare non-contact outs (e.g., caught stealing, pickoffs, etc.)
      TRUE ~ 0.00 
    )


  )
# At bat level data with PA index and season day

jrod_ab = jrod_data %>%
  arrange(game_date) %>%
  group_by(season, game_date) %>%
  mutate(pa_index = row_number()) %>%
  ungroup() %>%
  arrange(season, game_date, pa_index) %>%
  mutate(
    pa_id = row_number(),  
    season_day = as.numeric(game_date - as.Date(paste0(season, "-04-01")))
  )

# ------------------------------------------------------------------------------------
# -------------------------------------- ROLLING HMM FIT -----------------------------
# ------------------------------------------------------------------------------------

# Calculate a 25-PA rolling average of xwOBA
jrod_ab = jrod_ab %>%
  group_by(season) %>%
  arrange(pa_id) %>%
  mutate(
    xwoba_roll25 = rollmean(xwoba_continuous, k = 25, fill = NA, align = "right")
  ) %>%
  ungroup()

# HMMs cannot handle NA values, which will appear in the first 24 PAs of each season 
# due to the rolling window. We will filter them out for the model.
jrod_hmm_data = jrod_ab %>% filter(!is.na(xwoba_roll25))

# Recalculate plate appearances per season for the trimmed dataset
pa_counts_roll = jrod_hmm_data %>%
  group_by(season) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(n)

## ------------------------------------------------------------------------------------
## -------------------------------------- PRE HMM DIAGNOSTICS -----------------------------------
## ------------------------------------------------------------------------------------

par(mfrow = c(1, 3))
hist(jrod_ab$xwoba_continuous,
     breaks = 40, main = "xwOBA raw", xlab = "xwOBA",
     col = "steelblue", border = "white")
abline(v = median(jrod_ab$xwoba_continuous, na.rm = TRUE), col = "green", lwd = 2)
abline(v = mean(jrod_ab$xwoba_continuous, na.rm = TRUE),   col = "orange", lwd = 2)
legend("topright", c("Median", "Mean"), col = c("green", "orange"), lwd = 2, cex = 0.7)

# Square-root transform (preserves zeros, reduces right skew)
hist(sqrt(jrod_ab$xwoba_continuous),
     breaks = 40, main = "xwOBA sqrt-transformed", xlab = "sqrt(xwOBA)",
     col = "steelblue", border = "white")

# QQ plot of raw
qqnorm(jrod_ab$xwoba_continuous, main = "QQ plot — raw xwOBA")
qqline(jrod_ab$xwoba_continuous, col = "red")

par(mfrow = c(1, 1))

# 2. Numeric skew stats
cat("Skewness (raw):  ", skewness(jrod_ab$xwoba_continuous, na.rm = TRUE), "\n")
cat("Skewness (sqrt): ", skewness(sqrt(jrod_ab$xwoba_continuous), na.rm = TRUE), "\n")
cat("Kurtosis (raw):  ", kurtosis(jrod_ab$xwoba_continuous, na.rm = TRUE), "\n")

# 3. Outlier count (xwOBA > 1.5 are essentially all HRs)
cat("\nPAs with xwOBA > 1.5:", 
    sum(jrod_ab$xwoba_continuous > 1.5, na.rm = TRUE), 
    "(", round(mean(jrod_ab$xwoba_continuous > 1.5, na.rm = TRUE)*100, 1), "% of PAs)\n")

# 4. Zero-xwOBA PAs (outs with no contact value)  
cat("PAs with xwOBA = 0:", 
    sum(jrod_ab$xwoba_continuous == 0, na.rm = TRUE), "\n")

par(mfrow = c(1, 2))

hist(jrod_hmm_data$xwoba_roll25,
     breaks = 40, main = "25-PA Rolling xwOBA", xlab = "xwOBA",
     col = "steelblue", border = "white")

qqnorm(jrod_hmm_data$xwoba_roll25, main = "QQ plot — Rolling xwOBA")
qqline(jrod_hmm_data$xwoba_roll25, col = "red")

par(mfrow = c(1, 1))

cat("Skewness:", skewness(jrod_hmm_data$xwoba_roll25, na.rm = TRUE), "\n")
cat("Kurtosis:", kurtosis(jrod_hmm_data$xwoba_roll25, na.rm = TRUE), "\n")

shapiro.test(sample(jrod_hmm_data$xwoba_roll25, 500))

# ------------------------------------------------------------------------------------
# -------------------------------------- MODEL SELECTION (BIC/AIC) -------------------
# ------------------------------------------------------------------------------------

n_states_to_try = 2:5

model_selection_results = map_dfr(n_states_to_try, function(k) {
  hmm_k = depmix(xwoba_roll25 ~ 1, 
                 data    = jrod_hmm_data, 
                 nstates = k, 
                 family  = gaussian(), 
                 ntimes  = pa_counts_roll)
  
  fit_k = tryCatch(
    fit(hmm_k, verbose = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit_k)) {
    return(tibble(states = k, logLik = NA, AIC = NA, BIC = NA))
  }
  
  tibble(
    states = k,
    logLik = as.numeric(logLik(fit_k)),
    AIC    = AIC(fit_k),
    BIC    = BIC(fit_k)
  )
})

print(model_selection_results)

# Plot it
model_selection_results %>%
  pivot_longer(cols = c(AIC, BIC), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = states, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "gray40") +
  annotate("text", x = 3.1, y = max(model_selection_results$AIC, na.rm = TRUE), 
           label = "Your choice", hjust = 0, color = "gray40", size = 3.5) +
  scale_color_manual(values = c("AIC" = "#de2d26", "BIC" = "#3182bd")) +
  theme_minimal() +
  labs(
    title = "HMM Model Selection: AIC & BIC by Number of States",
    x     = "Number of Hidden States",
    y     = "Information Criterion",
    color = "Metric"
  ) +
  theme(plot.title = element_text(face = "bold", size = 14))


## ------------------------------------------------------------------------------------
## -------------------------------------- HMM FIT -----------------------------------
## ------------------------------------------------------------------------------------

# Fit a 3-state HMM on the rolling average
hmm_roll = depmix(xwoba_roll25 ~ 1, data = jrod_hmm_data, nstates = 3, family = gaussian(), ntimes = pa_counts_roll)
fit_roll = fit(hmm_roll, verbose = FALSE)

# Decode the states
jrod_hmm_data$state_roll = viterbi(fit_roll)$state


# ------------------------------------------------------------------------------------
# -------------------------------------- DYNAMIC STATE MAPPING -----------------------
# ------------------------------------------------------------------------------------

# 1. Find the mean rolling xwOBA for each state and assign labels dynamically
state_mapping = jrod_hmm_data %>%
  group_by(state_roll) %>%
  summarize(mean_val = mean(xwoba_roll25, na.rm = TRUE)) %>%
  arrange(mean_val) %>%
  mutate(state_name = c("Cold", "Average", "Hot")) 

# 2. Delete the old state_name column from the dataset to prevent a collision
jrod_hmm_data$state_name = NULL

# 3. Join the new labels back to the main dataset
jrod_hmm_data = jrod_hmm_data %>%
  left_join(state_mapping %>% dplyr::select(state_roll, state_name), by = "state_roll") %>%
  mutate(state_name = factor(state_name, levels = c("Cold", "Average", "Hot")))

# ------------------------------------------------------------------------------------
# -------------------------------------- STREAK VISUALIZATION ------------------------
# ------------------------------------------------------------------------------------
for (yr in c(2022, 2023, 2024, 2025)) {
  p = jrod_hmm_data %>%
    filter(season == yr) %>%
    ggplot(aes(x = game_date, y = xwoba_roll25, color = state_name)) +
    geom_line(aes(group = 1), color = "gray50", alpha = 0.4, linewidth = 0.8) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Cold" = "#3182bd", "Average" = "#bdbdbd", "Hot" = "#de2d26")) +
    theme_minimal() +
    labs(
      title = paste("Julio Rodríguez: True Performance Streaks (", yr, ")"),
      subtitle = "3-State HMM applied to a 25-PA Rolling xwOBA",
      x = "Game Date",
      y = "25-PA Rolling xwOBA",
      color = "HMM State"
    ) +
    theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")
  
  print(p)
}
