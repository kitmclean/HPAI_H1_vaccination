
# title: "Two-Virus SLIR Model with Vaccination: Simulation Workflow"
# author: "Kit D. McLean, Phd -- Postdoc Fellow, U of Michigan Dept of Epidemiology"

# This document simulates a two-virus SLIR model with vaccination.
# We examine how different vaccination rates (nu) affect epidemic dynamics of seasonal flu (H1) 
# and bovine-derived avian flu (H5).

###############################################################################
# 1. Setup: metadata and required packages
###############################################################################

# Load necessary packages
library(deSolve)
library(furrr)
library(future)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lhs)
library(ppcor)
library(scales)
library(GGally)
library(patchwork)
library(cowplot)
library(grid)


###############################################################################
# 2. Core model definition
#    - Parameters and base parameter list
#    - ODE system (SLIR_two_virus_vax)
#    - Initial conditions
#    - Simulation time settings
###############################################################################

#### Define Parameters and Initial Conditions

PopSize   <- 10000
CowPop    <- 280000 

H5_cow_prev_peak <- 0.40   # peak fraction of herd infected
InfCowPop_peak   <- CowPop * H5_cow_prev_peak

# Time-varying H5 prevalence in cattle: simple Gaussian epidemic
InfCowPop_t <- function(t, params) {
  t_peak <- 100
  sigma  <- 40
  
  # baseline peak prevalence = 40%
  baseline_peak <- 0.40
  
  # multiply by sampled factor (e.g., 0.5–1.5)
  peak_prev <- baseline_peak * params$cow_prev_factor
  
  InfCowPop_peak <- CowPop * peak_prev
  
  InfCowPop_peak * exp(-((t - t_peak)^2) / (2 * sigma^2))
}



# Contact rate between dairy workers and sick cattle
contact_rate_cow2human <- 1

# Transmission rates (initial placeholders; overwritten by calibration below)
beta1 <- 0.9          # initial guess for H1 (seasonal flu)
beta5 <- 1e-8         # initial guess for H5 (avian flu)

# Latent periods (1 / incubation period)
sigma1 <- 1/3   # H1
sigma5 <- 1/4   # H5

# Recovery rates (1 / infectious period)
gamma1 <- 1/5   # H1
gamma5 <- 1/9   # H5
gammac <- 1/9   # Co-infection

# Cross-immunity factors (1 == no effect)
chi15 <- 0.9      # H1 infection & recovery DECREASES susceptibility to H5
chi51 <- 0.9      # H5 infection & recovery DECREASES susceptibility to H1 

# Vaccine efficacies
VE1 <- 0.5    # H1 vaccine efficacy against H1 infection
VE5 <- 0      # H1 vaccine efficacy against H5 infection

# Vaccination recovery rate modifier (vaccinated recover faster)
delta1 <- 1.3   # H1 vaccinated recover from H1 30% faster
delta5 <- 1     # H1 vaccinated do not recover from H5 any faster

# Vaccination rate (will be overwritten)
nu <- 0

# Amplitude multiplier on cow epidemic peak (1 = baseline 40% prevalence)
cow_prev_factor <- 1

# Define a template parameter list
base_params <- list(beta1 = beta1, 
                    beta5 = beta5,
                    sigma1 = sigma1,
                    sigma5 = sigma5,
                    gamma1 = gamma1,
                    gamma5 = gamma5,
                    gammac = gammac,
                    chi15 = chi15,
                    chi51 = chi51,
                    VE1 = VE1,
                    VE5 = VE5,
                    delta1 = delta1,
                    delta5 = delta5,
                    nu = nu,
                    cow_prev_factor = cow_prev_factor
)  


#### 2.2 ODE system: SLIR_two_virus_vax

SLIR_two_virus_vax <- function(t, x, params) {
  
  with(as.list(c(as.list(x), params)), {
    
    N <- sum(x)
    
    # Forces of infection
    lambda1     <- beta1 * (I1S5 + I1L5 + I1I5a + I1I5b + I1R5 +
                              vax_I1S5 + vax_I1L5 + vax_I1I5a + vax_I1I5b + vax_I1R5) / N
    
    lambda5 <- beta5 * InfCowPop_t(t, params) * contact_rate_cow2human
    
    
    lambda1_vax <- (1 - VE1) * lambda1
    lambda5_vax <- (1 - VE5) * lambda5
    
    ## Unvaccinated
    dS1S5    <- -lambda1 * S1S5 - lambda5 * S1S5 - nu * S1S5 # 
    dL1S5    <- lambda1 * S1S5 - sigma1 * L1S5 - lambda5 * L1S5 #
    dI1S5    <- sigma1 * L1S5 - gamma1 * I1S5 - lambda5 * I1S5 #
    dR1S5    <- gamma1 * I1S5 - chi15 * lambda5 * R1S5 #
    dR1L5    <- chi15 * lambda5 * R1S5 - sigma5 * R1L5 #
    dR1I5    <- sigma5 * R1L5 - gamma5 * R1I5 + gamma1 * I1L5 #
    dI1L5    <- lambda5 * I1S5 + lambda5 * L1S5 - sigma5 * I1L5 - gamma1 * I1L5 #
    dI1I5a   <- sigma5 * I1L5 - gammac * I1I5a #
    dS1L5    <- lambda5 * S1S5 - sigma5 * S1L5 - lambda1 * S1L5 #
    dS1I5    <- sigma5 * S1L5 - gamma5 * S1I5 - lambda1 * S1I5 #
    dS1R5    <- gamma5 * S1I5 - chi51 * lambda1 * S1R5 #
    dL1R5    <- chi51 * lambda1 * S1R5 - sigma1 * L1R5 #
    dI1R5    <- sigma1 * L1R5 - gamma1 * I1R5 + gamma5 * L1I5 #
    dL1I5    <- lambda1 * S1L5 + lambda1 * S1I5 - sigma1 * L1I5 - gamma5 * L1I5 #
    dI1I5b   <- sigma1 * L1I5 - gammac * I1I5b #
    
    ## Vaccinated
    dvax_S1S5  <- nu * S1S5 - lambda1_vax * vax_S1S5 - lambda5_vax * vax_S1S5 #
    dvax_L1S5  <- lambda1_vax * vax_S1S5 - sigma1 * vax_L1S5 - lambda5_vax * vax_L1S5 #
    dvax_I1S5  <- sigma1 * vax_L1S5 - delta1 * gamma1 * vax_I1S5 - lambda5_vax * vax_I1S5 #
    dvax_R1S5  <- delta1 * gamma1 * vax_I1S5 - chi15 * lambda5_vax * vax_R1S5 #
    dvax_R1L5  <- chi15 * lambda5_vax * vax_R1S5 - sigma5 * vax_R1L5 #
    dvax_R1I5  <- sigma5 * vax_R1L5 - delta5 * gamma5 * vax_R1I5 + delta1 * gamma1 * vax_I1L5 #
    dvax_I1L5  <- lambda5_vax * vax_I1S5 + lambda5_vax * vax_L1S5 - sigma5 * vax_I1L5 - delta1 * gamma1 * vax_I1L5 # 
    dvax_I1I5a <- sigma5 * vax_I1L5 - gammac * vax_I1I5a #
    dvax_S1L5  <- lambda5_vax * vax_S1S5 - sigma5 * vax_S1L5 - lambda1_vax * vax_S1L5 #
    dvax_S1I5  <- sigma5 * vax_S1L5 - delta5 * gamma5 * vax_S1I5 - lambda1_vax * vax_S1I5 #
    dvax_S1R5  <- delta5 * gamma5 * vax_S1I5 - chi51 * lambda1_vax * vax_S1R5  #
    dvax_L1R5  <- chi51 * lambda1_vax * vax_S1R5 - sigma1 * vax_L1R5 #
    dvax_I1R5  <- sigma1 * vax_L1R5 - delta1 * gamma1 * vax_I1R5 + delta5 * gamma5 * vax_L1I5 #
    dvax_L1I5  <- lambda1_vax * vax_S1L5 + lambda1_vax * vax_S1I5 - sigma1 * vax_L1I5 - delta5 * gamma5 * vax_L1I5 #
    dvax_I1I5b <- sigma1 * vax_L1I5 - gammac * vax_I1I5b #
    
    ## Fully immune
    dR1R5   <- gamma1 * I1R5 + gamma5 * R1I5 + gammac * (I1I5a + I1I5b) +
      delta1 * gamma1 * vax_I1R5 + delta5 * gamma5 * vax_R1I5 + gammac * (vax_I1I5a + vax_I1I5b) #
    
    ## Cumulative coinfections
    dC_coinfection <- sigma5 * I1L5 + sigma1 * L1I5 +
      sigma5 * vax_I1L5 + sigma1 * vax_L1I5
    
    ## Cumulative H1 infections
    dC_H1 <- lambda1 * (S1S5 + S1L5 + S1I5 + S1R5) +
      lambda1_vax * (vax_S1S5 + vax_S1L5 + vax_S1I5 + vax_S1R5)
    
    ## Cumulative H5 infections
    dC_H5 <- lambda5 * (S1S5 + L1S5 + I1S5 + R1S5) +
      lambda5_vax * (vax_S1S5 + vax_L1S5 + vax_I1S5 + vax_R1S5)
    
    ## Cumulative vaccinations 
    dC_Vax <- nu * S1S5
    
    ## Return derivatives
    list(c(
      dS1S5, dL1S5, dI1S5, dR1S5, dR1L5, dR1I5,
      dI1L5, dI1I5a, dS1L5, dS1I5, dS1R5, dL1R5,
      dI1R5, dL1I5, dI1I5b,
      dvax_S1S5, dvax_L1S5, dvax_I1S5, dvax_R1S5,
      dvax_R1L5, dvax_R1I5, dvax_I1L5, dvax_I1I5a,
      dvax_S1L5, dvax_S1I5, dvax_S1R5, dvax_L1R5,
      dvax_I1R5, dvax_L1I5, dvax_I1I5b,
      dR1R5, dC_coinfection, dC_H1, dC_H5, dC_Vax
    ))
  })
}


#### Set Initial Conditions

# Initialize all compartments to zero
x0 <- numeric(35)

# Set initial conditions
x0[1] <- PopSize - 10   # S1S5: fully susceptible, unvaccinated
x0[3] <- 9              # I1S5: infected with H1 only
x0[10] <- 1             # S1I5: infected with H5 only

# State variable names
names(x0) <- c("S1S5","L1S5","I1S5","R1S5","R1L5","R1I5",
               "I1L5","I1I5a","S1L5","S1I5","S1R5","L1R5","I1R5",
               "L1I5","I1I5b","vax_S1S5","vax_L1S5","vax_I1S5","vax_R1S5",
               "vax_R1L5","vax_R1I5","vax_I1L5","vax_I1I5a","vax_S1L5","vax_S1I5",
               "vax_S1R5","vax_L1R5","vax_I1R5","vax_L1I5","vax_I1I5b","R1R5", 
               "C_coinfection", "C_H1", "C_H5", "C_Vax")


#### Simulation Settings

# Time steps
times <- seq(0, 200, by = 1)  # Daily for 200 days

season_length <- 200  # days 


###############################################################################
# 3. Calibration utilities and calibrated parameter values
#    - Helper functions to calibrate:
#        * beta1 to target H1 attack rate
#        * beta5 to target H5 attack rate
#        * nu to target final vaccination coverage
#    - Apply calibration and store calibrated values
###############################################################################

#### Helper to compute H5 attack rate for a given beta5 ----------------------

compute_H5_attack <- function(beta5_value,
                              base_params,
                              x0,
                              times) {
  
  params <- base_params
  params$beta5  <- beta5_value
  
  # Remove H1 and all vaccination effects to match calibration setup
  params$beta1  <- 0
  params$nu     <- 0
  params$VE1    <- 0
  params$VE5    <- 0
  params$delta1 <- 1
  params$delta5 <- 1
  
  # Fully susceptible workers, no initial H1 or H5
  x0_cal <- x0
  x0_cal[]       <- 0
  x0_cal["S1S5"] <- PopSize
  
  out <- ode(y = x0_cal, times = times,
             func = SLIR_two_virus_vax, parms = params)
  df  <- as.data.frame(out)
  
  tail(df$C_H5, 1) / PopSize
}

#### Calibration utilities ---------------------------------------------------

calibrate_beta1 <- function(target_attack,
                            base_params,
                            x0,
                            times,
                            beta_interval = c(0.05, 1)) {
  
  attack_diff <- function(beta1_guess) {
    params <- base_params
    params$beta1  <- beta1_guess
    params$beta5  <- 0
    params$nu     <- 0
    params$VE1    <- 0
    params$VE5    <- 0
    params$delta1 <- 1
    params$delta5 <- 1
    
    x0_cal <- x0
    x0_cal["S1S5"] <- PopSize - 10
    x0_cal["I1S5"] <- 10
    x0_cal["S1I5"] <- 0
    
    out <- ode(y = x0_cal, times = times,
               func = SLIR_two_virus_vax, parms = params)
    df  <- as.data.frame(out)
    
    tail(df$C_H1, 1) / PopSize - target_attack
  }
  
  uniroot(attack_diff, interval = beta_interval)$root
}

calibrate_beta5_attack <- function(target_attack_H5,
                                   base_params,
                                   x0,
                                   times,
                                   beta_interval = c(1e-10, 1e-8)
) {
  
  attack_diff <- function(beta5_guess) {
    compute_H5_attack(beta5_guess,
                      base_params = base_params,
                      x0          = x0,
                      times       = times) -
      target_attack_H5
  }
  
  uniroot(attack_diff, interval = beta_interval)$root
}

calibrate_nu <- function(target_coverage,
                         base_params,
                         x0,
                         times,
                         nu_interval = c(0, 1)) {
  
  coverage_diff <- function(nu_guess) {
    params <- base_params
    params$nu <- nu_guess
    
    out <- ode(y = x0, times = times,
               func = SLIR_two_virus_vax, parms = params)
    df  <- as.data.frame(out)
    
    tail(df$C_Vax, 1) / PopSize - target_coverage
  }
  
  uniroot(coverage_diff, interval = nu_interval)$root
}

#### Calibrate parameters ----------------------------------------------------

target_attack_H1 <- 0.25

beta1_calibrated <- calibrate_beta1(target_attack = target_attack_H1,
                                    base_params   = base_params,
                                    x0            = x0,
                                    times         = times)
## Sweep over a beta5 grid and compute attack rates

# Choose a grid of beta5 values around search interval
beta5_grid <- 10^seq(-12, -7, length.out = 25)

H5_attack_grid <- sapply(beta5_grid, compute_H5_attack,
                         base_params = base_params,
                         x0          = x0,
                         times       = times)

diagnostic_df <- data.frame(
  beta5  = beta5_grid,
  attack = H5_attack_grid
)


target_attack_H5 <- 0.07

beta5_calibrated <- calibrate_beta5_attack(
  target_attack_H5 = target_attack_H5,
  base_params      = base_params,
  x0               = x0,
  times            = times,
  beta_interval    = c(1e-10, 1e-8)
)

desired_coverages <- c(0, 0.05, 0.25, 0.60, 0.80)

nu_values <- sapply(desired_coverages, function(cov) {
  calibrate_nu(target_coverage = cov,
               base_params     = base_params,
               x0              = x0,
               times           = times,
               nu_interval     = c(0, 0.1))
})

nu_min_lhs <- min(nu_values)
nu_max_lhs <- max(nu_values)

base_params$beta1 <- beta1_calibrated
base_params$beta5 <- beta5_calibrated


###############################################################################
# Main Figure 1: Baseline dynamics → coinfection emergence
# Panel A: H1 and H5 incidence over time (baseline, nu=0)
# Panel B: H1–H5 coinfection incidence over time
###############################################################################

# --- 1) Run baseline simulation (no vaccination) ---
params_baseline <- base_params
params_baseline$nu <- 0  # baseline, no vaccination

out_base <- ode(
  y     = x0,
  times = times,
  func  = SLIR_two_virus_vax,
  parms = params_baseline
)

df_base <- as.data.frame(out_base)

# --- 2) Convert cumulative counters to DAILY incidence ---
# incidence[t] = C[t] - C[t-1]
inc_df <- df_base %>%
  dplyr::select(time, C_H1, C_H5, C_coinfection) %>%
  mutate(
    inc_H1   = c(NA, diff(C_H1)),
    inc_H5   = c(NA, diff(C_H5)),
    inc_coin = c(NA, diff(C_coinfection))
  ) %>%
  filter(!is.na(inc_H1)) %>%
  mutate(
    # express as "per 1,000 workers per day" for readability
    inc_H1_per1000   = inc_H1   / PopSize * 1000,
    inc_H5_per1000   = inc_H5   / PopSize * 1000,
    inc_coin_per1000 = inc_coin / PopSize * 1000
  )

# --- 3) Build Panel A (H1 and H5 incidence) ---
panelA_df <- inc_df %>%
  dplyr::select(time, inc_H1_per1000, inc_H5_per1000) %>%
  pivot_longer(
    cols = c(inc_H1_per1000, inc_H5_per1000),
    names_to = "virus",
    values_to = "incidence"
  ) %>%
  mutate(
    virus = recode(
      virus,
      inc_H1_per1000 = "Human H1 incidence",
      inc_H5_per1000 = "Cattle H5 prevalence (external reservoir)"
    )
  )

pA <- ggplot(panelA_df, aes(x = time, y = incidence, linetype = virus)) +
  geom_line(linewidth = 1) +
  labs(
    x = NULL,
    y = "",
    title = "A. Baseline H1 and H5 incidence (no vaccination)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

# NOTE: If H5 is visually dwarfed by H1 use a square-root transform 
# pA <- pA + scale_y_continuous(trans = "sqrt")

# --- 4) Build Panel B (coinfection incidence) ---
pB <- ggplot(inc_df, aes(x = time, y = inc_coin_per1000)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Day of season",
    y = "Human incidence (H1) and cattle infection prevalence (H5)",
    title = "B. Baseline H1–H5 coinfection incidence (no vaccination)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )

# --- 5) Combine into one dense main figure ---
fig1 <- pA / pB + plot_layout(heights = c(1, 1))

fig1

# # --- 6) Save ---
# ggsave("MainFigure1_baseline_dynamics_coinfection.png", fig1,
#        width = 7.2, height = 6.2, dpi = 600)
# 
# ggsave("MainFigure1_baseline_dynamics_coinfection.pdf", fig1,
#        width = 7.2, height = 6.2)



# --- Helper: run scenario and extract summary metrics ---
run_scenario_summary <- function(nu_value, target_cov_label) {
  params <- base_params
  params$nu <- nu_value

  out <- ode(y = x0, times = times, func = SLIR_two_virus_vax, parms = params)
  df  <- as.data.frame(out)

  # daily H1 prevalence = all infectious-with-H1 compartments / PopSize
  df <- df %>%
    mutate(
      H1_prev = (I1S5 + I1L5 + I1I5a + I1I5b + I1R5 +
                   vax_I1S5 + vax_I1L5 + vax_I1I5a + vax_I1I5b + vax_I1R5) / PopSize
    )
  

  data.frame(
    target_coverage    = target_cov_label,
    nu                = nu_value,
    achieved_coverage  = tail(df$C_Vax, 1) / PopSize,
    H1_attack_rate     = tail(df$C_H1, 1) / PopSize,
    H5_attack_rate     = tail(df$C_H5, 1) / PopSize,
    coinfections_total = tail(df$C_coinfection, 1),
    coinfections_per1000 = tail(df$C_coinfection, 1) / PopSize * 1000,
    peak_H1_prev       = max(df$H1_prev, na.rm = TRUE)   # <-- Panel B
  )
}


# --- Build scenario table across desired coverages ---
scenario_table <- do.call(rbind, lapply(seq_along(desired_coverages), function(i) {
  run_scenario_summary(nu_values[i], desired_coverages[i])
}))

# --- Add absolute and relative reductions vs baseline (coverage = 0) ---
baseline <- scenario_table$coinfections_total[1]

scenario_table <- scenario_table %>%
  dplyr::mutate(
    coinfections_averted_total = baseline - coinfections_total,
    coinfections_averted_per1000 = coinfections_averted_total / PopSize * 1000,
    coinfections_reduction_pct =
      100 * (baseline - coinfections_total) / baseline
  
  ) %>%
  dplyr::arrange(target_coverage)

scenario_table


# Main Figure 2: Effect of vaccination on H1 dynamics
# Panels:
#  A: H1 attack rate vs vaccination coverage
#  B: Peak H1 prevalence vs vaccination coverage

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)


# ---- Choose which coverage measure to plot on x-axis ----
# Use achieved coverage if it is very different from target; otherwise target is fine.
coverage_var <- "achieved_coverage"   # or "target_coverage"

# ---- Panel A data: attack rate summary ----
panelA_df <- scenario_table %>%
  mutate(coverage = .data[[coverage_var]]) %>%
  dplyr::select(coverage, H1_attack_rate) %>%
  arrange(coverage)

# ---- Panel B data: compute peak prevalence and CI from time series ----

coverage_var <- "achieved_coverage"  # or "target_coverage"

plot_df <- scenario_table %>%
  mutate(coverage = .data[[coverage_var]]) %>%
  arrange(coverage)

# Panel A: attack rate vs coverage
pA <- ggplot(plot_df, aes(x = coverage, y = H1_attack_rate)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "A",
    x = "Seasonal influenza vaccination coverage",
    y = "H1 attack rate (200-day cumulative)"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Panel B: peak prevalence vs coverage
pB <- ggplot(plot_df, aes(x = coverage, y = peak_H1_prev)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "B",
    x = "Seasonal influenza vaccination coverage",
    y = "Peak H1 prevalence"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

fig2 <- pA + pB + plot_layout(ncol = 2)

ggsave("MainFigure2_H1_vax_effects.png", fig2, width = 10, height = 4, dpi = 600)
ggsave("MainFigure2_H1_vax_effects.pdf", fig2, width = 10, height = 4)

fig2
