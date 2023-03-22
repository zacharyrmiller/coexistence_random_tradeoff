# Email: zrmiller@illinois.edu
# March 2022

### Supplementary figures

library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(deSolve)
library(truncnorm)

source("./cc_functions.R")
source("./distribution_definitions.R")



##### Fig. S2 | Distribution of X_i is approximately exponential #####

m <- 1
pool_sizes <- c(2^5, 2^7, 2^9, 2^11)

# initialize tibble for simulation results
results_x <- tibble(x = numeric(), distribution = character(), id = numeric(), n = numeric())

for(n in pool_sizes){
  
  # sample one random pool using each example distribution
  pool_uniform <- sample_uniform(n, m)
  pool_exponential <- sample_exponential(n, m)
  pool_pareto <- sample_pareto(n, m)
  pool_triangle <- sample_triangle(n, m)
  
  # get the colonization rates and thresholds for each persisting species
  cl_uniform <- get_attractor(pool_uniform, m, return_bound = TRUE)
  cl_exponential <- get_attractor(pool_exponential, m, return_bound = TRUE)
  cl_pareto <- get_attractor(pool_pareto, m, return_bound = TRUE)
  cl_triangle <- get_attractor(pool_triangle, m, return_bound = TRUE)
  
  # compute the differences X_i = c_i - l_{i-1} and normalize (multiply by n f(l_{i-1}))
  x_uniform <- (cl_uniform$c - head(cl_uniform$l, -1)) * n * dunif(head(cl_uniform$l, -1) - m)
  x_exponential <- (cl_exponential$c - head(cl_exponential$l, -1)) * n * dexp(head(cl_exponential$l, -1) - m)
  x_pareto <- (cl_pareto$c - head(cl_pareto$l, -1)) * n * ((2 * m^2)/ head(cl_pareto$l, -1)^3)
  x_triangle <- (cl_triangle$c - head(cl_triangle$l, -1)) * n * (2 * (head(cl_triangle$l, -1) - m))
  
  # add distribution names
  uniform_df <- data.frame(x = x_uniform, distribution = "Uniform")
  exponential_df <- data.frame(x = x_exponential, distribution = "Exponential")
  pareto_df <- data.frame(x = x_pareto, distribution = "Pareto")
  triangle_df <- data.frame(x = x_triangle, distribution = "Triangular")
  
  # combine results
  results_x <- results_x %>% add_row(rbind(uniform_df, exponential_df, pareto_df, triangle_df) %>%
                                       group_by(distribution) %>% 
                                       mutate(id = 1:n(), n = n))
}

# plot empirical CDF of (rescaled X_i) for each distribution and n, with standard exponential for comparison
p <- results_x %>% ggplot() + 
  aes(x = x, color = distribution) + 
  stat_ecdf(pad = FALSE) +
  geom_line(aes(y = 1 - exp(-x)), linetype = "dashed", color = "black") + # add standard exponential CDF for reference
  scale_color_brewer(palette = "Dark2") + 
  facet_grid(distribution~n, scales = "free") + 
  xlab(expression(X[i])) + 
  ylab("Cumulative probability") + 
  theme_bw() + 
  theme(legend.position = "none")

show(p)

##### Fig. S3 and S4 | Actual vs. approximate threshold values #####

m <- 1
pool_sizes <- c(2^5, 2^7, 2^9, 2^11)

# initialize tibble for simulation results
results_l <- tibble(l = numeric(), distribution = character(), method = character(), id = numeric(), n = numeric())

for(n in pool_sizes){
  
  # sample one random pool using each example distribution
  pool_uniform <- sample_uniform(n, m)
  pool_exponential <- sample_exponential(n, m)
  pool_pareto <- sample_pareto(n, m)
  pool_triangle <- sample_triangle(n, m)
  
  # get actual (exact) threshold values
  l_uniform <- get_attractor(pool_uniform, m, return_bound = TRUE)$l
  l_exponential <- get_attractor(pool_exponential, m, return_bound = TRUE)$l
  l_pareto <- get_attractor(pool_pareto, m, return_bound = TRUE)$l
  l_triangle <- get_attractor(pool_triangle, m, return_bound = TRUE)$l
  
  # get approximate threshold values (using linear niche shadow approximation)
  l_approximate_uniform <- get_approximate_attractor(pool_uniform, m, return_bound = TRUE)$l
  l_approximate_exponential <- get_approximate_attractor(pool_exponential, m, return_bound = TRUE)$l
  l_approximate_pareto <- get_approximate_attractor(pool_pareto, m, return_bound = TRUE)$l
  l_approximate_triangle <- get_approximate_attractor(pool_triangle, m, return_bound = TRUE)$l
  
  # add distribution and method names
  uniform_exact <- data.frame(l = l_uniform, distribution = "Uniform", method = "exact")
  exponential_exact <- data.frame(l = l_exponential, distribution = "Exponential", method = "exact")
  pareto_exact <- data.frame(l = l_pareto, distribution = "Pareto", method = "exact")
  triangle_exact <- data.frame(l = l_triangle, distribution = "Triangular", method = "exact")
  
  # add distribution and method names
  uniform_approximate <- data.frame(l = l_approximate_uniform, distribution = "Uniform", method = "approximate")
  exponential_approximate <- data.frame(l = l_approximate_exponential, distribution = "Exponential", method = "approximate")
  pareto_approximate <- data.frame(l = l_approximate_pareto, distribution = "Pareto", method = "approximate")
  triangle_approximate <- data.frame(l = l_approximate_triangle, distribution = "Triangular", method = "approximate")
  
  # combine results
  results_l <- results_l %>% add_row(rbind(uniform_exact, exponential_exact, pareto_exact, triangle_exact,
                                           uniform_approximate, exponential_approximate, pareto_approximate, triangle_approximate) %>% 
                                       group_by(distribution, method) %>% 
                                       mutate(id = 1:n(), n = n))
}

# plot actual vs. approximate for bounded distributions on standard scale
p_bounded <- results_l %>% group_by(distribution, n) %>%
  filter(distribution %in% c("Triangular", "Uniform")) %>% # bounded distributions
  pivot_wider(names_from = method, values_from = l) %>%
  mutate(missing = is.na(exact * approximate)) %>% # indicate which values don't have a corresponding prediction/observation
  replace(is.na(.), m) %>% # plot missing values along the corresponding axis
  ggplot() + 
  aes(x = exact, y = approximate, color = distribution, shape = missing) + 
  geom_abline(slope = 1) + # add 1:1 line for reference
  geom_point() + 
  facet_grid(distribution~n, scales = "free") + 
  scale_color_manual(values = brewer.pal(4, "Dark2")[3:4]) + 
  scale_shape_manual(values = c(19, 4)) + 
  xlab("Exact") +
  ylab("Approximate") + 
  theme_bw() + 
  theme(legend.position = "none")

show(p_bounded)

# plot actual vs. approximate for unbounded distributions on log-log scale
p_unbounded <- results_l %>% group_by(distribution, n) %>%
  filter(distribution %in% c("Exponential", "Pareto")) %>% # unbounded distributions
  filter(l > 1) %>%
  pivot_wider(names_from = method, values_from = l) %>%
  mutate(missing = is.na(exact * approximate)) %>% # indicate which values don't have a corresponding prediction/observation
  replace(is.na(.), m) %>% # plot missing values along the corresponding axis
  ggplot() + 
  aes(x = exact, y = approximate, color = distribution, shape = missing) + 
  geom_abline(slope = 1) + # add 1:1 line for reference
  geom_point() + 
  facet_grid(distribution~n, scales = "free") + 
  scale_x_log10() + scale_y_log10() + # log-log scales 
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(19, 4)) + 
  xlab("Exact") +
  ylab("Approximate") +
  theme_bw() + 
  theme(legend.position = "none")

show(p_unbounded)


##### Fig. S5 | Probability of persistence by competitive rank #####

m <- 1
n_reps <- 10^5 # number of simulations for each n
pool_sizes <- c(2^5, 2^7, 2^9, 2^11)

# initialize tibble for simulation results
results <- tibble(pool_size = numeric(),
                  distribution = character(),
                  rank = numeric(),
                  fraction = numeric())

for(n in pool_sizes){
  
  print(n)
  
  # initialize one counter for each rank (to track how many times species at each rank persist)
  counter_uniform <- rep(0, n)
  counter_exponential <- rep(0, n)
  counter_pareto <- rep(0, n)
  counter_triangle <- rep(0, n)
  
  for(i in 1:n_reps){
    
    # sample one random pool using each example distribution
    pool_uniform <- sort(sample_uniform(n, m))
    pool_exponential <- sort(sample_exponential(n, m))
    pool_pareto <- sort(sample_pareto(n, m))
    pool_triangle <- sort(sample_triangle(n, m))
    
    # get set of persisting species and increment counters for those species
    counter_uniform <- counter_uniform + (pool_uniform %in% get_attractor(pool_uniform, m, ordered = TRUE))
    counter_exponential <- counter_exponential + (pool_exponential %in% get_attractor(pool_exponential, m, ordered = TRUE))
    counter_pareto <- counter_pareto + (pool_pareto %in% get_attractor(pool_pareto, m, ordered = TRUE))
    counter_triangle <- counter_triangle + (pool_triangle %in% get_attractor(pool_triangle, m, ordered = TRUE))
  }
  
  # combine results
  results <- results %>% add_row(pool_size = n,
                                 distribution = rep(c("Uniform", "Exponential", "Pareto", "Triangular"), each = n),
                                 rank = rep(1:n, times = 4),
                                 fraction = c(counter_uniform, counter_exponential, counter_pareto, counter_triangle)/n_reps)
}

# plot probability of persistence by rank
p <- results %>%
  filter(rank > 1) %>% # exclude the first species (which always persists)
  ggplot() + 
  aes(x = rank, y = fraction, color = distribution) + 
  geom_hline(yintercept = 1/2, linetype = "dashed") + # add 1/2 theory line for reference
  geom_point() + 
  stat_smooth(method = "lm", formula = y~1, se = FALSE, fullrange = TRUE) + # show average probability of persistence across all ranks
  facet_grid(distribution~pool_size, scales = "free") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = "Dark2") +
  xlab("Rank") + 
  ylab("Probability of persistence") +
  theme_bw() + 
  theme(legend.position = "none")

show(p)


##### Fig. S6 | Probability of persistence by competitive rank for truncated normal distributions #####

# Define three normal distributions truncated at m (with different modes)
sample_normal_high <- function(n, m){
  return(rtruncnorm(n = n, a = m, b = Inf, mean = m + 3, sd = 1))
}

sample_normal_low <- function(n, m){
  return(rtruncnorm(n = n, a = m, b = Inf, mean = m + 2, sd = 1))
}

sample_normal_half <- function(n, m){
  return(rtruncnorm(n = n, a = m, b = Inf, mean = m, sd = 1))
}

m <- 1
n_reps <- 10^5 # number of simulations for each n
pool_sizes <- c(2^5, 2^7, 2^9, 2^11)

# initialize tibble for simulation rsults
results <- tibble(pool_size = numeric(),
                  distribution = character(),
                  rank = numeric(),
                  fraction = numeric())

for(n in pool_sizes){
  print(n)
  
  # initialize one counter for each rank (to track how many times species at each rank persist)
  counter_high <- rep(0, n)
  counter_low <- rep(0, n)
  counter_half <- rep(0, n)
  
  for(i in 1:n_reps){
    
    # sample one random pool using each example distribution
    pool_high <- sort(sample_normal_high(n, m))
    pool_low <- sort(sample_normal_low(n, m))
    pool_half <- sort(sample_normal_half(n, m))
    
    # get set of persisting species and increment counters for those species
    counter_high <- counter_high + (pool_high %in% get_attractor(pool_high, m, ordered = TRUE))
    counter_low <- counter_low + (pool_low %in% get_attractor(pool_low, m, ordered = TRUE))
    counter_half <- counter_half + (pool_half %in% get_attractor(pool_half, m, ordered = TRUE))
  }
  
  #combine results
  results <- results %>% add_row(pool_size = n,
                                 distribution = rep(c("Mode = m + 3", "Mode = m + 2", "Mode = m"), each = n),
                                 rank = rep(1:n, times = 3),
                                 fraction = c(counter_high, counter_low, counter_half)/n_reps)
}

p_b <- results %>%
  filter(rank > 1) %>% # exclude the first species (which always persists)
  ggplot() + 
  aes(x = rank, y = fraction, color = distribution) + 
  geom_hline(yintercept = 1/2, linetype = "dashed") + # add theory 1/2 line for reference
  geom_point() + 
  stat_smooth(method = "lm", formula = y~1, se = FALSE) + # show average probability of persistence across all ranks
  facet_grid(distribution~pool_size, scales = "free") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = brewer.pal(7, "Dark2")[5:7]) +
  xlab("Rank") + 
  ylab("Probability of survival") +
  theme_bw() + 
  theme(legend.position = "none")

# generate density plots for truncated normal distributions
x_seq <- seq(m, m + 6, by = 0.01)
df_norm <- data.frame(x = rep(x_seq, 3), 
                      distribution = rep(c("Mode = m + 3", "Mode = m + 2", "Mode = m"), 
                                         each = length(x_seq)),
                      density = c(dtruncnorm(x = x_seq, a = m, b = Inf, mean = m + 3, sd = 1),
                                  dtruncnorm(x = x_seq, a = m, b = Inf, mean = m + 2, sd = 1),
                                  dtruncnorm(x = x_seq, a = m, b = Inf, mean = m, sd = 1)),
                      dummy = " ") # add a dummy variable to facet by (this helps align plots a and b)

p_a <- df_norm %>% 
  ggplot() + 
  aes(x = x, y = density, color = distribution) + 
  geom_line() + 
  geom_vline(xintercept = m, linetype = "dashed") +
  facet_grid(distribution~dummy) + # facet by dummy variable
  scale_color_manual(values = brewer.pal(7, "Dark2")[5:7]) +
  xlab("Colonization rate (c)") + 
  ylab("Probability density") + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")

#combine plots
p_combine <- ggarrange(plotlist = list(p_a, p_b), ncol = 2, labels = "auto", widths = c(1, 3))

show(p_combine)


##### Fig. S7 | Distribution of gaps (spacings) for uniformly distributed colonization rates #####

n <- 64
n_reps <- 10^4

gaps <- c() # initialize vector of gaps
for(i in 1:n_reps){
  
  pool_uniform <- sample_uniform(n, m) # sample one pool
  survivors_uniform <- get_attractor(pool_uniform, m) # get set of persisting species
  gaps_uniform <- diff(survivors_uniform) # get gaps between colonization rates in the assembled metacommunity
  
  gaps <- append(gaps, gaps_uniform) # add to vector of gaps
}

gaps_df <- data.frame(gaps = gaps)
p <- gaps_df %>% ggplot() + 
  aes(x = gaps) + 
  geom_histogram(aes(y = ..density..), bins = 50, # show distribution of simulation results
                 fill = brewer.pal(7, "Dark2")[4],
                 alpha = 0.4) + 
  #geom_density(color = brewer.pal(7, "Dark2")[4], size = 1.2) + 
  geom_line(data = data.frame(gaps = seq(0, max(gaps_df), by = 0.001)) %>% # show theory distribution
              mutate(y = n^2 * gaps * exp(-n * gaps)), # Erlang density
            aes(y = y), color = "blue", size = 1.1) + 
  geom_line(data = data.frame(gaps = seq(0, max(gaps_df), by = 0.001)) %>% # show random sampling "null model"
              mutate(y = dbeta(gaps, 1, n/2)), 
            aes(y = y),  color = "black", linetype = "dashed", size = 1.1) + 
  xlim(c(0, quantile(gaps, 0.99))) + # truncate outliers for better visualization
  ylab("Density") + 
  xlab("Spacing between rates") + 
  theme_bw()

show(p)


##### Fig. S8 | De novo assembly over long times #####

m <- 1
n_reps <- 10^3 # how many replicate trajectories for each distribution
n_invasions <- 10^4 # how many invasion attempts for each trajectory

dist_list <- list(sample_uniform, sample_exponential, sample_triangle, sample_pareto) # list of example distributions
results <- matrix(nrow = 4 * n_reps * n_invasions, ncol = 5) # initialize matrix for results

for(j in 1:4){ # loop over 4 example distributions
  
  sample_function <- dist_list[[j]] # which distribution to sample from
  for(i in 1:n_reps){
    
    if(i %% 100 == 0) print(i)
    
    c_vec <- c() # initalize empty vector for colonization rates
    for(k in 1:n_invasions){
      
      invader <- sample_function(1, m) # sample a random invader
      c_vec <- append(c_vec, invader) # add the invader to the metacommunity
      c_vec <- get_attractor(c_vec, m) # get updated attractor following invasion
      success <- invader %in% c_vec # did the invasion succeed?
      
      # add to results
      results[(j - 1) * (n_reps * n_invasions) + (i - 1) * n_invasions + k, ] <- c(i, j, k, 
                                                                                   success, 
                                                                                   length(c_vec))
    }
  }
}

# tidy up
colnames(results) <- c("rep", "distribution", "invasion_attempt", "success", "n_sp")
results <- results %>% as_tibble() %>% mutate(distribution = recode(distribution,
                                                                    `1` = "Uniform", 
                                                                    `2` = "Exponential", 
                                                                    `3` = "Triangular", 
                                                                    `4` = "Pareto"))

# summarize simulation results
results_summary <- results %>% group_by(distribution, invasion_attempt) %>% 
  summarize(mean = mean(n_sp), # get the mean richness for each distribution after tau invasion attempts
            ul = mean + sd(n_sp), ll = mean - sd(n_sp)) %>% # get the standard deviation of richness at each timepoint
  mutate(worst_case = log(invasion_attempt / 2 + 1) - digamma(1)) # add theory lower bound

# plot trajectory summaries
p_without_replacement <- results_summary %>% 
  ggplot() + 
  aes(x = invasion_attempt, color = distribution, group = distribution, fill = distribution) + 
  geom_ribbon(aes(ymin = ll, ymax = ul), alpha = 0.2, colour = NA) + # add background ribbon to show /pm 1 standard deviation
  geom_line(aes(y = mean), size = 1.2) + # plot means
  geom_line(aes(y = worst_case), linetype = "dashed", color = "black") + # add theory lower bound
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_log10() + # plot with log scale
  xlab("Invasion attempts") + 
  ylab("Number of species") + 
  theme_classic() + 
  theme(legend.title = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        panel.background = element_rect(color = "black", size = 1),
        legend.position = c(0.2, 0.82))

show(p_without_replacement)


##### Figs. S9 - S12 | Equilibrium occupancies vs. competitive rank #####

m <- 1
n_reps <- 10^4 # number of simulations for each n
pool_sizes <- c(2^5, 2^7, 2^9)

# initialize tibble for simulation results
results <- tibble(pool_size = numeric(),
                  distribution = character(),
                  rank = numeric(),
                  mean_abun = numeric())

for(n in pool_sizes){
  
  print(n)
  
  # initialize one counter for each rank (to track how many times species at each rank persist)
  counter_uniform <- rep(0, n)
  counter_exponential <- rep(0, n)
  counter_pareto <- rep(0, n)
  counter_triangle <- rep(0, n)
  
  # initialize cumulative occupancy for each rank (to calculate mean occupancy across realizations)
  p_uniform <- rep(0, n)
  p_exponential <- rep(0, n)
  p_pareto <- rep(0, n)
  p_triangle <- rep(0, n)
  
  for(i in 1:n_reps){
    
    # sample one random pool using each example distribution
    pool_uniform <- sort(sample_uniform(n, m))
    pool_exponential <- sort(sample_exponential(n, m))
    pool_pareto <- sort(sample_pareto(n, m))
    pool_triangle <- sort(sample_triangle(n, m))
    
    # get set of persisting species
    survivors_uniform <- get_attractor(pool_uniform, m, ordered = TRUE)
    survivors_exponential <- get_attractor(pool_exponential, m, ordered = TRUE)
    survivors_pareto <- get_attractor(pool_pareto, m, ordered = TRUE)
    survivors_triangle <- get_attractor(pool_triangle, m, ordered = TRUE)
    
    # increment counters for each persisting species
    counter_uniform <- counter_uniform + (pool_uniform %in% survivors_uniform)
    counter_exponential <- counter_exponential + (pool_exponential %in% survivors_exponential)
    counter_pareto <- counter_pareto + (pool_pareto %in% survivors_pareto)
    counter_triangle <- counter_triangle + (pool_triangle %in% survivors_triangle)
    
    # for each persisting species, add occupancy to cumulative sum
    p_uniform[pool_uniform %in% survivors_uniform] <- p_uniform[pool_uniform %in% survivors_uniform] + find_eq_occupancies(survivors_uniform, m, ordered = TRUE)
    p_exponential[pool_exponential %in% survivors_exponential] <- p_exponential[pool_exponential %in% survivors_exponential] + find_eq_occupancies(survivors_exponential, m, ordered = TRUE)
    p_pareto[pool_pareto %in% survivors_pareto] <- p_pareto[pool_pareto %in% survivors_pareto] + find_eq_occupancies(survivors_pareto, m, ordered = TRUE)
    p_triangle[pool_triangle %in% survivors_triangle] <- p_triangle[pool_triangle %in% survivors_triangle] + find_eq_occupancies(survivors_triangle, m, ordered = TRUE)
  }
  
  # for each pool size, calculate the mean occupancy (conditional on persistence) for each rank
  results <- results %>% 
    add_row(pool_size = n,
            distribution = rep(c("Uniform", "Exponential", "Pareto", "Triangular"), each = n),
            rank = rep(1:n, times = 4),
            mean_abun = c(p_uniform, p_exponential, p_pareto, p_triangle) / c(counter_uniform, counter_exponential, counter_pareto, counter_triangle))
}


# plot results for uniform distribution (Fig. S9)
p_uniform <- results %>%
  filter(rank > 1, distribution == "Uniform") %>%
  ggplot() + aes(x = rank, y = mean_abun, color = distribution) + 
  geom_point() + 
  facet_wrap(pool_size~., scales = "free") +
  geom_line(aes(y = 1 / (pool_size * ((rank/pool_size) + m)^(3/2))), # add theory prediction (from Eq. S42)
            linetype = "dashed", color = "black") + 
  scale_color_manual(values = brewer.pal(4, "Dark2")[4]) + 
  xlab("Rank") + 
  ylab("Mean occupancy") +
  theme_bw() + 
  theme(legend.position = "none")

show(p_uniform)

# plot results for triangular distribution (Fig. S10)
p_triangular <- results %>%
  filter(rank > 1, distribution == "Triangular") %>%
  ggplot() + aes(x = rank, y = mean_abun, color = distribution) + 
  geom_point() + 
  facet_wrap(pool_size~., scales = "free") +
  geom_line(aes(y = 1 / (pool_size * 2 * sqrt(rank/pool_size) * (sqrt(rank/pool_size) + m)^(3/2))), # add theory prediction (from Eq. S42)
            linetype = "dashed", color = "black") + 
  scale_color_manual(values = brewer.pal(4, "Dark2")[3]) + 
  xlab("Rank") + 
  ylab("Mean occupancy") +
  theme_bw() + 
  theme(legend.position = "none")

show(p_triangular)

# plot results for Pareto distribution (Fig. S11)
p_pareto <- results %>%
  filter(rank > 1, distribution == "Pareto") %>%
  ggplot() + aes(x = rank, y = mean_abun, color = distribution) + 
  geom_point() + 
  facet_wrap(pool_size~., scales = "free") +
  geom_line(aes(y = 1 / (pool_size * 2 * sqrt(m) * (1 - (rank/pool_size))^(3/4))), # add theory prediction (from Eq. S42)
            linetype = "dashed", color = "black") + 
  scale_color_manual(values = brewer.pal(4, "Dark2")[2]) + 
  xlab("Rank") + 
  ylab("Mean occupancy") +
  theme_bw() + 
  theme(legend.position = "none")

show(p_pareto)

# plot results for exponential distribution (Fig. S12)
p_exponential <- results %>%
  filter(rank > 1, distribution == "Exponential") %>%
  ggplot() + aes(x = rank, y = mean_abun, color = distribution) + 
  geom_point() + 
  facet_wrap(pool_size~., scales = "free") +
  geom_line(aes(y = 1 / (pool_size * (1 - (rank/pool_size)) * (-log(1 - (rank/pool_size)) + m)^(3/2))), # add theory prediction (from Eq. S42)
            linetype = "dashed", color = "black") + 
  scale_color_manual(values = brewer.pal(4, "Dark2")[1]) + 
  xlab("Rank") + 
  ylab("Mean occupancy") +
  theme_bw() + 
  theme(legend.position = "none")

show(p_exponential)