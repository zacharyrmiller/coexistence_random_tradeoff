library(tidyverse)
library(RColorBrewer)
library(ggpubr)

source("cc_functions.R")

##### Distributions to use throughout #####

# Uniform with UL = m + 1
sample_uniform <- function(n, m){
  return(runif(n, min = m, max = m + 1))
}

# Exponential with rate = 1
sample_exponential <- function(n, m){
  return(rexp(n, rate = 1) + m)
}

# Pareto with alpha = 2
sample_pareto <- function(n, m){
  raw <- runif(n)
  transformed <- m / (raw^(1/2))
  return(transformed)
}

# Triangular distribution with mode and upper limit at m + 1
sample_triangle <- function(n, m){
  raw <- runif(n)
  transformed <- sqrt(runif(raw)) + m
  return(transformed)
}

##### Fig. 1 | Scaling of expected richness with pool size #####

m <- 1
n_reps <- 10^5 # number of simulations for each n
pool_sizes <- unique(c(seq(2, 20, by = 2), 2^(1:10)))

results <- tibble(n = numeric(),
                  rep = numeric(),
                  Uniform = numeric(),
                  Exponential = numeric(),
                  Pareto = numeric(),
                  Triangular = numeric())

for(n in pool_sizes){
  print(n)
  for(i in 1:n_reps){
    
    pool_uniform <- sample_uniform(n, m)
    pool_exponential <- sample_exponential(n, m)
    pool_pareto <- sample_pareto(n, m)
    pool_triangle <- sample_triangle(n, m)
    
    num_survive_uniform <- length(get_attractor(pool_uniform, m))
    num_survive_exponential <- length(get_attractor(pool_exponential, m))
    num_survive_pareto <- length(get_attractor(pool_pareto, m))
    num_survive_triangle <- length(get_attractor(pool_triangle, m))
    
    results <- results %>% add_row(n = n,
                                   rep = i,
                                   Uniform = num_survive_uniform,
                                   Exponential = num_survive_exponential,
                                   Pareto = num_survive_pareto,
                                   Triangular = num_survive_triangle)
  }
}

results_tidy <- results %>% pivot_longer(cols = -c(n, rep))


# Log-log main plot
p_expectation <- results_tidy %>% 
  filter(n %in% 2^(1:10)) %>%
  group_by(n, name) %>% 
  summarize(mean = mean(value)) %>%
  ggplot() + 
  aes(x = n, y = mean, group = name, colour = name, shape = name) + 
  #annotate("text", x = 8, y = 7, label = "*", size = 10) + 
  #annotate("text", x = 64, y = 48, label = "*", size = 10) + 
  #annotate("text", x = 512, y = 384, label = "*", size = 10) + 
  geom_point(size = 3, stroke = 1.2) +
  scale_x_log10() + 
  scale_y_log10() + 
  scale_color_brewer(palette = "Dark2") + 
  geom_abline(slope = 1, intercept = -log10(2), linetype = "dashed") + 
  theme_classic() +
  theme(legend.title = element_blank(),
        panel.background = element_rect(color = "black", size = 1),
        legend.position = c(0.25, 0.75)) +
  xlab("Number of species in pool (n)") + 
  ylab("Mean number of coexisting species (k)")

# Standard scale (low n) inset
p_inset <- results_tidy %>% 
  filter(n <= 20) %>%
  group_by(n, name) %>% 
  summarize(mean = mean(value)) %>%
  ggplot() + 
  aes(x = n, y = mean, group = name, colour = name, shape = name) + 
  geom_point(size = 1.5, stroke = 1.5) +
  geom_abline(slope = 1/2, intercept = 0, linetype = "dashed") + 
  expand_limits(x = 0, y = 0) + 
  scale_x_continuous(breaks = c(0, 10, 20)) + 
  scale_y_continuous(breaks = c(0, 5, 10, 15)) +
  scale_color_brewer(palette = "Dark2") +  
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        panel.background = element_rect(color = "black", size = 1))

# combine plots
comb_plot <- p_expectation + annotation_custom(
  ggplotGrob(p_inset), 
  xmin = 1.7, xmax = 3, ymin = 0, ymax = 1.3
)

show(comb_plot)

##### Fig. 2 | Histograms showing binomial convergence #####
# using simulation results from previous section

p_histograms <- results_tidy %>% 
  filter(n %in% c(2^3, 2^6, 2^9)) %>% 
  group_by(n, value, name) %>% 
  count(value) %>% group_by(n, name) %>% 
  mutate(`Relative density` = nn / max(nn),
         theory = dbinom(value, n, 1/2) / dbinom(round(n/2), n, 1/2),
         n = paste("n =", n)) %>% 
  ggplot() + 
  aes(x = value, y = `Relative density`, colour = name, fill = name) + 
  geom_bar(stat = "identity", width = 1) + 
  geom_line(aes(y = theory), size = 1, color = "black") + 
  facet_grid(fct_rev(name) ~ fct_rev(n), scales = "free") + 
  scale_y_continuous(breaks = c(0, 0.5, 1)) + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_color_brewer(palette = "Dark2") + 
  xlab("Number of surviving species (k)") + 
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", size = 0.75))

show(p_histograms)


##### Fig. 3 | Distribution of survivor colonization rates, and distribution of gaps #####

n <- 64
n_reps <- 1000

results_uniform <- results_exponential <- results_pareto <- results_triangle <- tibble(value = numeric(),
                                                                                       type = character())
for(i in 1:n_reps){
  
  pool_uniform <- sample_uniform(n, m)
  pool_exponential <- sample_exponential(n, m)
  pool_pareto <- sample_pareto(n, m)
  pool_triangle <- sample_triangle(n, m)
  
  survivors_uniform <- get_attractor(pool_uniform, m)
  survivors_exponential <- get_attractor(pool_exponential, m)
  survivors_pareto <- get_attractor(pool_pareto, m)
  survivors_triangle <- get_attractor(pool_triangle, m)
  
  gaps_uniform <- diff(survivors_uniform)
  gaps_exponential <- diff(survivors_exponential)
  gaps_pareto <- diff(survivors_pareto)
  gaps_triangle <- diff(survivors_triangle)
  
  control_gaps_uniform <- diff(sort(sample(pool_uniform, length(survivors_uniform))))
  control_gaps_exponential <- diff(sort(sample(pool_exponential, length(survivors_exponential))))
  control_gaps_pareto <- diff(sort(sample(pool_pareto, length(survivors_pareto))))
  control_gaps_triangle <- diff(sort(sample(pool_triangle, length(survivors_triangle))))
  
  skip_gaps_uniform <- diff(sort(pool_uniform)[seq(1, n, by = 2)])
  skip_gaps_exponential <- diff(sort(pool_exponential)[seq(1, n, by = 2)])
  skip_gaps_pareto <- diff(sort(pool_pareto)[seq(1, n, by = 2)])
  skip_gaps_triangle <- diff(sort(pool_triangle)[seq(1, n, by = 2)])
                                                 
  results_uniform <- results_uniform %>% add_row(value = c(pool_uniform, survivors_uniform, gaps_uniform, 
                                                           control_gaps_uniform, skip_gaps_uniform),
                                                 type = c(rep("pool", length(pool_uniform)),
                                                          rep("survivors", length(survivors_uniform)),
                                                          rep("actual", length(gaps_uniform)),
                                                          rep("control", length(control_gaps_uniform)),
                                                          rep("skip", length(skip_gaps_uniform))))
  results_exponential <- results_exponential %>% add_row(value = c(pool_exponential, survivors_exponential, gaps_exponential, 
                                                           control_gaps_exponential, skip_gaps_exponential),
                                                 type = c(rep("pool", length(pool_exponential)),
                                                          rep("survivors", length(survivors_exponential)),
                                                          rep("actual", length(gaps_exponential)),
                                                          rep("control", length(control_gaps_exponential)),
                                                          rep("skip", length(skip_gaps_exponential))))
  results_pareto <- results_pareto %>% add_row(value = c(pool_pareto, survivors_pareto, gaps_pareto, 
                                                           control_gaps_pareto, skip_gaps_pareto),
                                                 type = c(rep("pool", length(pool_pareto)),
                                                          rep("survivors", length(survivors_pareto)),
                                                          rep("actual", length(gaps_pareto)),
                                                          rep("control", length(control_gaps_pareto)),
                                                          rep("skip", length(skip_gaps_pareto))))
  results_triangle <- results_triangle %>% add_row(value = c(pool_triangle, survivors_triangle, gaps_triangle, 
                                                           control_gaps_triangle, skip_gaps_triangle),
                                                 type = c(rep("pool", length(pool_triangle)),
                                                          rep("survivors", length(survivors_triangle)),
                                                          rep("actual", length(gaps_triangle)),
                                                          rep("control", length(control_gaps_triangle)),
                                                          rep("skip", length(skip_gaps_triangle))))
}

results_uniform <- results_uniform %>% mutate(name = "Uniform")
results_exponential <- results_exponential %>% mutate(name = "Exponential")
results_pareto <- results_pareto %>% mutate(name = "Pareto")
results_triangle <- results_triangle %>% mutate(name = "Triangular")

results_all <- rbind(results_uniform, results_exponential, results_pareto, results_triangle) %>%
  mutate(kind = ifelse(type == "pool" | type == "survivors", "c", "gaps"))


# CDF of colonization rates in pool and in survivors
p_marginal <- results_all %>% filter(kind == "c") %>%
  ggplot() + 
  aes(x = value, colour = name, 
      linetype = reorder(type, desc(type)), 
      alpha = reorder(type, desc(type))) +
  stat_ecdf(pad = "false", size = 1.2) + 
  coord_cartesian(xlim = c(1, 3)) +
  xlab("Colonization rate (c)") + 
  ylab("Cumulative probability") + 
  theme_classic() + 
  scale_alpha_discrete(range = c(0.5, 1), guide = "none") + 
  scale_linetype_discrete(guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.title = element_blank(),
        panel.background = element_rect(color = "black", size = 1),
        legend.position = c(0.75, 0.3),
  )
  
# Density plots of gaps in random sub-samples from pool vs. survivors vs. "skip-one" sub-sample from pool
p_gaps <- results_all %>% filter(kind == "gaps", type != "skip") %>%
  ggplot() + 
  aes(x = log10(value)) +
  geom_density(aes(fill = name, colour = name, alpha = type, linetype = type), size = 1) +
  geom_density(data = results_all %>% filter(kind == "gaps", type == "skip"),
               aes(x = log10(value)),
               linetype = "dotted", size = 1) + 
  #geom_line(aes(y = theory)) + 
  #geom_line(data = tibble(x = log10(seq(0, 1, by = 0.0001)),
  #                 y = dbeta(seq(0, 1, by = 0.0001), 2, n-2) * seq(0, 1, by = 0.0001) * log(10),
  #                 name = "Uniform"),
  #          aes(x = x, y = y), 
  #          linetype = "dashed", size = 1, color = "blue") + 
  # geom_line(data = tibble(x = log(seq(0, 1, by = 0.0001)),
  #                         y = dbeta(seq(0, 1, by = 0.0001), 1, n/2 - 1) * seq(0, 1, by = 0.0001)),
  #           aes(x = x, y = y), 
  #           linetype = "dotted", size = 1) +
  coord_cartesian(xlim = c(-3, -0.5)) + 
  scale_x_continuous(breaks = c(-3, -2, -1),
                     labels = c(expression(10^-3), expression(10^-2), expression(10^-1))) + 
  scale_alpha_discrete(range = c(0.1, 0.2), guide = "none") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_color_brewer(palette = "Dark2") + 
  facet_wrap(.~name, scales = "free") + 
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", size = 1)) + 
  xlab("Spacing between rates") + 
  ylab("Density")

comb_plot <- ggarrange(plotlist = list(p_marginal, p_gaps), ncol = 2, labels = "auto", widths = c(1, 1.25))

show(comb_plot)
