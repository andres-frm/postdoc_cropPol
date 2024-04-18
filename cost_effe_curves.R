pks <- c('tidyverse', 'rethinking', 'rstan', 'magrittr', 'cmdstanr',
         'ggdag', 'dagitty', 'readxl', 'brms', 'cowplot', 'parallel', 
         'compiler')

sapply(pks, library, character.only = T)

options(mc.core = parallel::detectCores())

source('functions_mod_diagnostics.r')

sim_mod <- 
  dagify(Y ~ FW + nF, 
         nF ~ FloD + Pfs + percF, 
         FloD ~ percF,
         FW ~ PollD,
         Pfs ~ nVis,
         PollD ~ nVis + wild,
         nVis ~ apis, 
         apis ~ propF + apisN + timeF,
         BHd ~ BHs,
         propF ~ BHs + BHd,
         exposure = 'BHd', 
         outcome = 'Y')

coordinates(sim_mod) <- 
  list(
    y = c(BHs = 1, propF = 0, BHd = -1, apis = 0, 
          apisN = 1, timeF = -1,  nVis = 0, PollD = -1, wild = -2, Pfs = 1,
          nF = 1, FW = -1, FloD = 2, Y = 0, percF = 1.5), 
    
    x = c(BHs = 1, propF = 2, BHd = 1, apis = 3, 
          apisN = 3,  timeF = 3, nVis = 4, PollD = 5, wild = 5, Pfs = 5, 
          nF = 6, FW = 6, FloD = 6, Y = 7, percF = 7)
  )

labels <- tibble(lab = c('Pre-pollination', 'Pollination', 'Production'), 
                 x = c(2, 4.5, 6.5), 
                 y = rep(2.75, 3))

ggdag(sim_mod) +
  geom_rect(aes(xmin = 0.5, xmax = 3.5, ymin = -2.5, ymax = 2.5), 
            fill = 'tan1', alpha = 0.01) +
  geom_rect(aes(xmin = 3.5, xmax = 5.5, ymin = -2.5, ymax = 2.5), 
            fill = 'lightpink', alpha = 0.01) +
  geom_rect(aes(xmin = 5.5, xmax = 7.5, ymin = -2.5, ymax = 2.5), 
            fill = 'blue4', alpha = 0.01) +
  annotate(geom = 'text', x = labels$x, y = labels$y, label = labels$lab) +
  geom_dag_node(color = 'lightblue') +
  geom_dag_text(color = 'black') +
  theme_classic() +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())

# ========== Hive density and strength ===========

# Function to simulate hives, of population N and around 40% of foragers 
# individuals. Number of bees is estimated from a NegBin distribution
# and the proportion of foragers with a beta distribution.

farrar_law <- 
  tribble(~population, ~workers, ~proportion,
          10000, 2000, 0.2,
          20000, 5000, 0.25, 
          30000, 10000, 0.3,
          40000, 20000, 0.5,
          50000, 30000, 0.6,
          60000, 39000, 0.65)

farrar_law %$% plot(population, proportion, type = 'l')

plot(density(rbeta(1e3, 2.5, 8)))

hives_ha <- function(n_hives = 25, # n hives to be simulated
                     mu_pop = 1e4, # mu parameter NegBin distribution. N Bees per hive
                     sigma_pop = 10, # dispersion parameter NegBin distribution. N Bees per hive
                     beta_1 = 2.5, # parameter 1 beta distribution. Proportion of foragers
                     beta_2 = 8, # parameter 1 beta distribution. Proportion of foragers
                     seed = 1){
  
  hive_ha <- n_hives # hives per ha
  
  set.seed(seed)
  beehive_strength <- rnbinom(hive_ha, size = sigma_pop, mu = mu_pop)
  
  set.seed(seed)
  prop_foragers <- rbeta(hive_ha, beta_1, beta_2)
  
  set.seed(seed)
  foragers_bees <- round(beehive_strength * sample(prop_foragers, hive_ha))
  
  return(foragers_bees)
}

hives_ha(1, seed = rpois(1, 100))

# ============ Foraging pattern ========

trip_frequency_HQ <- rnorm(5e4, 
                           mean(c(1.23, 1.21)), 
                           mean(c(0.06, 0.15))) # high quality hives

trip_span_HQ <- rnorm(5e4, 
                      mean(c(1667.7, 2683.9)), 
                      mean(c(126.4, 250.7))) # low quality hives

trips_bee_HQ <- trip_frequency_HQ * ((8*60)/90) # number of trips per day

foraging_time_HQ <- trip_span_HQ * trips_bee_HQ # foraging time in seconds


trip_frequency_LQ <- rnorm(5e4, 
                           mean(c(1, 1.4)), 
                           mean(c(0, 0.16))) # low quality hives

trip_span_LQ <- rnorm(5e4, 
                      mean(c(1210.8, 1680.4)), 
                      mean(c(157.6, 212.6))) # low quality hives

trips_bee_LQ <- trip_frequency_LQ * ((8*60)/90) # number of trips per day

foraging_time_LQ <- trip_span_LQ * trips_bee_LQ # foraging time in seconds


par(mfrow = c(2, 2), mar = c(4.2, 4.2, 1, 1))

plot(density(trip_frequency_HQ), col = 'gold3', main = '', ylim = c(0, 5.5),
     xlab = 'Number of foraging trips (trips/90 min)', lwd = 2)
lines(density(trip_frequency_LQ), col = 'cyan4', main = '',
      xlab = 'Number of foraging trips (trips/90 min)', lwd = 2)

plot(density(trip_span_HQ/60), col = 'gold3', main = '', xlim = c(10, 50),
     xlab = 'Trip span (min)', lwd = 2)
lines(density(trip_span_LQ/60), main = '', col = 'cyan4',
      xlab = 'Trip span (min)', lwd = 2)

plot(density(trips_bee_HQ), main = '', col = 'gold3',
     xlab = 'Daily trip per bee', lwd = 2, ylim = c(0, 1))
lines(density(trips_bee_LQ), main = '', col = 'cyan4',
      xlab = 'Daily trip per bee (min)', lwd = 2)

plot(density((foraging_time_HQ/60)/60), col = 'gold3', main = '', ylim = c(0, 1.1),
     xlab = 'Daily foraging time per bee (h)', lwd = 2, xlim = c(1, 6.5))
lines(density((foraging_time_LQ/60)/60), main = '', col = 'cyan4',
      xlab = 'Daily foraging time per bee (h)', lwd = 2)

par(mfrow = c(1, 1))

#(feral_pop <- (quantile(1:sum(beehive_strength), 0.1)))

# ======== Time of visits ======

visit_span <- readRDS('honeybee_visit_span.rds')

visit_span <- visit_span$A_mellifera

plot(density(rlnorm(1e3, exp(0.5), exp(0.1))))
plot(density(rgamma(1e3, exp(2)/2, 1/2)))

cat(file = 'honeybee_visit.stan', 
    '
    data{
      int N;
      vector[N] visit;
      array[N] int individual;
    }
    
    parameters{
      vector[N] z_alpha;
      real<lower = 0> mu;
      real<lower = 0> sigma;
      real<lower = 0> sigma1;
    }
    
    transformed parameters{
      vector[N] alpha;
      alpha = mu + z_alpha * sigma1;
    }
    
    model{
      vector[N] p;
      mu ~ normal(5, 1.5);
      sigma ~ exponential(1);
      sigma1 ~ exponential(1);
      z_alpha ~ normal(0, 1);
    
      for (i in 1:N) {
        p[i] = alpha[individual[i]];
        p[i] = exp(p[i]);
      }
      
      visit ~ gamma(p/sigma, 1/sigma);
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p;
      array[N] real ppcheck;
      
      for (i in 1:N) {
        p[i] = alpha[individual[i]];
        p[i] = exp(p[i]);
      }
    
      for (i in 1:N) log_lik[i] = gamma_lpdf(visit[i] | p[i]/sigma, 1/sigma);
    
      ppcheck = gamma_rng(p/sigma, 1/sigma);
    }
    ')

file <- paste(getwd(), '/honeybee_visit.stan', sep = '')
fit_visit_honeybee <- cmdstan_model(file, compile = T)

mod_visit_honeybee <- 
  fit_visit_honeybee$sample(
    data = list(visit = visit_span, 
                N = length(visit_span), 
                individual = 1:length(visit_span)),
    iter_sampling = 30e3, 
    iter_warmup = 500, 
    chains = 3, 
    parallel_chains = 3, 
    thin = 3, 
    seed = 123, 
    refresh = 500
  )

(output_mod_visit_honeybee <- mod_visit_honeybee$summary())

mod_diagnostics(mod_visit_honeybee, output_mod_visit_honeybee)

par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))
for (i in 1:9) trace_plot(mod_visit_honeybee, output_mod_visit_honeybee$variable[i], 3)
par(mfrow = c(1, 1))


ppcheck_visit_honeybee <- mod_visit_honeybee$draws(variables = 'ppcheck', 
                                                   format = 'matrix')

plot(density(ppcheck_visit_honeybee[1, ]), ylim = c(0, 0.09), lwd = 0.1, 
     main = '', xlab = 'Honeybee time of floral visit (s)')
for (i in 1:100) lines(density(ppcheck_visit_honeybee[i, ]), lwd = 0.1)
lines(density(visit_span), col = 'red', lwd = 2)

time_visit_honeybee <- as.vector(ppcheck_visit_honeybee)

plot(density(time_visit_honeybee), col = 'tan1', lwd = 3)

# ======== foraging trip per day ======

# function to estimate the number of foraging trip per individual i in 
# hive j, depending on the quality of the hive.

N_foraging_trips <- 
  
  function(iter, # bees per hive 
           time_foragin, # effective time to forage during the day
           time_per_visit, # posterior distribution of time spend per visit
           seed = 123) { 
    
    set.seed(seed)
    time_foragin <- sample(time_foragin, iter, T)
    
    t1 <- Sys.time()
    visits_bee <- 
      sapply(1:iter, FUN = 
               function(x) {
                 
                 time <- 0
                 vis <- 0
                 
                 while (time < time_foragin[[x]]) {
                   time <- time + sample(time_per_visit, 1)
                   vis <- vis + 1
                 }
                 
                 vis
                 
               })
    message(paste('Execution time', Sys.time() - t1))
    
    return(visits_bee)
    
  }

t1 <- Sys.time()
visits_day <- N_foraging_trips(iter = 2e4, 
                               time_foragin = foraging_time_LQ, 
                               time_per_visit = time_visit_honeybee)
Sys.time() - t1

t1 <- Sys.time()
visits_day_HQ <- N_foraging_trips(iter = 2e4, 
                               time_foragin = foraging_time_HQ, 
                               time_per_visit = time_visit_honeybee)
Sys.time() - t1


plot(density(visits_day), main = '', xlab = 'Number of floral visits of\n a single forager honeybee bee per day',
     lwd = 4, col = 'lightblue', xlim = c(395, 2000))
lines(density(visits_day_HQ), lwd = 4, col = 'tan1')

# ====== Plants ha, flower density and percentage ======= 

# ====== 1. fruit set =====

dat_fs <- readRDS('fruit_set.rds')

dat_fs <- dat_fs[, c("plant_id", "farm_id", "year_id",
                     "beehive", "treatment", "flowers", "fruits")]

dat_fs$treatment <- ifelse(dat_fs$treatment == 1, 'close', 'open')

dat_fs |> 
  ggplot(aes(as.factor(treatment), fruits)) +
  geom_boxplot()

dat_fs2 <- lapply(6:7, 
                  function(x) {
                    t <- read_xlsx('tesis_pablo.xlsx', 
                                   shee = x, 
                                   col_names = T, 
                                   na = 'NA')[, 1:5]
                    
                    t$locality <- 'entre_rios'
                    t$treatment <- 'open'
                    
                    colnames(t) <- c('farm_id', 'plant_id', 'branch_id',
                                     'flowers', 'fruits', 'locality_id', 
                                     'treatment')

                    t$branch_id <- t %$% paste(farm_id,
                                              plant_id,
                                              branch_id,
                                              sep = '')

                    t$plant_id <- t %$% paste(farm_id,
                                              plant_id,
                                              sep = '')
                    
                    t
                  })

names(dat_fs2) <- paste('y', c(2016, 2021), sep = '')
dat_fs2$y2016$year_id <- '2016'
dat_fs2$y2021$year_id <- '2021'

dat_fs$branch_id <- 1
dat_fs$locality_id <- 'tucuman'


dat_fs <- dat_fs[, c("farm_id", "plant_id", 
                     'branch_id', "flowers", 
                     "fruits", "locality_id", 
                     "treatment", "year_id")]

dat_fs$branch_id <- paste(dat_fs$farm_id, dat_fs$plant_id, 
                          dat_fs$branch_id, sep = '_')
dat_fs$plant_id <- paste(dat_fs$farm_id, 
                         dat_fs$plant_id, sep = '_')


dat_fs <- rbind(dat_fs, dat_fs2$y2016, dat_fs2$y2021)

dat_fs <- na.omit(dat_fs)

dat_fs <- 
  lapply(dat_fs, 
         function(x) {
           if (is.character(x)) as.numeric(as.factor(x))
           else x
         })

unlist(lapply(dat_fs, function(x) sum(is.na(x))))

dat_fs$N <- length(dat_fs$plant_id)
dat_fs$N_site <- length(unique(dat_fs$locality_id))
dat_fs$N_farm <- length(unique(dat_fs$farm_id))
dat_fs$N_plant <- length(unique(dat_fs$plant_id))
dat_fs$N_branch <- length(unique(dat_fs$branch_id))
dat_fs$N_treatment <- length(unique(dat_fs$treatment))
dat_fs$N_year <- length(unique(dat_fs$year_id))


plot(density(rnorm(1e3, 50, 25)))

cat(file = 'fruit_set.stan', 
    "
    data{
      int N;
      int N_site;
      int N_farm;
      int N_plant;
      //int N_branch;
      int N_treatment;
      int N_year;
      array[N] int flowers;
      array[N] int fruits;
      array[N] int year_id;
      array[N] int locality_id;
      array[N] int farm_id;
      array[N] int plant_id;
      //array[N] int branch_id;
      array[N] int treatment;
    }
    
    parameters{
      vector[N_treatment] t;
      
      matrix[N_treatment, N_year] Z_year;
      cholesky_factor_corr[N_treatment] R_year;
      vector<lower = 0>[N_treatment] sigma_year;
      
      matrix[N_treatment, N_site] Z_locality;
      cholesky_factor_corr[N_treatment] R_locality;
      vector<lower = 0>[N_treatment] sigma_locality;
    
      matrix[N_treatment, N_farm] Z_farm;
      cholesky_factor_corr[N_treatment] R_farm;
      vector<lower = 0>[N_treatment] sigma_farm;
    
      matrix[N_treatment, N_plant] Z_plant;
      cholesky_factor_corr[N_treatment] R_plant;
      vector<lower = 0>[N_treatment] sigma_plant;
      
      //vector[N_branch] branch;
      
    }
    
    transformed parameters{
      matrix[N_year, N_treatment] year;
      matrix[N_site, N_treatment] locality;
      matrix[N_farm, N_treatment] farm;
      matrix[N_plant, N_treatment] plant;
    
      year = (diag_pre_multiply(sigma_year, R_year) * Z_year)';
      locality = (diag_pre_multiply(sigma_locality, R_locality) * Z_locality)';
      farm = (diag_pre_multiply(sigma_farm, R_farm) * Z_farm)';
      plant = (diag_pre_multiply(sigma_plant, R_plant) * Z_plant)';
    }
    
    model{
      vector[N] p;
      t ~ normal(0, 1);
    
      to_vector(Z_year) ~ normal(0, 0.5);
      R_year ~ lkj_corr_cholesky(2);
      sigma_year ~ exponential(1);
    
      to_vector(Z_locality) ~ normal(0, 0.5);
      R_locality ~ lkj_corr_cholesky(2);
      sigma_locality ~ exponential(1);
    
      to_vector(Z_farm) ~ normal(0, 0.5);
      R_farm ~ lkj_corr_cholesky(2);
      sigma_farm ~ exponential(1);
    
      to_vector(Z_plant) ~ normal(0, 0.5);
      R_plant ~ lkj_corr_cholesky(2);
      sigma_plant ~ exponential(1);
      //branch ~ normal(0, 1);
    
      for (i in 1:N) {
        p[i] = t[treatment[i]] + 
                locality[locality_id[i], treatment[i]] + 
                farm[farm_id[i], treatment[i]] +
                plant[plant_id[i], treatment[i]] +
                year[year_id[i], treatment[i]]; 
        
        p[i] = inv_logit(p[i]);
      }
      
      fruits ~ binomial(flowers, p);
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] p1;
      array[N] int ppcheck;
      matrix[N_treatment, N_treatment] Rho_year;
      matrix[N_treatment, N_treatment] Rho_locality;
      matrix[N_treatment, N_treatment] Rho_farm;
      matrix[N_treatment, N_treatment] Rho_plant;
    
      Rho_year = multiply_lower_tri_self_transpose(R_year);
      Rho_locality = multiply_lower_tri_self_transpose(R_locality);
      Rho_farm = multiply_lower_tri_self_transpose(R_farm);
      Rho_plant = multiply_lower_tri_self_transpose(R_plant);
    
      for (i in 1:N) {
        p1[i] = t[treatment[i]] + 
                locality[locality_id[i], treatment[i]] + 
                farm[farm_id[i], treatment[i]] +
                plant[plant_id[i], treatment[i]] +
                year[year_id[i], treatment[i]]; 
        
        p1[i] = inv_logit(p1[i]);
      }
      
      for (i in 1:N) log_lik[i] = binomial_lpmf(fruits[i] | flowers[i], p1[i]);
    
      ppcheck = binomial_rng(flowers, p1);
    }
    ")

file <- paste(getwd(), '/fruit_set.stan', sep = '')

fit_fruit_set <- cmdstan_model(file, compile = T)

mod_fs_exp <- 
  fit_fruit_set$sample(
    data = dat_fs, 
    chains = 3, 
    parallel_chains = 3, 
    iter_warmup = 500, 
    iter_sampling = 4e3,
    thin = 3,
    seed = 123,
    refresh = 200
  )

output_mod_fruitset <- mod_fs_exp$summary() 

mod_diagnostics(mod_fs_exp, output_mod_fruitset)

trace_plot(mod_fs_exp, output_mod_fruitset$variable[1], 3)

par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))
for (i in 2:10) {
  trace_plot(mod_fs_exp, output_mod_fruitset$variable[i], 3)
}
par(mfrow = c(1, 1))

post_fruitset <- mod_fs_exp$draws(format = 'df')

colnames(post_fruitset)

post_fruitset <- as.matrix(post_fruitset)

post_fruitset <- 
  list(treatment = post_fruitset[, grep('^t', colnames(post_fruitset))], 
       locality = post_fruitset[, grep('^locality', colnames(post_fruitset))], 
       farm = post_fruitset[, grep('^farm', colnames(post_fruitset))], 
       plant = post_fruitset[, grep('^plant', colnames(post_fruitset))], 
       #branch = post_fruitset[, grep('^branch', colnames(post_fruitset))], 
       year = post_fruitset[, grep('^year', colnames(post_fruitset))])

ppcheck_fruit_set <- mod_fs_exp$draws(variables = 'ppcheck', format = 'matrix')

plot(density(ppcheck_fruit_set[1, ]), ylim = c(0, 0.025), lwd = 0.1, 
     main = '', xlab = 'Number of fruit')
for (i in 1:100) lines(density(ppcheck_fruit_set[i, ]), lwd = 0.1)
lines(density(dat_fs$fruits), col = 'red', lwd = 2)


fruit_set <- 
  inv_logit(
    post_fruitset$treatment[, 2, drop = T] +
      apply(post_fruitset$locality, 1, mean) + 
      apply(post_fruitset$farm, 1, mean) +
      apply(post_fruitset$plant, 1, mean) +
      apply(post_fruitset$year, 1, mean) 
  )

plot(density(fruit_set), main = '', 
     xlab = 'blueberry fruit set', col = 'lightblue', lwd = 4)

mean(fruit_set)


# ======= 2. fruit produces  ====

fruit_plant <- readRDS('fruit_plant.rds')

fruit_plant$plant_id <- as.factor(paste(fruit_plant$farm, fruit_plant$plant, sep = '_'))

fruit_plant$plot <- as.factor(paste(fruit_plant$farm, fruit_plant$plot, sep = '_'))

fruit_plant <- fruit_plant[, c("farm", "plot", 'plant_id', "total_fruts")]

fruit_plant <- lapply(fruit_plant, function(x) if(is.factor(x)) as.numeric(x) else(x))

fruit_plant$N <- length(fruit_plant$farm)
fruit_plant$N_farm <- length(unique(fruit_plant$farm))
fruit_plant$N_plot <- length(unique(fruit_plant$plot))
fruit_plant$N_plant <- length(unique(fruit_plant$plant))
fruit_plant$total_fruts <- round(fruit_plant$total_fruts)

plot(density(rnbinom(1e3, size = 2, mu = exp(8))))

cat(file = 'fruits_plant.stan', 
    "
    data{
      int N;
      int N_farm;
      int N_plot;
      int N_plant;
      array[N] int plant_id;
      array[N] int total_fruts;
      array[N] int farm;
      array[N] int plot;
    }
    
    parameters{
      vector[N_plant] z_plant;
      real mu_plant;
      real<lower = 0> sigma_plant;
    
      vector[N_farm] z_farm;
      real mu_farm;
      real<lower = 0> sigma_farm;
    
      vector[N_plot] z_plot;
      real mu_plot;
      real<lower = 0> sigma_plot;
    
      real<lower = 0> scale;
    }
    
    transformed parameters{
      vector[N_plant] theta;
      vector[N_farm] alpha;
      vector[N_plot] tau;
      theta = mu_plant + z_plant * sigma_plant;
      alpha = mu_farm + z_farm * sigma_farm;
      tau = mu_plot + z_plot * sigma_plot;
      
    }
    
    model{
      vector[N] mu;
    
      mu_plant ~ normal(7, 2);
      z_plant ~ normal(0, 1);
      sigma_plant ~ exponential(1);
    
      mu_farm ~ normal(0, 1);
      z_farm ~ normal(0, 1);
      sigma_farm ~ exponential(1);
      
      mu_plot ~ normal(0, 1);
      z_plot ~ normal(0, 1);
      sigma_plot ~ exponential(1);
    
      scale ~ exponential(1);
    
      for (i in 1:N) {
        mu[i] = theta[plant_id[i]] + alpha[farm[i]] + 
                tau[plot[i]];
        mu[i] = exp(mu[i]);
      }
    
      total_fruts ~ neg_binomial_2(mu, scale);  
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] mu_1;
      array[N] int ppcheck;
    
      for (i in 1:N) {
        mu_1[i] = theta[plant_id[i]] + alpha[farm[i]] + 
                tau[plot[i]];
        mu_1[i] = exp(mu_1[i]);
      }
    
      for (i in 1:N) log_lik[i] = neg_binomial_2_lpmf(total_fruts[i] | mu_1[i], scale);
    
      ppcheck = neg_binomial_2_rng(mu_1, scale);
      
    }
    
    ")

file <- paste(getwd(), '/fruits_plant.stan', sep = '')

fit_fruit_plant <- cmdstan_model(file, compile = T)

mod_fruit_plant <- 
  fit_fruit_plant$sample(
    data = fruit_plant, 
    chains = 3, 
    parallel_chains = 3, 
    iter_sampling = 2e3, 
    iter_warmup = 500, 
    thin = 3, 
    seed = 123, 
    refresh = 200
  )

output_mod_tot_fru <- mod_fruit_plant$summary()

mod_diagnostics(mod_fruit_plant, output_mod_tot_fru)

trace_plot(mod_fruit_plant, output_mod_tot_fru$variable[1], 3)

par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))
for (i in 2:10) trace_plot(mod_fruit_plant, output_mod_tot_fru$variable[i], 3)
par(mfrow = c(1, 1))

post_tot_fruit <- mod_fruit_plant$draws(format = 'df')

post_tot_fruit <- 
  list(plant = post_tot_fruit[, grep('theta', colnames(post_tot_fruit))],
       farm = post_tot_fruit[, grep('alpha', colnames(post_tot_fruit))], 
       plot = post_tot_fruit[, grep('tau', colnames(post_tot_fruit))], 
       scale = post_tot_fruit$scale)

ppcheck_tot_fru <- mod_fruit_plant$draws(variables = 'ppcheck', format = 'matrix')

plot(density(ppcheck_tot_fru[1, ]), ylim = c(0, 0.00025), main = '', 
     xlab = 'Fruit produced per plant', lwd = 0.1)
for (i in 1:100) lines(density(ppcheck_tot_fru[i, ]), lwd = 0.1)
lines(density(fruit_plant$total_fruts), col = 'red', lwd = 2)


total_fruts <- as.vector(ppcheck_tot_fru)


# ======== 3. flowers produced by plant ========

set.seed(123)
total_flowers <- 
  total_fruts + 
  (total_fruts * 
     (1-sample(fruit_set, length(total_fruts), T)))

total_flowers <- round(total_flowers)

lines(density(total_fruts), col = 'red')

plot(density(total_flowers), main = '', 
     xlab = 'Flowers per blueberry bush', xlim = c(0, 20e3), 
     lwd = 4, col = 'lightblue') 


# ================== 4. Plants ha ======


p_01ha <- round((1501* # plantas of emerald in per ha (2011 in Citromax)
                   1)/0.45)


# =============== Simulating honeybee visits ========

simulated_visits <- 
  function(p_ha, # plants per ha
           flowers_plant, # total flowers per plant
           beta1 = 6, #flowering percentage (par 1 beta distribution) 
           beta2 = 4, #flowering percentage (par 2 beta distribution) 
           visits_bee, # visits per honeybee individual
           bees_hive, # number of hives and honeybees per hive
           hive_aggregate = F, # if T same plants are visited by all hives 
                               # if F each hive has its own plants
           short = F, # if TRUE return averaged values per plant instead of total 
                      # flowers
           seed = 123) {
    
    mins <- Sys.time()
    
    set.seed(seed)
    plants <- sample(flowers_plant, p_ha, replace = T) # plants per ha 
    
    set.seed(seed)
    flowering_perc <- rbeta(p_ha, beta1, beta2) # flowering percentage parm
    
    plants <- round(plants * flowering_perc) # flowering percentage
    
    
    plants <- lapply(plants, FUN = 
                       function(x) {
                         rep(0, x)
                       }) # number of flowers per plant 
    
    
    names(plants) <- paste('plant', 1:p_ha, sep = '')
    
    sum_flowers <- sum(unlist(lapply(plants, length), use.names = F))
    
    flowers <- rep(0, sum_flowers) # total flowers in 1 ha
    
    foragers_bees <- bees_hive # number of bees per hive
    
    foragers_visits <- lapply(foragers_bees, FUN = 
                                function(x) {
                                  sample(visits_bee, x, replace = T)
                                }) # trips per bee per hive
    
    if (hive_aggregate) {
      
      flowers <- 
        lapply(foragers_visits, FUN = 
                 function (hive) {
                   
                   
                   for (day in 1:2) { # number of pollination days (blueberry days receptivity)
                     
                     message(paste('Starting day', day, 'of pollination'))
                     
                     for (trips_bee in seq_along(hive)) { # trips per bee
                       
                       flower_id <- sample.int(sum_flowers, hive[trips_bee]) # visited flowers
                       
                       for (k in seq_along(flower_id)) flowers[flower_id[k]] <- flowers[flower_id[k]] + 1
                       
                     }
                     
                   }
                   
                   flowers
                 })
      
      plants <- 
        lapply(flowers, FUN = # iteration across flowers visited by bees from hive i
                 function(visits_hive_i) {
                   
                   lapply(plants, FUN = 
                            function(plant_i) {
                              
                              sample(visits_hive_i, length(plant_i)) # assign flower visits to each plant
                              
                            })
                   
                 })
      
      i <- 1
      while (i < length(plants)) {
        
        for (j in 1:p_ha) plants[[i+1]][[j]] <- plants[[i]][[j]] + plants[[i+1]][[j]]
        
        i <- i + 1
      }
      
    } else {
      
      flowers <- 
        lapply(foragers_visits, FUN = 
                 function (hive) {
                   
                   
                   for (day in 1:2) { # number of pollination days (blueberry days receptivity)
                     
                     message(paste('Starting day', day, 'of pollination'))
                     
                     for (trips_bee in seq_along(hive)) { # trips per bee
                       
                       flower_id <- sample.int(sum_flowers, hive[trips_bee]) # visited flowers
                       
                       for (k in seq_along(flower_id)) flowers[flower_id[k]] <- flowers[flower_id[k]] + 1
                       
                     }
                     
                   }
                   
                   flowers
                 })
      
      plants <- 
        lapply(flowers, FUN = # iteration across flowers visited by bees from hive i
                 function(visits_hive_i) {
                   
                   lapply(plants, FUN = 
                            function(plant_i) {
                              
                              sample(visits_hive_i, length(plant_i)) # asigns flower visits to each plant
                              
                            })
                   
                 })
      
    }
    
    names(plants) <- paste('Hive', 1:length(plants), sep = '')
    
    if (short) {
      
      message('Starting averagin visits across plants')
      
      plants_mu <- 
        lapply(plants, FUN = 
                 function(x) {
                   
                   t <- unlist(lapply(x, mean), use.names = F)
                   
                   tibble(plant = paste('plant', 1:length(t), sep = ''), 
                          visits = t)
                   
                 })
      
      for (i in seq_along(plants_mu)) {
        plants_mu[[i]]$n_hives <- paste(i)
      }
      
      plants_mu <- do.call('rbind', plants_mu)
      
      plants_mu <- 
        plants_mu |> 
        group_by(n_hives) |> 
        transmute(mu = median(visits), 
                  li = quantile(visits, 0.025),
                  ls = quantile(visits, 0.975)) |> 
        unique()
      
      message(paste('Toral execution time:'))
      print(Sys.time() - mins)
      return(plants_mu)
      
    } else {
      
      message(paste('Toral execution time:'))
      print(Sys.time() - mins)
      return(plants)
      
    }
    
  }

simulated_visits <- cmpfun(simulated_visits)



vis_LQ <- simulated_visits(p_ha = p_01ha,
                      flowers_plant = total_flowers, 
                      visits_bee = visits_day, 
                      bees_hive = hives_ha(20, mu_pop = 10e3), 
                      hive_aggregate = T)



vis_HQ <- simulated_visits(p_ha = p_01ha,
                         flowers_plant = total_flowers, 
                         visits_bee = visits_day_HQ, 
                         bees_hive = hives_ha(20, mu_pop = 20e3), 
                         hive_aggregate = T)


par(mfrow = c(2, 2), mar = c(4.2, 4.2, 1.5, 1.5))

plot(density(vis_LQ$Hive1[[1]]), ylim = c(0, 4), xlim = c(-0.5, 5),
     xlab = paste('Honeybee visits per flower at crop level\n (density:', 
                  1, ' hive per ha)'), 
     lwd = 0.3, main = '', col = 1)

for (j in 1:100) {
  lines(density(vis_LQ$Hive1[[j]]), col = 1, lwd = 0.3) 
}

for (j in 1:100) {
  lines(density(vis_HQ$Hive1[[j]]), col = 2, lwd = 0.3) 
}

text(x = c(4, 4), y = c(3, 3.25), c('Low quality', 'Hight quality'), 
     col = 1:2)


plot(density(vis_LQ$Hive7[[1]]), ylim = c(0, 0.45), xlim = c(-0.5, 20),
     xlab = paste('Honeybee visits per flower at crop level\n (density:', 
                  7, ' hive per ha)'), 
     lwd = 0.3, main = '', col = 1)

for (j in 1:100) {
  lines(density(vis_LQ$Hive7[[j]]), col = 1, lwd = 0.3) 
}

for (j in 1:100) {
  lines(density(vis_HQ$Hive7[[j]]), col = 2, lwd = 0.3) 
}

text(x = c(15, 15), y = c(0.35, 0.4), c('Low quality', 'Hight quality'), 
     col = 1:2)

plot(density(vis_LQ$Hive14[[1]]), ylim = c(0, 0.21), xlim = c(-0.5, 35),
     xlab = paste('Honeybee visits per flower at crop level\n (density:', 
                  14, ' hive per ha)'), 
     lwd = 0.3, main = '', col = 1)

for (j in 1:100) {
  lines(density(vis_LQ$Hive14[[j]]), col = 1, lwd = 0.3) 
}

for (j in 1:100) {
  lines(density(vis_HQ$Hive14[[j]]), col = 2, lwd = 0.3) 
}

text(x = c(28, 28), y = c(0.15, 0.17), c('Low quality', 'Hight quality'), 
     col = 1:2)

plot(density(vis_LQ$Hive20[[1]]), ylim = c(0, 0.15), xlim = c(-0.5, 50),
     xlab = paste('Honeybee visits per flower at crop level\n (density:', 
                  20, ' hive per ha)'), 
     lwd = 0.3, main = '', col = 1)

for (j in 1:100) {
  lines(density(vis_LQ$Hive20[[j]]), col = 1, lwd = 0.3) 
}

for (j in 1:100) {
  lines(density(vis_HQ$Hive20[[j]]), col = 2, lwd = 0.3) 
}

text(x = c(35, 35), y = c(0.11, 0.12), c('Low quality', 'Hight quality'), 
     col = 1:2)


par(mfrow = c(1, 1))


cluster <- makeCluster(detectCores() - 1)

clusterExport(cluster, c('simulated_visits', 'total_flowers', 
                         'visits_day', 'hives_ha', 'p_01ha', 'visits_day_HQ'))

clusterEvalQ(cluster, {
  pks <- c('tidyverse', 'rethinking', 'magrittr', 'cmdstanr', 'parallel')
  
  sapply(pks, library, character.only = T)
})

t1 <- Sys.time()
p <- parLapply(cluster, 1:100, fun = 
                 function(i) {
                   print(paste('Starting sim', i))
                   x <- simulated_visits(p_ha = p_01ha,
                                         flowers_plant = total_flowers, 
                                         visits_bee = visits_day, 
                                         bees_hive = hives_ha(20, seed = i+500), 
                                         hive_aggregate = T, 
                                         short = T) 
                   x$sim <- paste('sim', i, sep = '')
                   x
                 })
Sys.time() - t1

p <- do.call('rbind', p) 

t1 <- Sys.time()
p_HQ <- parLapply(cluster, 1:100, fun = 
                 function(i) {
                   
                   x <- simulated_visits(p_ha = p_01ha,
                                         flowers_plant = total_flowers, 
                                         visits_bee = visits_day_HQ, 
                                         bees_hive = hives_ha(20, mu_pop = 20e3,
                                                              seed = i+500), 
                                         hive_aggregate = T, 
                                         short = T) 
                   x$sim <- paste('sim', i, sep = '')
                   x
                 })
Sys.time() - t1

stopCluster(cluster)
rm(list = 'cluster')

p_HQ <- do.call('rbind', p_HQ)

p$quality <- 'low'
p_HQ$quality <- 'hight'

vis_hives <- rbind(p, p_HQ)

unique(vis_hives$sim)

vis_hives |> 
  ggplot(aes(as.numeric(n_hives), mu, shape = sim,
             ymin = li, ymax = ls, color = quality)) +
  geom_line(alpha = 0.3) + #geom_ribbon(alpha = 0.2) +
  labs(x = 'Hives per blueberry ha', 
       y = 'Average flower visits per flower\n at crop level') +
  scale_shape_manual(values = rep(1, 100)) +
  scale_color_manual(values = c('lightblue3', 'tan1')) +
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid = element_blank())



# ====== 5. Hoenybee pollen deposition ====
apisSVP <- readRDS('honeybee_svpd_data.rds')

apisSVP <- apisSVP[-which.max(apisSVP)] # removing value with 155 pollen deposition

plot(density(rnbinom(1e4, size = 2, mu = 15)), 
     main = expression(mu['honeybee pollen deposition'~'~'~'NegBinom(size = 2, mu = 15)']), 
     xlab = 'Single visit pollen deposition')

dat_apis <- list(pollen = apisSVP, 
                 N = length(apisSVP), 
                 ind_id = 1:length(apisSVP), 
                 N_ind = length(apisSVP))

cat(file = 'model_apis_SVP.stan', 
    '
    data{
      int N;
      int N_ind;
      array[N] int pollen;
      array[N] int ind_id;
    }
    
    parameters{
      vector[N_ind] z_alpha;
      real mu;
      real<lower = 0> sigma;
      real<lower = 0> phi;
    }
    
    transformed parameters{
      vector[N_ind] alpha;
      alpha = mu + z_alpha * sigma;
    }
    
    model {
      vector[N] lambda;
      z_alpha ~ normal(0, 1);
      mu ~ normal(3, 0.5);
      sigma ~ exponential(1);
      phi ~ exponential(1);
      
      for (i in 1:N){
        lambda[i] = alpha[ind_id[i]];
        lambda[i] = exp(lambda[i]);
      }
      
      pollen ~ neg_binomial_2(lambda, phi);
    }
    
    generated quantities{
      vector[N] log_lik;
      array[N] int ppcheck;
      vector[N] lambda1;
    
      for (i in 1:N){
        lambda1[i] = alpha[ind_id[i]];
        lambda1[i] = exp(lambda1[i]);
      }
    
      for (i in 1:N) log_lik[i] = neg_binomial_2_lpmf(pollen[i] | lambda1[i], phi);
    
      ppcheck = neg_binomial_2_rng(lambda1, phi);
      
    }
    ')

file <- paste(getwd(), '/model_apis_SVP.stan', sep = '')
fit_apis_svp <- cmdstan_model(file, compile = T)

model_apis_svp <- 
  fit_apis_svp$sample(
    data = dat_apis, 
    chains = 3, 
    iter_sampling = 10e3, 
    iter_warmup = 500, 
    parallel_chains = 3, 
    thin = 3, 
    refresh = 200, 
    seed = 123
  )

out_put_apis_svp <- model_apis_svp$summary()

mod_diagnostics(model_apis_svp, out_put_apis_svp)

model_apis_svp$summary()

par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))
for (i in 1:9) trace_plot(model_apis_svp, out_put_apis_svp$variable[i], 3)
par(mfrow = c(1, 1))

post_apis_svp <- model_apis_svp$draws(format = 'df')[-1]
post_apis_svp$mu <- exp(post_apis_svp$mu)

pp_check_svp <- model_apis_svp$draws(variables = 'ppcheck', format = 'matrix')
  
plot(density(pp_check_svp[1, ]), lwd = 0.1, main = '', 
     xlab = 'Single visit pollen deposition', ylim = c(0, 0.05))
for (i in 1:100) lines(density(pp_check_svp[i, ]), lwd = 0.1)
lines(density(dat_apis$pollen), col = 'red', lwd = 3)

set.seed(1234)
sim_visits <- 
  apply(pp_check_svp[c(sample(1:nrow(pp_check_svp), size = 100, replace = F)),], 1, 
        simplify = 'list', FUN = 
          function(x) {
            
            df_2 <- tibble()
            
            for (j in 1:10) {
              
              sum_vis <- 
                sapply(1:2e3, FUN = 
                         function(k) {
                           sum(sample(x, size = j, T))
                         })
              
              df_ <- tibble(poll_ac = sum_vis,
                            vis = j)
              
              df_2 <- rbind(df_2, df_)
              
            }
            df_2[sample(1:nrow(df_2), size = 200, replace = F), ]
          })

names(sim_visits) <- paste('sim', 1:length(sim_visits))

for (i in seq_along(sim_visits)) {
  sim_visits[[i]]$sim <- rep(paste('sim', i), nrow(sim_visits[[i]]))
}

sim_visits <- do.call('rbind', sim_visits)

sim_visits %$% plot(vis + rnorm(1e4, 0.1, 0.1), poll_ac, cex = 0.1, 
                    xlab = 'Number of honeybee visits', 
                    ylab = 'Pollen deposition', 
                    ylim = c(0, 600))
points(1:10, 
       tapply(sim_visits$poll_ac, sim_visits$vis, FUN = mean), 
       col = 'red', pch = 15)
abline(h = c(112, 274), lty = 2, col = 'red')
text(x = 2, y = 300, 'Optimal pollination')

beta <- rnorm(20, 3, 0.5)
alpha <- rnorm(20, 20, 5)

sim_visits <- lapply(sim_visits[1:2], function(x) x)

sim_visits$N <- length(sim_visits$poll_ac)

cat(file = 'slope_svp_apis.stan', 
    '
    data{
      int N;
      array[N] int poll_ac;
      array[N] int vis;
    }
    
    parameters{
      real alpha;
      real beta;
      real<lower = 0> sigma;
    }
    
    model{
      vector[N] lambda;
      alpha ~ normal(20, 5);
      beta ~ normal(3, 0.5);
      sigma ~ exponential(1);
    
      for (i in 1:N) {
        lambda[i] = alpha + beta*vis[i];
        lambda[i] = exp(lambda[i]);
      }
    
      poll_ac ~ neg_binomial_2(lambda, sigma);
    }
    
    generated quantities{
      vector[N] lambda1;
      array[N] int ppcheck;
      vector[N] log_lik;
     
      for (i in 1:N) {
        lambda1[i] = alpha + beta*vis[i];
        lambda1[i] = exp(lambda1[i]);
      }
    
      for (i in 1:N) log_lik[i] = neg_binomial_2_lpmf(poll_ac[i] | lambda1[i], sigma);
    
      ppcheck = neg_binomial_2_rng(lambda1, sigma);
    }
    
    ')

file <- paste(getwd(), '/slope_svp_apis.stan', sep = '')

fit_svp_slope <- cmdstan_model(file, compile = T)

mod_svp_slope <- 
  fit_svp_slope$sample(
    data = sim_visits, 
    chains = 3,
    parallel_chains = 3, 
    iter_sampling = 2e3, 
    iter_warmup = 500, 
    thin = 3,
    refresh = 200, 
    seed = 123
  )

# output_mod_slope <- mod_svp_slope$summary()
# 
# mod_diagnostics(mod_svp_slope, output_mod_slope)

post_svp_slope <- mod_svp_slope$draws(variables = c('alpha', 'beta'), format = 'df')
colnames(post_svp_slope)
#for (i in 2:3) post_svp_slope[[i]] <- exp(post_svp_slope[[i]])

# par(mfrow = c(2, 2))
# for (i in 1:4) trace_plot(mod_svp_slope, output_mod_slope$variable[i], 3)
# par(mfrow = c(1, 1))

ppcheck_slope <- 
  sapply(1:length(sim_visits$vis), FUN = 
           function(i) {
             
             x <- sim_visits$vis[i]
             
             lambda <- 
               post_svp_slope$alpha +
               post_svp_slope$beta*x
             
             rnbinom(100, 
                     mu = exp(lambda),
                     size = exp(post_svp_slope$alpha))
           })
ppcheck_slope2 <- mod_svp_slope$draws(variables = 'ppcheck', format = 'matrix')

plot(NULL, lwd = 0.1, main = '', 
     xlab = 'Single visit pollen deposition', 
     xlim = c(0, 600), ylim = c(0, 0.007), ylab = 'Density')
for (i in 1:50) lines(density(ppcheck_slope[i, ]), lwd = 0.1)
for (i in 1:50) lines(density(ppcheck_slope2[i, ]), lwd = 0.1, col = 'lightblue3')
lines(density(sim_visits$poll_ac), col = 'red')

rm(list = c('ppcheck_slope2'))
# here I have to fit the poisson model and then 


# ========= 6. Visit to pollen ========

pollen_fruit <- readRDS('datos_experimento.rds')
pollen_fruit <- pollen_fruit[, c('planta', "tratamiento", "fruto_diam", 
                                 "carga_poli", "carga_poli2")]
pollen_surpass <- read_xlsx('D:/github_repos/cap4_phd/all_data.xlsx', 
                      sheet = 7, na = 'NA')[, -2]

quantile(pollen_fruit[pollen_fruit$tratamiento == 'l', ]$carga_poli, 
         probs = seq(0.1, 1, by = 0.1))

dat_asymptote_pollen <- as_tibble(pollen_fruit[pollen_fruit$tratamiento == 'l', ])

dat_asymptote_pollen$farm <- 'sta_lu'
dat_asymptote_pollen$plant_id <- dat_asymptote_pollen %$% paste(farm, planta, sep = '')

pollen_surpass$plant <- pollen_surpass %$% paste(farm, plant, sep = '')
pollen_surpass$cultivar <- 'eme'

dat_asymptote_pollen <- 
  dat_asymptote_pollen |>
  select(farm, plant_id, carga_poli) |> 
  mutate(cultivar = 'sch')

colnames(dat_asymptote_pollen) <- colnames(pollen_surpass)

dat_asymptote_pollen <- rbind(dat_asymptote_pollen, 
                              pollen_surpass[!is.na(pollen_surpass$no_polen),])

levels(as.factor(dat_asymptote_pollen$cultivar))

dat_asymptote_pollen <- 
  lapply(dat_asymptote_pollen, function(x) if(!is.double(x)) as.integer(as.factor(x)) else(x))

dat_asymptote_pollen$N <- length(dat_asymptote_pollen$farm)
dat_asymptote_pollen$N_cultivar <- 2
dat_asymptote_pollen$N_plant <- max(dat_asymptote_pollen$plant)
dat_asymptote_pollen$N_farm <- max(dat_asymptote_pollen$farm)

plot(density(rnbinom(1e3, size = 2, mu = exp(5))))
cat(file = 'asymptote_pollen1.stan', 
    '
    data{
      int N;
      int N_cultivar;
      int N_plant;
      int N_farm;
      array[N] int farm;
      array[N] int plant;
      array[N] int no_polen;
      array[N] int cultivar;
    }
    
    parameters{
      vector[N_plant] z_alpha;
      real mu_alpha;
      real<lower = 0> sigma_alpha;
    
      vector[N_farm] z_theta;
      real mu_theta;
      real<lower = 0> sigma_theta;
    
      vector[N_cultivar] z_tau;
      real mu_tau;
      real<lower = 0> sigma_tau;
    
      real<lower = 0> scale;
    }
    
    transformed parameters{
      vector[N_plant] alpha;
      vector[N_farm] theta;
      vector[N_cultivar] tau;
    
      alpha = mu_alpha + z_alpha*sigma_alpha;
      tau = mu_tau + z_tau * sigma_tau;
      theta = mu_theta + z_theta * sigma_theta;
    }
    
    model{
      vector[N] mu;
      z_alpha ~ normal(0, 1);
      mu_alpha ~ normal(0, 1);
      sigma_alpha ~ exponential(1);
    
      z_theta ~ normal(0, 1);
      mu_theta ~ normal(0, 1);
      sigma_theta ~ exponential(1);
    
      z_tau ~ normal(0, 1);
      mu_tau ~ normal(5, 0.5);
      sigma_tau ~ exponential(1);
    
      scale ~ exponential(1);
    
      for (i in 1:N){
        mu[i] = tau[cultivar[i]] + theta[farm[i]] + alpha[plant[i]];
        mu[i] = exp(mu[i]);
      }
    
      no_polen ~ neg_binomial_2(mu, scale);
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] mu1;
      array[N] int ppcheck;
    
      for (i in 1:N){
        mu1[i] = tau[cultivar[i]] + theta[farm[i]] + alpha[plant[i]];
        mu1[i] = exp(mu1[i]);
      }
    
      for (i in 1:N) log_lik[i] = neg_binomial_2_lpmf(no_polen[i] | mu1[i], scale);
      
      ppcheck = neg_binomial_2_rng(mu1, scale);
    }
    ')


file <- paste(getwd(), '/asymptote_pollen1.stan', sep = '')
fit_asymtote_pollen <- cmdstan_model(file, compile = T)

mod_asymptote_pollen <- 
  fit_asymtote_pollen$sample(
    data = dat_asymptote_pollen, 
    chains = 3, 
    parallel_chains = 3, 
    iter_sampling = 4e3,
    iter_warmup = 500, 
    thin = 3, 
    refresh = 500, 
    seed = 123
  )

out_asymptote_pollen <- mod_asymptote_pollen$summary()

mod_diagnostics(mod_asymptote_pollen, out_asymptote_pollen)

ppcheck_asymptote_pol <- mod_asymptote_pollen$draws(variables = 
                                                      'ppcheck', 
                                                    format = 'matrix')

plot(density(ppcheck_asymptote_pol[1, ]), ylab = 'Density', 
     xlab = 'Pollen load', main = '', ylim = c(0, 0.009))
for (i in 1:100) lines(density(ppcheck_asymptote_pol[i, ], lwd = 0.1))
lines(density(dat_asymptote_pollen$no_polen), col = 'red', lwd = 2)

post_asymptote <- 
  mod_asymptote_pollen$draws(variables = 
                               c('alpha', 'theta', 'tau', 'scale'), 
                             format = 'df')

post_asymptote <- 
  list(plant = post_asymptote[, grepl('alpha', colnames(post_asymptote)), ], 
       farm = post_asymptote[, grepl('theta', colnames(post_asymptote)), ], 
       cultivar = post_asymptote[, grepl('tau', colnames(post_asymptote)), ],
       scale = post_asymptote[, grepl('scale', colnames(post_asymptote)), ])



sch <- with(post_asymptote, {
                 cultivar[, 2, drop = T] +
                   apply(plant, 1, mean) + 
                   apply(farm, 1, mean)
               })
  
sch <- replicate(100, 
                 rnbinom(1000, size = post_asymptote$scale$scale, 
                         mu = exp(sch)))

eme <- with(post_asymptote, {
                 cultivar[, 1, drop = T] +
                 apply(plant, 1, mean) + 
                 apply(farm, 1, mean)
              })

eme <- replicate(100, 
                 rnbinom(1000, size = post_asymptote$scale$scale, 
                         mu = exp(eme)))

plot(NULL, xlim = c(0, 450), ylim = c(0, 1), 
     xlab = 'Pollen load', ylab = 'eCDF')
for (i in 1:100) lines(ecdf(sch[, i]), col = 'lightblue')
for (i in 1:100) lines(ecdf(eme[, i]), col = 'tan')
abline(v = c(225, 250, 275), lty = c(2, 3, 2), col = 'red') # 

mu_asintota <- 250
asymptote <- rnorm(length(post_svp_slope$beta), mu_asintota, 15)
slope <- post_svp_slope$beta

plot(NULL, xlim = c(0, 30), ylim = c(0, 320), xlab = 'Honeybee visits', 
     ylab = 'Single visit pollen deposition')
for (i in 1:200) curve(asymptote[i]*(1-exp(-slope[i]*x)), 
                       add = T, lwd = 0.1)
curve(mean(asymptote)*(1-exp(-mean(slope)*x)), 
      add = T, lwd = 2, col = 'red')

pollen_deposition_fun <- function(x, 
                                  mu_est = T, 
                                  beta = slope, 
                                  theta = asymptote) {
  if (mu_est) {
    a <- mean(theta)
    b <- mean(beta)
    
    return(round(a * (1-exp(-b*x))))
  } else {
    a <- sample(theta, 1, T)
    b <- sample(beta, 1, T)
    
    return(round(a * (1-exp(-b*x))))
  }
}

pollen_deposition_fun <- cmpfun(pollen_deposition_fun)

pollen_deposition_LQ <- 
  lapply(vis_LQ, FUN = 
           function(x) {
             lapply(x, pollen_deposition_fun, mu_est = F)
           })

pollen_deposition_HQ <- 
  lapply(vis_HQ, FUN = 
           function(x) {
             lapply(x, pollen_deposition_fun, mu_est = F)
           })


par(mfrow = c(2, 2))

plot(density(pollen_deposition_LQ$Hive1[[1]]), ylim = c(0, 0.06), 
     xlim = c(-20, 170),
     xlab = paste('Pollen deposition per flower at farm level', 
                  '(1 hive per ha)'), 
     lwd = 0.3, main = '', col = 1)

for (j in 1:100) {
  lines(density(pollen_deposition_LQ$Hive1[[j]]), col = 1, lwd = 0.3) 
}

for (j in 1:100) {
  lines(density(pollen_deposition_HQ$Hive1[[j]]), col = 2, lwd = 0.3) 
}


plot(density(pollen_deposition_LQ$Hive7[[1]]), ylim = c(0, 0.025), 
     xlim = c(-10, 300),
     xlab = paste('Pollen deposition per flower at farm level', 
                  '(7 hives per ha)'), 
     lwd = 0.3, main = '', col = 1)

for (j in 1:100) {
  lines(density(pollen_deposition_LQ$Hive7[[j]]), col = 1, lwd = 0.3) 
}

for (j in 1:100) {
  lines(density(pollen_deposition_HQ$Hive7[[j]]), col = 2, lwd = 0.3) 
}


plot(density(pollen_deposition_LQ$Hive14[[1]]), ylim = c(0, 0.025), 
     xlim = c(-10, 300),
     xlab = paste('Pollen deposition per flower', 
                  '(14 hives per ha)'), 
     lwd = 0.3, main = '', col = 1)

for (j in 1:100) {
  lines(density(pollen_deposition_LQ$Hive14[[j]]), col = 1, lwd = 0.3) 
}

for (j in 1:100) {
  lines(density(pollen_deposition_HQ$Hive14[[j]]), col = 2, lwd = 0.3) 
}


plot(density(pollen_deposition_LQ$Hive20[[1]]), ylim = c(0, 0.05), 
     xlim = c(-10, 300),
     xlab = paste('Pollen deposition per flower', 
                  '(20 hives per ha)'), 
     lwd = 0.4, main = '', col = 1)

for (j in 1:100) {
  lines(density(pollen_deposition_LQ$Hive20[[j]]), col = 1, lwd = 0.3) 
}

for (j in 1:100) {
  lines(density(pollen_deposition_HQ$Hive20[[j]]), col = 2, lwd = 0.3) 
}


par(mfrow = c(1, 1))

rm(list = c('pollen_deposition_HQ', 'pollen_deposition_LQ',
            'vis_HQ', 'vis_LQ'))

crop_pollination <- function(p_ha, # plants per ha
                             flowers_plant, # total flowers per plant
                             beta1 = 6, #flowering percentage (par 1 beta distribution) 
                             beta2 = 4, #flowering percentage (par 2 beta distribution) 
                             visits_bee, # visits per honeybee individual
                             bees_hive, # number of hives and honeybees per hive
                             hive_aggregate = T, # if T same plants are visited by all hives 
                             # if F each hive has its own plants
                             # flowers
                             short = F,
                             seed = 123) {
  
  message('Starting floral visits')
  t1 <- Sys.time()
  
  visits <- 
    simulated_visits(p_ha,
                     flowers_plant,
                     beta1,
                     beta2,
                     visits_bee,
                     bees_hive,
                     hive_aggregate,
                     short = F, 
                     seed)
  
  message('Starting pollen deposition')
  PD <- 
    lapply(visits, FUN = 
             function(x) {
               lapply(x, pollen_deposition_fun, mu_est = F)
             })
  
  if (short == F) {
    
    PD
    
  } else {
    
    PD <- 
      lapply(PD, FUN = 
               function(x) {
                 
                 t <- unlist(lapply(x, mean), use.names = F)
                 
                 tibble(plant = paste('plant', 1:length(t), sep = ''), 
                        pollen = t)
                 
               })
    
    for (i in seq_along(PD)) {
      PD[[i]]$n_hives <- paste(i)
    }
    
    PD <- do.call('rbind', PD)
    
    PD <- 
      PD |> 
      group_by(n_hives) |> 
      transmute(mu = median(pollen), 
                li = quantile(pollen, 0.025),
                ls = quantile(pollen, 0.975)) |> 
      unique()
    
    message('Total execution time:')
    print(Sys.time() - t1)
    return(PD)
    
  }
  
}

crop_pollination <- cmpfun(crop_pollination)

cluster <- makeCluster(8)

clusterExport(cluster, c('simulated_visits', 'pollen_deposition_fun', 
                         'total_flowers', 'visits_day', 'hives_ha', 
                         'p_01ha', 'visits_day_HQ', 'crop_pollination', 
                         'asymptote', 'slope'))

clusterEvalQ(cluster, {
  pks <- c('tidyverse', 'magrittr', 'cmdstanr',
           'ggdag', 'dagitty', 'parallel')
  
  sapply(pks, library, character.only = T)
})

t <- Sys.time()
pollen_LQ <- 
  parLapply(cluster, 1:100, fun = 
              function(i) {
                
                x <- crop_pollination(p_ha = p_01ha,
                                      flowers_plant = total_flowers, 
                                      visits_bee = visits_day, 
                                      bees_hive = hives_ha(20, seed = i+500), 
                                      hive_aggregate = T, 
                                      short = T)
                
                x$sim <- paste('sim', i, sep = '')
                x
              })
Sys.time() - t

pollen_LQ <- do.call('rbind', pollen_LQ) 


t <- Sys.time()
pollen_HQ <- 
  parLapply(cluster, 1:100, fun = 
              function(i) {
                
                x <- crop_pollination(p_ha = p_01ha,
                                      flowers_plant = total_flowers, 
                                      visits_bee = visits_day_HQ, 
                                      bees_hive = hives_ha(20, mu_pop = 20e3, 
                                                           seed = i+500), 
                                      hive_aggregate = T, 
                                      short = T)
                x$sim <- paste('sim', i, sep = '')
                x
              })
Sys.time() - t

stopCluster(cluster)
rm(list = 'cluster')

pollen_HQ <- do.call('rbind', pollen_HQ)



pollen_LQ$quality <- 'low'
pollen_HQ$quality <- 'hight'

pollen_hives <- rbind(pollen_LQ, pollen_HQ)

plot_vis_hive <- 
  vis_hives |> 
  ggplot(aes(as.numeric(n_hives), mu, shape = sim,
             ymin = li, ymax = ls, color = quality)) +
  geom_line(alpha = 1, linewidth = 0.15) + #geom_ribbon(alpha = 0.2) +
  labs(x = 'Hives per blueberry ha', 
       y = 'Average flower visits per flower\n at crop level') +
  scale_shape_manual(values = rep(1, 100)) +
  scale_color_manual(values = c('lightblue3', 'tan1')) +
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid = element_blank())

plot_pollen_hive <- 
  pollen_hives |> 
  ggplot(aes(as.numeric(n_hives), mu, shape = sim,
             ymin = li, ymax = ls, color = quality)) +
  geom_line(linewidth = 0.15) + #geom_ribbon(alpha = 0.2) +
  labs(x = 'Hives per blueberry ha', 
       y = 'Average pollen deposition per flower\n at crop level') +
  scale_shape_manual(values = rep(1, 100)) +
  scale_color_manual(values = c('lightblue3', 'tan1')) +
  geom_hline(yintercept = c(112, 274), linetype = 2, color = 'red') +
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid = element_blank())

plot_grid(plot_vis_hive, plot_pollen_hive, ncol = 2)

ggsave('simulation_visit_pollen.jpg', width = 16, height = 10, units = 'cm', dpi = 700)

# ============= 7. Pollen to fruit =====

# ===== Contrasts fruit size cultivars =====

snow <- as_tibble(readRDS('datos_experimento.rds'))
snow <- snow[, c("planta", "tratamiento", "carga_poli", "fruto_diam")]
snow$farm <- 'sta_lu'
snow <- snow[snow$tratamiento == 'l', ]
snow$plant_id <- snow %$% paste(planta, tratamiento, sep = '')
snow$variedad <- 'sch'
snow1 <- snow[, c("farm", "plant_id", "variedad", "fruto_diam")]
colnames(snow1) <- c('farm', 'plant_id', 'cultivar', 'fruit_diam')

emerald <- readRDS('fruit_size.rds')
emerald <- emerald[, c("year", "farm", "plant", 
                       "fruit_diameter", "fruit_weight")]
emerald$plant_id <- emerald %$% paste(plant, farm, sep = '')
for (i in 1:2) emerald[[i]] <- as.factor(emerald[[i]])
emerald$variedad <- 'eme'
emerald1 <- emerald[, c('farm', "plant_id", "variedad", "fruit_diameter")]
colnames(emerald1) <- colnames(snow1)

pri_sj <- as_tibble(read.csv2('calidad_frutoPRI_SJ.csv', header = T, sep = ';'))
pri_sj$plant_id <- pri_sj %$% paste(finca, planta, sep = '')
pri_sj1 <- pri_sj[, c("finca", "plant_id", "variedad", "diamF")]
colnames(pri_sj1) <- colnames(snow1)

fruit_size_cultivars <- rbind(snow1, emerald1, pri_sj1)

fruit_size_cultivars[] <- 
  lapply(fruit_size_cultivars, function(x) if(!is.numeric(x)) as.factor(x) else(x))

fruit_size_cultivars <- 
  fruit_size_cultivars[fruit_size_cultivars$fruit_diam > 2, ]

dat_fru_size <- 
  lapply(fruit_size_cultivars[!is.na(fruit_size_cultivars$fruit_diam), ], 
         function(x) if(!is.numeric(x)) as.integer(x) else(x))

lapply(dat_fru_size, function(x) sum(is.na(x)))

dat_fru_size$N <- length(dat_fru_size$farm)
dat_fru_size$N_plants <- max(dat_fru_size$plant_id)
dat_fru_size$N_cultivars <- max(dat_fru_size$cultivar)
dat_fru_size$N_farms <- max(dat_fru_size$farm)

cat(file = 'fruit_size_comp.stan', 
    '
    data{
      int N;
      int N_cultivars;
      int N_plants;
      int N_farms;
      vector[N] fruit_diam;
      array[N] int cultivar;
      array[N] int plant_id;
      array[N] int farm;
    }
    
    parameters{
      vector[N_cultivars] z_alpha;
      real mu_alpha;
      real<lower = 0> sigma_alpha;
    
      vector[N_plants] z_theta;
      real mu_theta;
      real<lower = 0> sigma_theta;
    
      vector[N_farms] z_tau;
      real mu_tau;
      real<lower = 0> sigma_tau;
    
      real<lower = 0> sigma;
    }
    
    transformed parameters{
      vector[N_cultivars] alpha;
      vector[N_plants] theta;
      vector[N_farms] tau;
    
      alpha = mu_alpha + z_alpha * sigma_alpha;
      theta = mu_theta + z_theta * sigma_theta;
      tau = mu_tau + z_tau * sigma_tau;
    }
    
    model{
      vector[N] mu;
      sigma ~ exponential(1);
    
      z_alpha ~ normal(0, 1);
      mu_alpha ~ normal(15, 2.5);
      sigma_alpha ~ exponential(1);
    
      z_theta ~ normal(0, 1);
      mu_theta ~ normal(0, 1);
      sigma_theta ~ exponential(1);
    
      z_tau ~ normal(0, 1);
      mu_tau ~ normal(0, 1);
      sigma_tau ~ exponential(1);
    
      for (i in 1:N){
        mu[i] = alpha[cultivar[i]] + theta[plant_id[i]] + tau[farm[i]];
      }
    
      fruit_diam ~ normal(mu, sigma);
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] mu1;
      array[N] real ppcheck;
    
      for (i in 1:N){
        mu1[i] = alpha[cultivar[i]] + theta[plant_id[i]] + tau[farm[i]];
      }
    
      for (i in 1:N) log_lik[i] = normal_lpdf(fruit_diam[i] | mu1[i], sigma);
    
      ppcheck = normal_rng(mu1, sigma);
    }
    ')

file <- paste(getwd(), '/fruit_size_comp.stan', sep = '')
fit_fruit_size <- cmdstan_model(file, compile = T)

mod_fruit_size <- 
  fit_fruit_size$sample(
    data = dat_fru_size,
    iter_sampling = 4e3,
    iter_warmup = 500,
    chains = 3,
    parallel_chains = 3,
    thin = 3,
    refresh = 500,
    seed = 123
  )

out_fruit_size <- mod_fruit_size$summary()

mod_diagnostics(mod_fruit_size, out_fruit_size)

ppcheck_fruit_size <- mod_fruit_size$draws('ppcheck', format = 'matrix')

plot(density(ppcheck_fruit_size[1, ]), xlab = 'fruit size (mm)', 
     ylab = 'Density', main = '', lwd = 0.1, ylim = c(0, 0.25))
for (i in 1:500) lines(density(ppcheck_fruit_size[i, ]), lwd = 0.1)
lines(density(dat_fru_size$fruit_diam), col = 'red', lwd = 1.5)

post_fruit_size <- mod_fruit_size$draws(c('alpha', 'theta', 'tau', 'sigma'), 
                                        format = 'df')

post_fruit_size <- 
  list(
    cultivar = post_fruit_size[, grepl('alpha', colnames(post_fruit_size))],
    plant = post_fruit_size[, grepl('theta', colnames(post_fruit_size))],
    farm = post_fruit_size[, grepl('tau', colnames(post_fruit_size))],
    sigma = post_fruit_size[, grepl('sigma', colnames(post_fruit_size))]
  )

levels(fruit_size_cultivars$cultivar)

fruits_sim <- 
  lapply(1:4, FUN = 
           function(x) {
             
             cult <- levels(fruit_size_cultivars$cultivar)[x]
             
             mu <- 
               with(post_fruit_size, 
                    {
                      cultivar[, x, drop = T] +
                        apply(plant, 1, mean) +
                        apply(farm, 1, mean)
                    })
             replicate(1e3, rnorm(1e3, mu, post_fruit_size$sigma$sigma))
           })

names(fruits_sim) <- levels(fruit_size_cultivars$cultivar)

fruits_sim <- lapply(c(1, 2, 4), FUN = 
                       function(x) {
                         
                         sapply(1:100, FUN = 
                                  function(i) {
                                    
                                    fruits_sim[[x]][, i] - 
                                      fruits_sim[[3]][, i]
                                    
                                  }, simplify = 'array')
                       })

names(fruits_sim) <- levels(fruit_size_cultivars$cultivar)[-3]

par(mfrow = c(2, 2), mar = c(4.2, 4.2, 1, 1))

for (i in 1:3) {
  
  plot(density(fruits_sim[[i]][, 1]), lwd = 0.1, 
       xlab = paste('contrast sch - ', names(fruits_sim)[i]), 
       main = '', col = i, ylim = c(0, 0.25))
  
  for (j in 1:100) lines(density(fruits_sim[[i]][, j]), col = i, lwd = 0.1)
  
}
par(mfrow = c(1, 1))

fruits_sim <- 
  lapply(fruits_sim, FUN = 
           function(x) apply(x, 2, mean))

plot(NULL, lwd = 2, ylim = c(0, 8), xlim = c(-0.5, 2.2),
     xlab = paste('contrast sch - cultivars'), 
     main = '', ylab = 'Density')
for (i in 1:3) lines(density(fruits_sim[[i]]), lwd = 2, col = i)

# ======= pollen fruit function ======

snow <- readRDS('datos_experimento.rds')

colnames(snow)[10] <- 'fruit_sch'
snow$fruit_eme <- snow$fruit_sch + sample(fruits_sim$eme, nrow(snow), T)
snow$fruit_sj <- snow$fruit_sch + sample(fruits_sim$sj, nrow(snow), T)
snow$fruit_pri <- snow$fruit_sch + sample(fruits_sim$pri, nrow(snow), T)

snow <- 
  as_tibble(snow[, c(6, grep('fruit', colnames(snow)))])

par(mfrow = c(2, 2), mar = c(4.2, 4.2, 1, 1))
for (i in 2:5) plot(snow$carga_poli, snow[, i, drop = T], ylab = 'Fruit size (mm)', 
                    xlab = 'Pollen deposition', col = i, main = colnames(snow)[i])
par(mfrow = c(1, 1))

snow$carga_poli_z <- as.vector(scale(snow$carga_poli))

dat_pollen_fun <- lapply(snow, function(x) x)
dat_pollen_fun$N <- nrow(snow)
dat_pollen_fun$carga_poli2 <- dat_pollen_fun$carga_poli_z^2

x <- dat_pollen_fun %$% seq(min(carga_poli), max(carga_poli), by = 1)
a <- 17
b <- 100
m <- 0.05
plot(x, a / (1 + b*exp(-x*m)), type = 'l')

cat(file = 'pollen_function.stan', 
    '
    data{
      int N;
      vector[N] carga_poli;
      vector[N] fruit_sch;
      vector[N] fruit_eme;
      vector[N] fruit_pri;
      vector[N] fruit_sj;
    }
    
    parameters{
      real alpha_sch;
      real phi1_sch;
      real phi2_sch;
      real<lower = 0> sigma_sch;
    
      real alpha_pri;
      real phi1_pri;
      real phi2_pri;
      real<lower = 0> sigma_pri;
    
      real alpha_eme;
      real phi1_eme;
      real phi2_eme;
      real<lower = 0> sigma_eme;
    
      real alpha_sj;
      real phi1_sj;
      real phi2_sj;
      real<lower = 0> sigma_sj;
    }
    
    model{
    
      // model snowchaser
      vector[N] mu_sch;
      alpha_sch ~ normal(17, 2);
      phi1_sch ~ normal(0, 1);
      phi2_sch ~ normal(0, 1);
      sigma_sch ~ exponential(1);
      
      for (i in 1:N) { 
        mu_sch[i] = alpha_sch/(1+exp(-(phi1_sch + phi2_sch * carga_poli[i])));
      }
    
      fruit_sch ~ normal(mu_sch, sigma_sch);
    
      // model emerald
      vector[N] mu_eme;
      alpha_eme ~ normal(17, 2);
      phi1_eme ~ normal(0, 1);
      phi2_eme ~ normal(0, 1);
      sigma_eme ~ exponential(1);
      
      for (i in 1:N) { 
        mu_eme[i] = alpha_eme/(1+exp(-(phi1_eme + phi2_eme * carga_poli[i])));
      }
    
      fruit_eme ~ normal(mu_eme, sigma_eme);
    
      // model primadonna
      vector[N] mu_pri;
      alpha_pri ~ normal(17, 2);
      phi1_pri ~ normal(0, 1);
      phi2_pri ~ normal(0, 1);
      sigma_pri ~ exponential(1);
      
      for (i in 1:N) { 
        mu_pri[i] = alpha_pri/(1+exp(-(phi1_pri + phi2_pri * carga_poli[i])));
      }
    
      fruit_pri ~ normal(mu_pri, sigma_pri);
    
    
      // model san joaquin
      vector[N] mu_sj;
      alpha_sj ~ normal(17, 2);
      phi1_sj ~ normal(0, 1);
      phi2_sj ~ normal(0, 1);
      sigma_sj ~ exponential(1);
      
      for (i in 1:N) { 
        mu_sj[i] = alpha_sj/(1+exp(-(phi1_sj + phi2_sj * carga_poli[i])));
      }
    
      fruit_sj ~ normal(mu_sj, sigma_sj);

    }
    
    generated quantities{
      array[N] real ppcheck_sch;
      vector[N] mu_sch1;
      array[N] real ppcheck_pri;
      vector[N] mu_pri1;
      array[N] real ppcheck_eme;
      vector[N] mu_eme1;
      array[N] real ppcheck_sj;
      vector[N] mu_sj1;
      
      for (i in 1:N) { 
        mu_sch1[i] = alpha_sch/(1+exp(-(phi1_sch + phi2_sch * carga_poli[i])));
      }
    
      ppcheck_sch = normal_rng(mu_sch1, sigma_sch);
    
      for (i in 1:N) { 
        mu_eme1[i] = alpha_eme/(1+exp(-(phi1_eme + phi2_eme * carga_poli[i])));
      }
    
      ppcheck_eme = normal_rng(mu_eme1, sigma_eme);
    
      for (i in 1:N) { 
        mu_pri1[i] = alpha_pri/(1+exp(-(phi1_pri + phi2_pri * carga_poli[i])));
      }
    
      ppcheck_pri = normal_rng(mu_pri1, sigma_pri);
    
      for (i in 1:N) { 
        mu_sj1[i] = alpha_sj/(1+exp(-(phi1_sj + phi2_sj * carga_poli[i])));
      }
    
      ppcheck_sj = normal_rng(mu_sj1, sigma_sj);
    }
    ')

file <- paste(getwd(), '/pollen_function.stan', sep = '')
fit_fun_pollen <- cmdstan_model(file, compile = T)

mod_fun_pollen <- 
  fit_fun_pollen$sample(
    data = dat_pollen_fun, 
    iter_sampling = 10e3,
    iter_warmup = 500,
    chains = 3,
    parallel_chains = 3,
    thin = 3,
    refresh = 500,
    seed = 123
  )

out_mod_fun <- mod_fun_pollen$summary() 
out_mod_fun |> print(n = 20)

par(mfrow = c(3, 3))
for (i in 2:10) trace_plot(mod_fun_pollen, out_mod_fun$variable[i], 3)
par(mfrow = c(1, 1))

ppcheck_fun_poll <- mod_fun_pollen$draws(c('ppcheck_sch',
                                           'ppcheck_eme',
                                           'ppcheck_pri',
                                           'ppcheck_sj'), format = 'df')

ppcheck_fun_poll <- 
  list(sch = ppcheck_fun_poll[, grep('sch', colnames(ppcheck_fun_poll))],
       eme = ppcheck_fun_poll[, grep('eme', colnames(ppcheck_fun_poll))],
       sj = ppcheck_fun_poll[, grep('sj', colnames(ppcheck_fun_poll))],
       pri = ppcheck_fun_poll[, grep('pri', colnames(ppcheck_fun_poll))])

ppcheck_fun_poll <- lapply(ppcheck_fun_poll, as.matrix)

names(dat_pollen_fun)

par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))

for (i in 1:4) {
  plot(density(dat_pollen_fun[[i+1]]), col = 'red', main = '',
       xlab = paste('Fruit size', names(ppcheck_fun_poll)[i]), 
       lwd = 0.1)
  for (j in 1:100) lines(density(ppcheck_fun_poll[[i]][j, ]), lwd = 0.1)
  lines(density(dat_pollen_fun[[i+1]]), col = 'red', lwd = 1.5)
}

par(mfrow = c(1, 1))

post_functions <- mod_fun_pollen$draws(out_mod_fun$variable[2:17], 
                                       format = 'df')

post_functions <- 
  list(sch = post_functions[, grep('sch', colnames(post_functions))],
       eme = post_functions[, grep('eme', colnames(post_functions))],
       sj = post_functions[, grep('sj', colnames(post_functions))],
       pri = post_functions[, grep('pri', colnames(post_functions))])

z_x <- dat_pollen_fun %$% seq(min(carga_poli), max(carga_poli), 
                              length.out = 1e3)

par(mfrow = c(2, 2))
for (i in 1:4) {
  
  plot((mean(post_functions[[i]][, 1, drop = T])) /
         (1 + exp(-(mean(post_functions[[i]][, 2, drop = T]) +
                    mean(post_functions[[i]][, 3, drop = T]) * z_x))) ~ z_x, 
       lwd = 0.1, type = 'l', xlim = c(0, 400), 
       ylim = c(0, 20), xlab = 'Pollen deposition (z-scores)', 
       ylab = paste('Fruit diameter', names(post_functions)[i]), col = i)
  
  for (j in 1:100) {
    post_functions$sch %$% 
      lines(post_functions[[i]][, 1, drop = T][[j]] /
              (1 + exp(-(post_functions[[i]][, 2, drop = T][[j]] +
                         post_functions[[i]][, 3, drop = T][[j]] * z_x))) ~ z_x, 
            lwd = 0.1, type = 'l', col = i)
  }
  
  points(dat_pollen_fun$carga_poli, dat_pollen_fun[[i+1]], lwd = 2)
  
}

par(mfrow = c(1, 1))

plot_functions <- 
  lapply(post_functions, FUN = 
           function(x) {
             
             a <- x[, grep('alpha', colnames(x)), drop = T]
             phi1 <- x[, grep('phi1', colnames(x)), drop = T]
             phi2 <- x[, grep('phi2', colnames(x)), drop = T]
             s <- x[, grep('sigma', colnames(x)), drop = T]
             
             est <- sapply(1:length(z_x), FUN = 
                             function(i) {
                               
                               mu <- a / (1 + exp(-(phi1 + phi2 * z_x[i]))) 
                               mu
                               
                             }, simplify = 'array')
             
             sim <- sapply(1:length(z_x), FUN = 
                             function(i) {
                               
                               mu <- a / (1 + exp(-(phi1 + phi2 * z_x[i])))
                               rnorm(1e3, mu, s)
                               
                             }, simplify = 'array')
             
             sim <- apply(sim, 2, FUN = 
                            function(j) {
                              tibble(li = quantile(j, 0.025), 
                                     ls = quantile(j, 0.975))
                            }, simplify = 'list')
             
             sim <- do.call('rbind', sim)
             
             sim$mu <- apply(est, 2, mean)
             sim$x <- z_x
             
             
             as_tibble(sim)
             
           })

for (i in seq_along(plot_functions)) 
  plot_functions[[i]]$cultivar <- names(plot_functions)[i]

do.call('rbind', plot_functions) |> 
  ggplot(aes(x, mu, ymin = li, ymax = ls)) +
  geom_ribbon(aes(fill = cultivar), alpha = 0.25) +
  geom_line(aes(color = cultivar)) +
  labs(x = 'Stigmatic pollen load (z-scores)', 
       y = 'Predicted fruit diameter (mm)') +
  theme_bw() +
  theme(panel.grid = element_blank())

predict_fruit_size <- function(x, 
                               cultivar = 'sch', 
                               mean_est = TRUE,
                               seed = 123) {
  if (length(x) <= 1) {
    stop(message('x length must be greater than 1'))
  }
  df <- post_functions[[grep(cultivar, names(post_functions))]]
  a <- df[, 1, drop = T]
  phi1 <- df[, 2, drop = T]
  phi2 <- df[, 3, drop = T]
  sigma <- df[, 4, drop = T]
  zeros <- which(x == 0)
  non_zeros <- which(x > 0)
  x_temp <- x[non_zeros]
  
  n <- length(x_temp)
  
  if (mean_est) {
    set.seed(seed)
    a <- sample(a, n, T)
    set.seed(seed)
    phi1 <- sample(phi1, n, T)
    set.seed(seed)
    phi2 <- sample(phi2, n, T)
    mu <- a / (1 + exp(-(phi1 + phi2 * x_temp)))
    
    out <- mu
    
  } else {
    set.seed(seed)
    a <- sample(a, n, T)
    set.seed(seed)
    phi1 <- sample(phi1, n, T)
    set.seed(seed)
    phi2 <- sample(phi2, n, T)
    set.seed(seed)
    sigma <- sample(sigma, n, T)
    
    mu <- a / (1 + exp(-(phi1 + phi2 * x_temp)))
    
    set.seed(seed)
    out <- rnorm(n, mu, sigma)
  }
  x[zeros] <- 0
  x[non_zeros] <- out
  return(x)
}

predict_fruit_size <- cmpfun(predict_fruit_size)


plot(0:500, predict_fruit_size(0:500, cultivar = 'sch', mean_est = F), col = 2)
points(0:500, predict_fruit_size(0:500, cultivar = 'sch', mean_est = T), 
       pch = 16)

# ======= corr size ~ weight ====

fruit_size <- readRDS('fruit_size.rds')

fruit_size <- fruit_size[fruit_size$fruit_diameter >5, ]

fruit_size <- na.omit(fruit_size[, -ncol(fruit_size)])

fruit_size$plant_id <- fruit_size %$% paste(farm, plant, sep = '_')
fruit_size$farm <- as.factor(fruit_size$farm)
fruit_size$plant_id <- as.factor(fruit_size$plant_id)
fruit_size$year <- as.factor(fruit_size$year)

summary(fruit_size)

dat_size_weight <- 
  lapply(fruit_size[, c("year", "farm", "fruit_weight", 
                        "fruit_diameter", "plant_id")], 
         function(x) if(is.factor(x)) as.integer(x) else(x))

dat_size_weight$N <- length(dat_size_weight$year)
dat_size_weight$N_year <- max(dat_size_weight$year)
dat_size_weight$N_farm <- max(dat_size_weight$farm)
dat_size_weight$N_plant <- max(dat_size_weight$plant_id)

cat(file = 'diameter_weight.stan', 
    '
    data{
      int N;
      int N_year;
      int N_farm;
      int N_plant;
      array[N] int year;
      array[N] int farm;
      array[N] int plant_id;
      vector[N] fruit_weight;
      vector[N] fruit_diameter;
    }
    
    parameters{
      vector[N_plant] z_alpha;
      real mu_alpha;
      real<lower = 0> sigma_alpha;
    
      vector[N_farm] z_theta;
      real mu_theta;
      real<lower = 0> sigma_theta;
    
      vector[N_year] z_tau;
      real mu_tau;
      real<lower = 0> sigma_tau;
    
      real beta;
      real<lower = 0> sigma;
    }
    
    transformed parameters{
      vector[N_plant] alpha;
      vector[N_farm] theta;
      vector[N_year] tau;
      
      alpha = mu_alpha + z_alpha * sigma_alpha;
      theta = mu_theta + z_theta * sigma_theta;
      tau = mu_tau + z_tau * sigma_tau;
    }
    
    model{
      vector[N] mu;
      beta ~ lognormal(0, 1);
      sigma ~ exponential(1);
      
      z_alpha ~ normal(0, 1);
      mu_alpha ~ normal(1.5, 0.25);
      sigma_alpha ~ exponential(1);
    
      z_theta ~ normal(0, 1);
      mu_theta ~ normal(0, 1);
      sigma_theta ~ exponential(1);
    
      z_tau ~ normal(0, 1);
      mu_tau ~ normal(0, 1);
      sigma_tau ~ exponential(1);
    
      for (i in 1:N) {
        mu[i] = alpha[plant_id[i]] + theta[farm[i]] +
                tau[year[i]] + exp(beta*fruit_diameter[i]);
      }
    
      fruit_weight ~ normal(mu, sigma);
    }
    
    generated quantities{
      vector[N] log_lik;
      vector[N] mu1;
      array[N] real ppcheck;
    
      for (i in 1:N) {
        mu1[i] = alpha[plant_id[i]] + theta[farm[i]] +
                 tau[year[i]] + exp(beta*fruit_diameter[i]);
      }
    
      for (i in 1:N) log_lik[i] = normal_lpdf(fruit_weight[i] | mu1[i], sigma);
      
      ppcheck = normal_rng(mu1, sigma);
    }
    ')

file <- paste(getwd(), '/diameter_weight.stan', sep = '')

fit_size_weight <- cmdstan_model(file, compile = T)

mod_size_weight <- 
  fit_size_weight$sample(
    data = dat_size_weight, 
    iter_sampling = 3e3,
    iter_warmup = 500,
    chains = 3,
    parallel_chains = 3,
    thin = 3, 
    refresh = 1e3,
    seed = 123
  )


out_mod_size_weight <- mod_size_weight$summary(variables = 
                                                 c('alpha', 'theta', 'tau', 'beta', 
                                                   'sigma'))

mod_diagnostics(mod_size_weight, out_mod_size_weight)

ppcheck_size_weight <- mod_size_weight$draws('ppcheck', format = 'matrix')

plot(density(ppcheck_size_weight[1, ]), lwd = 0.1, 
     xlab = 'Fruit weight (g)', main = '',  ylim = c(0, 1))
for (i in 1:500) lines(density(ppcheck_size_weight[i, ]), 
                       lwd = 0.1)
lines(density(dat_size_weight$fruit_weight), 
      col = 'red', lwd = 2)

post_size_weight <- mod_size_weight$draws(c('alpha', 'theta', 'tau', 'beta', 
                                            'sigma'),
                                          format = 'df')

post_size_weight <- 
  list(plant = post_size_weight[, grep('alpha', colnames(post_size_weight))],
       farm = post_size_weight[, grep('theta', colnames(post_size_weight))],
       year = post_size_weight[, grep('tau', colnames(post_size_weight))],
       beta = post_size_weight[, grep('beta', colnames(post_size_weight))],
       sigma = post_size_weight[, grep('sigma', colnames(post_size_weight))])

predict_weight <- 
  function(x, mu_est = T, seed = 123) {
    
    zeros <- which(x == 0)
    non_zeros <- which(x > 0)
    
    n <- length(non_zeros)
    
    x_tem <- x[non_zeros]
    
    df <- lapply(post_size_weight[-c(4:5)], FUN = 
                   function(i) {
                     x <- rowMeans(i)
                     x <- rep(x, 20)
                   })
    
    df$beta <- rep(post_size_weight$beta$beta, 20)
    df$sigma <- rep(post_size_weight$sigma$sigma, 20)
    
    if (mu_est) {
      df <- lapply(df, FUN =
                     function(vec){
                       vec[1:n]
                     })
      sigma <- df[[5]]
      
      mu <- df[[1]] + df[[2]] + df[[3]] + exp(df[[4]]*x_tem)
      out <- mu
      
      if (sum(out < 0) >= 1) {
        n2 <- sum(out < 0)
        set.seed(seed)
        out[which(out < 0 )] <- rgamma(n2, 0.03/0.01, 1/0.01)
      }
      
    } else {
      
      df <- lapply(df, FUN =
                     function(vec){
                       vec[1:n]
                     })
      sigma <- df[[5]]
      
      mu <- df[[1]] + df[[2]] + df[[3]] + exp(df[[4]]*x_tem)
      
      set.seed(seed)
      out <- rnorm(n, mu, sigma)
      if (sum(out < 0) >= 1) {
        n2 <- sum(out < 0)
        set.seed(seed)
        out[which(out < 0 )] <- rgamma(n2, 0.03/0.01, 1/0.01)
      }
    }
    
    x[zeros] <- 0
    x[non_zeros] <- out
    return(x)
  }
predict_weight <- cmpfun(predict_weight)


t <- crop_pollination(p_ha = p_01ha,
                      flowers_plant = total_flowers, 
                      visits_bee = visits_day, 
                      bees_hive = hives_ha(10, seed = i+500), 
                      hive_aggregate = T, 
                      short = F)

size <- predict_fruit_size(t$Hive10$plant1, mean_est = T) 
weight <- predict_weight(size, mu_est = T)


par(mfrow = c(1, 3), mar = c(4, 4, 1, 1))
plot(density(t$Hive10$plant1))
plot(density(predict_fruit_size(t$Hive10$plant1, mean_est = F)), 
     col = 'red', lwd = 2, xlab = 'Fruit diameter mm', main = '')
lines(density(size), lwd = 2)
plot(density(predict_weight(size, mu_est = F)), 
     xlab = 'Fruit weight g', main = '', lwd = 2, col = 'red')
lines(density(weight), lwd = 2)
par(mfrow = c(1, 1))

xx <- sample(0:20, 10e3, T)

par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
plot(xx, predict_weight(xx, mu_est = F), 
     xlab = 'Fruit diameter (mm)', 
     ylab = 'Predicted weight (g)', 
     main = expression('Incluring'~sigma~'parameter'))
plot(xx, predict_weight(xx, mu_est = T), 
     xlab = 'Fruit diameter (mm)', 
     ylab = 'Predicted weight (g)', 
     main = 'average estimate')
par(mfrow = c(1, 1))

dat_size_weight %$% plot(fruit_diameter, fruit_weight)
lines(0:25, predict_weight(0:25, mu_est = T), col = 'red')

# ====== 8. Crop production =======

crop_yield <- function(p_ha, # plants per ha
                       flowers_plant, # total flowers per plant
                       beta1 = 6, #flowering percentage (par 1 beta distribution) 
                       beta2 = 4, #flowering percentage (par 2 beta distribution) 
                       visits_bee, # visits per honeybee individual
                       bees_hive, # number of hives and honeybees per hive
                       hive_aggregate = T, # if T same plants are visited by all hives 
                       # if F each hive has its own plants
                       # flowers
                       seed = 123, 
                       fruit_diameter = F,
                       average_diameter = T,
                       average_weight = T, 
                       cultivar = 'sch') {
  
  message('Starting floral visits')
  t1 <- Sys.time()
  
  pollen_deposition <- 
    crop_pollination(p_ha,
                     flowers_plant,
                     beta1,
                     beta2,
                     visits_bee,
                     bees_hive,
                     hive_aggregate,
                     short = F, 
                     seed)
  
  message('Starting production simulation (t per ha)')
  
  production_plant <-
    lapply(pollen_deposition, FUN =
               function(colmena) {
                 
                 mu <- lapply(colmena, FUN =
                                function(planta) {
                                  t <- predict_fruit_size(planta, cultivar = cultivar,
                                                          mean_est = average_diameter)
                                  t1 <- predict_weight(t, mu_est = average_weight)
                                  
                                  tibble(kg_plant = sum(t1)/1e3,
                                         mu_fruit_size = mean(t),
                                         sd_fruit_size = sd(t),
                                         mu_fruit_weight = mean(t1),
                                         sd_fruit_weight = sd(t1))
                                })
                 
                 mu <- do.call('rbind', mu)
                 mu$plant <- paste('plant', 1:length(colmena))
                 mu
               })
  
  for (i in seq_along(production_plant)) production_plant[[i]]$hives <- i
  
  production_ha <-
    lapply(production_plant, FUN =
               function(ha) {
                 tibble(t = sum(ha$kg_plant)/1e3)
               })
  for (i in seq_along(production_ha)) production_ha[[i]]$hives <- i
  
  production_ha <- do.call('rbind', production_ha)
  
  output <-
    list(production_ha = production_ha,
         production_ha_plant = production_plant)
  
  return(output)
}


cluster <- makeCluster(8)

clusterExport(cluster, c('simulated_visits', 'pollen_deposition_fun',
                         'total_flowers', 'visits_day', 'hives_ha',
                         'p_01ha', 'visits_day_HQ', 'crop_pollination',
                         'asymptote', 'slope', 'crop_yield',
                         'predict_fruit_size', 'predict_weight',
                         'post_functions', 'post_size_weight'))

clusterEvalQ(cluster, {
  pks <- c('magrittr', 'parallel', "compiler", "tidyverse")

  sapply(pks, library, character.only = T)
})


t <- Sys.time()
t_ha_LQ <- parLapply(cluster, 1:50, fun = 
                    function(i) {
                      message(paste('simulation', i))
                      
                      l <- crop_yield(p_ha = p_01ha,
                                      flowers_plant = total_flowers, 
                                      visits_bee = visits_day, 
                                      bees_hive = hives_ha(20, 
                                                           mu_pop = 1e4,
                                                           seed = 500+i), 
                                      hive_aggregate = T,
                                      average_diameter = F, 
                                      average_weight = F, 
                                      cultivar = 'sch')
                      
                      l$production_ha$sim <- paste('sim', i, sep = '')
                      l$production_ha$type <- 'LQ'
                      l$production_ha_plant <- 
                        do.call('rbind', l$production_ha_plant)
                      l$production_ha_plant$sim <- paste('sim', i, sep = '')
                      l$production_ha_plant$sim <- paste('sim', i, sep = '')
                      l$production_ha_plant$type <- 'LQ'
                      l
                    })

Sys.time() - t

names(t_ha_LQ) <- paste('sim', 1:length(t_ha_LQ), sep = '')

save.image('all_results.RData')

t <- Sys.time()
t_ha_HQ <- parLapply(cluster, 1:50, fun = 
                    function(i) {
                      message(paste('simulation', i))
                      
                      l <- crop_yield(p_ha = p_01ha,
                                      flowers_plant = total_flowers, 
                                      visits_bee = visits_day_HQ, 
                                      bees_hive = hives_ha(20, 
                                                           mu_pop = 2e4,
                                                           seed = 500+i), 
                                      hive_aggregate = T,
                                      average_diameter = F, 
                                      average_weight = F, 
                                      cultivar = 'sch')
                      
                      l$production_ha$sim <- paste('sim', i, sep = '')
                      l$production_ha$type <- 'HQ'
                      l$production_ha_plant <- 
                        do.call('rbind', l$production_ha_plant)
                      l$production_ha_plant$sim <- paste('sim', i, sep = '')
                      l$production_ha_plant$type <- 'HQ'
                      l
                    })

Sys.time() - t

stopCluster(cluster)
rm(list = "cluster")

names(t_ha_HQ) <- paste('sim', 1:length(t_ha_HQ), sep = '')

save.image('all_results.RData')

sims_lq <- lapply(t_ha_LQ, FUN = 
                 function(x) {
                   x$production_ha
                 })

sims_lq <- do.call('rbind', sims_lq)

sims_hq <- lapply(t_ha_HQ, FUN = 
                    function(x) {
                      x$production_ha
                    })

sims_hq <- do.call('rbind', sims_hq)


sims <- rbind(sims_lq, sims_hq)


sims |> 
  ggplot(aes(hives, t, linetype = sim, color = type)) +
  geom_line() +
  scale_color_manual(values = c('tan1', 'lightblue')) +
  scale_linetype_manual(values = rep(1, 60)) +
  theme_bw() +
  labs(y = 't ha', x = 'Hive density (ha)') +
  theme(legend.position = 'none')

sims |> 
  ggplot(aes(as.numeric(hives), t, linetype = sim, color = type)) +
  geom_jitter(size = 0.5) +
  scale_color_manual(values = c('tan1', 'lightblue')) +
  scale_shape_manual(values = rep(1, 60)) +
  theme_bw() +
  labs(y = 't ha', x = 'Hive density (ha)') +
  theme(legend.position = 'none')

sims_lq2 <- lapply(t_ha_LQ, FUN = 
                    function(x) {
                      x$production_ha_plant
                    })

sims_lq2 <- do.call('rbind', sims_lq2)

sims_lq2$type <- 'LQ'

sims_hq2 <- lapply(t_ha_HQ, FUN = 
                    function(x) {
                      x$production_ha_plant
                    })

sims_hq2 <- do.call('rbind', sims_hq2)


sims2 <- rbind(sims_lq2, sims_hq2)

sims2$hives <- as.factor(sims2$hives)

sims2 |> 
  ggplot(aes(hives, kg_plant, fill = type)) +
  geom_boxplot(outlier.alpha = 0) +
  lims(y = c(0, 15)) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  theme_bw() +
  theme(legend.position = 'none')
  

sims2.2 <- 
  sims2 |> 
  group_by(hives, type) |> 
  transmute(mu_t = median(kg_plant), 
            li_t = quantile(kg_plant, 0.025), 
            ls_t = quantile(kg_plant, 0.975), 
            mu_sd_f = median(sd_fruit_size), 
            li_sd_f = quantile(sd_fruit_size, 0.025), 
            ls_sd_f = quantile(sd_fruit_size, 0.975)) |> 
  unique()

sims2.2 |> 
  ggplot(aes(as.numeric(hives), y = mu_t, 
             ymin = li_t, ymax = ls_t, color = type, 
             fill = type)) +
  geom_ribbon(alpha = 0.3) +
  geom_line()

sims2 |> 
  ggplot(aes(hives, mu_fruit_size, fill = type)) +
  geom_boxplot(outlier.alpha = 0) +
  lims(y = c(0, 15)) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  theme_bw() +
  theme(legend.position = 'none')


sims2 |> 
  ggplot(aes(hives, sd_fruit_size, fill = type)) +
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  theme_bw() +
  theme(legend.position = 'none')

sims2 |> 
  ggplot(aes(hives, sd_fruit_weight, fill = type)) +
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  theme_bw() +
  theme(legend.position = 'none')



summary(sims$t)







t <- Sys.time()
t_ha_LQ <- crop_yield(p_ha = p_01ha,
                      flowers_plant = total_flowers, 
                      visits_bee = visits_day, 
                      bees_hive = hives_ha(1,
                                           mu_pop = 1e4,
                                           seed = 500), 
                      hive_aggregate = T,
                      average_diameter = F, 
                      average_weight = F)
Sys.time() - t

t <- Sys.time()
t_ha_HQ <- crop_yield(p_ha = p_01ha,
                      flowers_plant = total_flowers, 
                      visits_bee = visits_day_HQ, 
                      bees_hive = hives_ha(5, 
                                           mu_pop = 2e4,
                                           seed = 500), 
                      hive_aggregate = T,
                      average_diameter = F, 
                      average_weight = F)
Sys.time() - t

t_ha_LQ$production_ha %$% plot(hives, t, type = 'l', col = 'red', 
                               xlim = c(0, 10), ylim = c(0, 20))
t_ha_HQ$production_ha %$% lines(hives, t, type = 'l', col = 'blue')



####### continuar con las figuras de los outputs





pollen_deposition_LQ <- 
  lapply(pollen_deposition_LQ, FUN = 
           function(x) {
             lapply(x, FUN = 
                      function(j) {
                        
                        fruto <- sapply(j, FUN = 
                                          function(k) {
                                            if(k == 0) 0
                                            else production_function(k)
                                          })
                        sum(fruto)/1e3
                      })
           })

pollen_deposition_HQ <- 
  lapply(pollen_deposition_HQ, FUN = 
           function(x) {
             lapply(x, FUN = 
                      function(j) {
                        
                        fruto <- sapply(j, FUN = 
                                          function(k) {
                                            if(k == 0) 0
                                            else production_function(k)
                                          })
                        sum(fruto)/1e3
                      })
           })


pollen_deposition_HQ <- 
  lapply(1:length(pollen_deposition_HQ), FUN = 
           function(x) {
             
             t <- unlist(pollen_deposition_HQ[x], use.names = F)
             
             tibble(kg = t, 
                    hive = paste(x), 
                    quality = 'Hight')
             
           })

pollen_deposition_LQ <- 
  lapply(1:length(pollen_deposition_LQ), FUN = 
           function(x) {
             
             t <- unlist(pollen_deposition_LQ[x], use.names = F)
             
             tibble(kg = t, 
                    hive = paste(x), 
                    quality = 'Hight')
             
           })

pollen_deposition_LQ <- do.call('rbind', pollen_deposition_LQ)
pollen_deposition_LQ$quality <- 'Low'
pollen_deposition_HQ <- do.call('rbind', pollen_deposition_HQ)

production <- rbind(pollen_deposition_HQ, 
                    pollen_deposition_LQ)

production$hive <- as.factor(production$hive)
production$hive <- factor(production$hive, 
                          levels = as.character(sort(as.numeric(levels(production$hive)))))

production <- split(production, list(production$quality, 
                                     production$hive))

production$Hight.1

production <- lapply(production, FUN = 
                       function(x) {
                         
                         v <- vector('double', 2e3)
                         
                         for (i in 1:2e3) {
                           v[[i]] <- sum(sample(x$kg, p_01ha, replace = T))
                         }
                         
                         v <- v/1e3
                         
                         tibble(mu = mean(v), 
                                li = quantile(v, 0.025), 
                                ls = quantile(v, 0.975), 
                                hive = x$hive[[1]],
                                quality = x$quality[[1]])
                       })

production <- do.call('rbind', production)

production |> 
  ggplot() +
  geom_line(data = production[production$quality == 'Low', ], 
            aes(hive, mu, group = 1), color = 'tan1') +
  geom_vline(xintercept = c(4, 9, 10, 13), linetype = 2, col = 'red') +
  geom_point(data = production[production$quality == 'Low', ], 
             aes(hive, mu), color = 'tan1') +
  geom_errorbar(data = production[production$quality == 'Low', ], 
                aes(hive, mu, ymin = li, ymax = ls), width = 0, 
                color = 'tan1') +
  geom_line(data = production[production$quality == 'Hight', ], 
            aes(hive, mu, group = 1), color = 'lightblue3') +
  geom_point(data = production[production$quality == 'Hight', ], 
             aes(hive, mu), color = 'lightblue3') +
  geom_errorbar(data = production[production$quality == 'Hight', ], 
                aes(hive, mu, ymin = li, ymax = ls), width = 0, 
                color = 'lightblue3') +
  labs(x = 'Hive density per blueberry ha', y = 'Average production per ha (t)') +
  theme_bw() +
  theme(legend.position = 'top', 
        panel.grid = element_blank())

ggsave('production_preliminary.jpg', width = 14, height = 8, 
       units = 'cm', dpi = 700)

quantile(unlist(pollen_deposition_HQ$Hive5, use.names = F))
















optimal_model <- bf(fruto_diam ~ carga_poli + carga_poli2 + tratamiento)

get_prior(optimal_model, data = pollen_fruit, family = "Gaussian")

prior_optMod <- c(set_prior('normal(5, 3)', class = 'Intercept'),
                  set_prior('normal(0, 0.2)', class = 'b', coef = 'carga_poli'),
                  set_prior('normal(0, 0.2)', class = 'b', coef = 'carga_poli2'),
                  set_prior('exponential(1)', class = 'sigma'))

optim_model <- brm(formula = optimal_model, data = pollen_fruit, family = gaussian,
                   prior = prior_optMod, future = T, chains = 3,
                   thin = 3, iter = 6e3, warmup =  500, seed = 222)

summary(optim_model)

coefs_optim <- fixef(optim_model)

a <- coefs_optim[3, 1]
b <- coefs_optim[2, 1]

(x_value <- -b/(2*a)) # function vertex
(y_value <- (a*(x_value)^2) + (b*x_value) + coefs_optim[1, 1]) # asíntot

# =====================================


x <- seq(0,120, 1)
#y <- 250/(1+18*exp(-1*x))
a <- 250 # slope 
b1 <- 0.5 # 
y <- a*(1-exp(-b1*x))
plot(x, y, type="l",main="Benefit and cost curve\n honeybee - blueberry interaction")
points(7, y[x==7], col = 'red', pch = 19)
a2 <- 150 # asymptote
b2 <- 0.1 # slope factor
g <- 60 # response to x a the middle point 
y1 <- a2/(1+exp(b2*(g-x)))
lines(x,y1,ylim=c(0,140),type="l",main="four-parameter logistic", col = 'red') #two-parameter logistic
lines(x, y - y1 ,ylim=c(0,140),type="l",main="four-parameter logistic", 
      col = 'red', lty = 2)


# ====== modelo sim =====

x <- seq(0,120, 1)
#y <- 250/(1+18*exp(-1*x))
a <- 200 # slope 
b1 <- 0.2 # 
y <- a*(1-exp(-b1*x))
plot(x, y, type="l", ylim = c(0, 250),
     main="Benefit and cost curve\n honeybee - blueberry interaction")
points(7, y[x==7], col = 'red', pch = 19)
a2 <- 100 # asymptote
b2 <- 0.1 # slope factor
g <- 60 # response to x a the middle point 
y1 <- a2/(1+exp(b2*(g-x)))
lines(x,y1,ylim=c(0,140),type="l",main="four-parameter logistic", col = 'red') #two-parameter logistic
lines(x, y - y1 ,ylim=c(0,140),type="l",main="four-parameter logistic", 
      col = 'red', lty = 2)

x

plot(rnbinom(length(x), size = 1, 
             mu = a2/(1+exp(b2*(g-x)))) ~ x)

pred <- rnbinom(length(x), size = 0.5, 
                mu = a2/(1+exp(b2*(g-x))))

cos_pois <- rpois(length(x), lambda = a2/(1+exp(b2*(g-x))))

benef_pois <- rpois(length(x), lambda = a*(1-exp(-b1*x)))

plot(benef_pois ~ x)
points(cos_pois ~ x, col = 'red')
points((benef_pois - cos_pois) ~ x, col = 'green')


cos_NB <- rnbinom(length(x), size = 0.5, 
                mu = a2/(1+exp(b2*(g-x))))

benef_NB <- rnbinom(length(x), mu = a*(1-exp(-b1*x)),
                    size = 0.5)

plot(benef_NB ~ x)
points(cos_NB ~ x, col = 'red')
points((benef_NB-cos_NB) ~ x, col = 'green')


# ========== Modelo tasa de visitas =======








