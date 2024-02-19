pks <- c('tidyverse', 'rethinking', 'rstan', 'magrittr', 'cmdstanr',
         'ggdag', 'dagitty', 'readxl')
sapply(pks, library, character.only = T)


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

# ========== Hive density and strensth ===========

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


par(mfrow = c(2, 2))

plot(density(trip_frequency_HQ), col = 'gold3', main = '', ylim = c(0, 5.5),
     xlab = 'Number of foraging trips (trips/90 min)', lwd = 2)
lines(density(trip_frequency_LQ), col = 'cyan4', main = '',
      xlab = 'Number of foraging trips (trips/90 min)', lwd = 2)

plot(density(trip_span_HQ/60), col = 'gold3', main = '', xlim = c(10, 50),
     xlab = 'Trip span (min)', lwd = 2)
lines(density(trip_span_LQ/60), main = '', col = 'cyan4',
      xlab = 'Trip span (min)', lwd = 2)

plot(density(trips_bee_HQ), main = '', col = 'gold3',
     xlab = 'Daily trip per bee (min)', lwd = 2, ylim = c(0, 1))
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
plot(density(rgamma(1e3, 10/2, 1/2)))

cat(file = 'honeybee_visit.stan', 
    '
    data{
      int N;
      vector[N] visit;
    }
    
    parameters{
      real<lower = 0> mu;
      real<lower = 0> sigma;
    }
    
    model{
      mu ~ normal(10, 2.5);
      sigma ~ exponential(1);
      
      visit ~ gamma(mu/sigma, 1/sigma);
    }
    ')

file <- paste(getwd(), '/honeybee_visit.stan', sep = '')
fit_visit_honeybee <- cmdstan_model(file, compile = T)

mod_visit_honeybee <- 
  fit_visit_honeybee$sample(
    data = list(visit = visit_span, 
                N = length(visit_span)),
    iter_sampling = 2e3, 
    iter_warmup = 500, 
    chains = 3, 
    parallel_chains = 3, 
    thin = 3, 
    seed = 123, 
    refresh = 500
  )

mod_visit_honeybee$summary()

output_mod_visit_honeybee <- mod_visit_honeybee$summary()

trace_plot <- function (fit, par, n_chains) {
  
  d <- as_mcmc.list(fit)
  
  d <- lapply(d, FUN = 
                function(x) {
                  x <- as.data.frame(unclass(as.matrix(x)))
                  x$iter <- 1:nrow(x)
                  x
                })
  
  for (i in 1:n_chains) d[[i]]$chain <- i
  
  plot(d[[1]][, 'iter', drop = T], d[[1]][, par, drop = T], 
       type = 'l', col = 1, main = '', ylab = par, 
       xlab = 'Iteration')
  for (i in 2:n_chains) {
    lines(d[[i]][, 'iter', drop = T], d[[i]][, par, drop = T], 
          col = i, main = '', ylab = par, 
          xlab = 'Iteration')
  }
  
} 

post_visit_honeybee <- mod_visit_honeybee$draws(format = 'df')

par(mfrow = c(1, 2))
for (i in 2:3) trace_plot(mod_visit_honeybee, output_mod_visit_honeybee$variable[i], 3)
par(mfrow = c(1, 1))

ppcheck_visit_honeybee <- 
  sapply(1:length(visit_span), FUN = 
           function(x) {
             
             rgamma(1e3, post_visit_honeybee$mu[x]/post_visit_honeybee$sigma[x],
                    1/post_visit_honeybee$sigma[x])
             
           })

plot(density(ppcheck_visit_honeybee[1, ]), ylim = c(0, 0.09), lwd = 0.1, 
     main = '', xlab = 'Honeybee time of floral visit (s)')
for (i in 1:100) lines(density(ppcheck_visit_honeybee[i, ]), lwd = 0.1)
lines(density(visit_span), col = 'red', lwd = 2)


time_visit_honeybee <- as.vector(ppcheck_visit_honeybee)

lines(density(time_visit_honeybee), col = 'tan1', lwd = 3)

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


plot(density(visits_day), main = '', xlab = 'Number of floral visits of\n a single forager honeybee bee per day')

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
    '
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
      vector[N_year] year;
      vector[N_site] locality;
      vector[N_farm] farm;
      vector[N_plant] plant;
      //vector[N_branch] branch;
      vector[N_treatment] t;
    }
    
    model{
      vector[N] p;
      year ~ normal(0, 1);
      locality ~ normal(0, 1);
      farm ~ normal(0, 1);
      plant ~ normal(0, 1);
      //branch ~ normal(0, 1);
      t ~ normal(0, 1);
    
      for (i in 1:N) {
        p[i] = t[treatment[i]] + 
                locality[locality_id[i]] + 
                farm[farm_id[i]] +
                plant[plant_id[i]] +
                //branch[branch_id[i]] +
                year[year_id[i]]; 
        
        p[i] = inv_logit(p[i]);
      }
      
      fruits ~ binomial(flowers, p);
    }
    ')

file <- paste(getwd(), '/fruit_set.stan', sep = '')

fit_fruit_set <- cmdstan_model(file, compile = T)

mod_fs_exp <- 
  fit_fruit_set$sample(
    data = dat_fs, 
    chains = 3, 
    parallel_chains = 3, 
    iter_warmup = 500, 
    iter_sampling = 2e3,
    thin = 3,
    seed = 123,
    refresh = 200
  )

output_mod_fruitset <- mod_fs_exp$summary() 

output_mod_fruitset |> print(n = 464)

trace_plot(mod_fs_exp, output_mod_fruitset$variable[1], 3)

par(mfrow = c(2, 4))
for (i in 2:9) {
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
       year = post_fruitset[, grep('^year', colnames(post_fruitset))],
       BH = post_fruitset[, grep('^BH', colnames(post_fruitset))])

ppcheck_fruit_set <- 
  sapply(1:length(dat_fs$plant_id), FUN = 
           function(x) {
             
             p <- with(dat_fs, 
                       {
                         t <- treatment[x]
                         l <- locality_id[x]
                         f <- farm_id[x]
                         p <- plant_id[x]
                         #b <- branch_id[x]
                         y <- year_id[x]
                         
                         
                         post_fruitset$treatment[, t, drop = T] +
                           post_fruitset$locality[, l, drop = T] +
                           post_fruitset$farm[, f, drop = T] +
                           post_fruitset$plant[, p, drop = T] +
                           #post_fruitset$branch[, b, drop = T] +
                           post_fruitset$year[, y, drop = T] 
                         
                       })
             
             rbinom(1e3, dat_fs$flowers[x], inv_logit(p))
             
           })

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
     xlab = 'blueberry fruit set')

mean(fruit_set)


# ======= 2. fruit produces  ====

fruit_plant <- readRDS('fruit_plant.rds')

fruit_plant$plot <- as.factor(paste(fruit_plant$farm, fruit_plant$plot, sep = '_'))

fruit_plant <- fruit_plant[, c("farm", "plot", "total_fruts")]

fruit_plant <- lapply(fruit_plant, function(x) if(is.factor(x)) as.numeric(x) else(x))

fruit_plant$N <- length(fruit_plant$farm)
fruit_plant$N_farm <- length(unique(fruit_plant$farm))
fruit_plant$N_plot <- length(unique(fruit_plant$plot))
fruit_plant$total_fruts <- round(fruit_plant$total_fruts)

plot(density(rnbinom(1e3, size = 2, mu = exp(8))))

cat(file = 'fruits_plant.stan', 
    '
    data{
      int N;
      int N_farm;
      int N_plot;
      array[N] int total_fruts;
      array[N] int farm;
      array[N] int plot;
    }
    
    parameters{
      vector[N_farm] alpha;
      vector[N_plot] tau;
      real<lower = 0> scale;
    }
    
    model{
      vector[N] mu;
      alpha ~ normal(7, 2);
      tau ~ normal(0, 0.5);
      scale ~ exponential(1);
    
      for (i in 1:N) {
        mu[i] = alpha[farm[i]] + tau[plot[i]];
        mu[i] = exp(mu[i]);
      }
    
      total_fruts ~ neg_binomial_2(mu, scale);  
    }
    ')


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
output_mod_tot_fru |> print(n = 21)

trace_plot(mod_fruit_plant, output_mod_tot_fru$variable[1], 3)

par(mfrow = c(1, 3))
for (i in 2:4) trace_plot(mod_fruit_plant, output_mod_tot_fru$variable[i], 3)
par(mfrow = c(1, 1))

post_tot_fruit <- mod_fruit_plant$draws(format = 'df')

post_tot_fruit <- 
  list(farm = post_tot_fruit[, grep('alpha', colnames(post_tot_fruit))], 
       plot = post_tot_fruit[, grep('tau', colnames(post_tot_fruit))], 
       scale = post_tot_fruit$scale)

ppcheck_tot_fru <- 
  sapply(1:length(fruit_plant$farm), FUN = 
           function(x) {
             
             mu <- with(fruit_plant, 
                        {
                          f <- farm[x]
                          p <- plot[x]
                          
                          post_tot_fruit$farm[, f, drop = T] +
                            post_tot_fruit$plot[, p, drop = T]
                        })
             
             rnbinom(1e3, size = post_tot_fruit$scale, mu = exp(mu))
             
           })

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
     xlab = 'Flowers per blueberry bush', xlim = c(0, 20e3)) 


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
           seed = 123) {
    
    set.seed(seed)
    plants <- sample(flowers_plant, p_ha, replace = T) # plants per ha 
    
    set.seed(seed)
    flowering_perc <- rbeta(p_ha, beta1, beta2) # flowering percentage parm
    
    plants <- round(plants * flowering_perc) # flowering percentage
    
    try <- plants
    
    plants <- lapply(plants, FUN = 
                       function(x) {
                         rep(0, x)
                       })
    
    
    names(plants) <- paste('plant', 1:p_ha, sep = '')
    
    foragers_bees <- bees_hive
    
    foragers_visits <- lapply(foragers_bees, FUN = 
                                function(x) {
                                  sample(visits_bee, x, replace = T)
                                })
    
    if (hive_aggregate) {
      
      plants <- 
        lapply(foragers_visits, FUN = 
                 function(colmena) {
                   
                   for (k in 1:2) { # number of pollination days (blueberry days receptivity)
                     message(paste('Starting day', k, 'of pollination'))
                     
                     for (j in seq_along(colmena)) { # bees on beehive 
                       
                       for (i in 1:colmena[j]){ # foraging trips per bee
                         
                         plant_i <- sample(1:length(plants), 1) # plant i
                         
                         flower_i <- sample(1:length(plants[[plant_i]]), 1) # flower j
                         
                         plants[[plant_i]][flower_i] <- plants[[plant_i]][flower_i] + 1 # visit to ij
                       }
                     }
                   }
                   
                   plants
                   
                 })
      
      i <- 1
      while (i < length(plants)) {
        
        for (j in 1:p_ha) plants[[i+1]][[j]] <- plants[[i]][[j]] + plants[[i+1]][[j]]
        
        i <- i + 1
      }
      
    } else {
      
      plants <- 
        lapply(foragers_visits, FUN = 
                 function(colmena) {
                   
                   for (k in 1:2) { # number of pollination days (blueberry days receptivity)
                     message(paste('Starting day', k, 'of pollination'))
                     
                     for (j in seq_along(colmena)) { # bees on beehive 
                       
                       for (i in 1:colmena[j]){ # foraging trips per bee
                         
                         plant_i <- sample(1:length(plants), 1) # plant i
                         
                         flower_i <- sample(1:length(plants[[plant_i]]), 1) # flower j
                         
                         plants[[plant_i]][flower_i] <- plants[[plant_i]][flower_i] + 1 # visit to ij
                       }
                     }
                   }
                   
                   plants
                   
                 })
      
    }
    
    names(plants) <- paste('Hive', 1:length(plants), sep = '')
    
    return(plants)
    
  }


t1 <- Sys.time()
p <- simulated_visits(p_ha = p_01ha,
                      flowers_plant = total_flowers, 
                      visits_bee = visits_day, 
                      bees_hive = rep(hives_ha(1), 20), 
                      hive_aggregate = T)
Sys.time() - t1

unlist(lapply(p$Hive1, sum), use.names = F)

sum(p$Hive1$plant1)

sum(p$Hive2$plant1) - sum(p$Hive1$plant1)

sum(p$Hive20$plant1) - sum(p$Hive19$plant1)

length(p)

par(mfrow = c(2, 2))
for (i in seq_along(p)) {
  
  plot(density(p[[i]][[1]]), ylim = c(0, 6), xlim = c(-0.5, 6),
       xlab = paste('Honeybee visits per flower\n (density:', 
                    i, ' hive per ha)'), 
       lwd = 0.1, main = '', col = i)
  
  for (j in 1:100) {
    lines(density(p[[i]][[j]]), col = i, lwd = 0.1) 
  }
}
par(mfrow = c(1, 1))





p$output$Hive1$plant5

p$frame

length(p$Hive1$plant1)
t <- 1

par(mfrow = c(2, 2))

repeat {
  
  if (t > 4) break
  
  plot(density(p[[t]][[t]]), ylim = c(0, 5), xlim = c(-0.5, 6), 
       main = paste('Hive', t), xlab = 'Honeybee visits per flower\n (density: 1 hive per ha)')
  for (i in seq_along(plants$Hive1)) {
    lines(density(p[[t]][[i]]), lwd = 0.01)
  }
  lines(density(unlist(lapply(plants[[t]], median))), col = 'red', lwd = 2)
  
  t <- t + 1
  
}
par(mfrow = c(1, 1))












# ======= second half of the simulation

# ====== Hoenybee pollen deposition ====
apisSVP <- readRDS('honeybee_svpd_data.rds')

apisSVP <- apisSVP[-which.max(apisSVP)] # removing 155 pollen deposition

tibble(x = c(0, 70)) |> 
  ggplot(aes(x)) +
  stat_function(fun = dlnorm, args = c(3, 1)) +
  labs(x = NULL, y = NULL, 
       title = expression(mu['honeybee pollen deposition'~'~'~'LogNormal(3, 1)']))

dat_apis <- list(visits = apisSVP)

plot(density(rnbinom(1e3, mu = exp(3), exp(1))))

model_apisNB <- 
  ulam(
    alist(
      visits ~ dgampois(lambda, phi),
      log(lambda) <- mu,
      mu ~ dnorm(3, 1),
      phi ~ dexp(1)
    ), data = dat_apis, chains = 3, cores = 3, warmup = 500, iter = 2e3
  )

precis(model_apisNB, depth = 2)
exp(precis(model_apisNB, depth = 2))

post_apis_NB <- extract.samples(model_apisNB)

pp_simNB <- sapply(1:100, simplify = 'array', FUN = 
                     function(i) {
                       with(post_apis_NB, 
                            {rnbinom(1e3, mu = exp(mu[i]), size = exp(phi[i]))})
                     })

model_apisPois <- 
  ulam(
    alist(
      visits ~ dpois(lambda),
      log(lambda) <- mu,
      mu ~ dnorm(3, 1)
    ), data = dat_apis, chains = 3, cores = 3, warmup = 500, iter = 2e3
  )

precis(model_apisPois, depth = 2)
exp(precis(model_apisPois, depth = 2))

post_apisPois <- extract.samples(model_apisPois)

pp_simPois <- sapply(1:100, FUN = 
                       function(i) {
                         with(post_apisPois, 
                              {rpois(1e3, exp(mu[i]))})
                       })
plot(density(rlnorm(1e3, log(17), log(2))))

model_apislog <- 
  ulam(
    alist(
      visits ~ dlnorm(mu, sigma),
      mu <- alpha,
      alpha ~ dnorm(17, 1),
      sigma ~ dexp(1)
    ), data = dat_apis, chains = 3, 
    cores = 3, warmup = 500, iter = 2e3
  )

precis(model_apislog, depth = 2)
exp(precis(model_apislog, depth = 2))

post_apislog <- extract.samples(model_apislog)

pp_simLOG <- sapply(1:100, FUN = 
                      function(i) {
                        with(post_apislog, 
                             {rlnorm(1e3, alpha[i], sigma[i])})
                      })

plot(density(apisSVP), col = 'tan1', lwd = 2, ylim = c(0, 0.08), main = '')
for (i in 1:100) lines(density(pp_simNB[, i]), lwd = 0.1)
for (i in 1:100) lines(density(pp_simPois[, i]), lwd = 0.1, col = 'red')
for (i in 1:100) lines(density(pp_simLOG[, i]), lwd = 0.1, col = 'green')
lines(density(apisSVP), col = 'tan1', lwd = 3)

# mejor el modelo gamma-possion (binomial negativa)



exp(mean(post_apis_NB$mu)) # Average pollen deposition
exp(quantile(post_apis_NB$mu, probs = c(0.025, 0.975))) # Credibility intervals

df <- tibble()
df2 <- tibble()


set.seed(1234)
for (j in 1:10) {
  
  sum_pol <- vector("double", 5000)
  for (i in 1:5000) {
    sum_pol[[i]] <- sum(sample(apisSVP, j, T))
  }
  
  df <- tibble(poll_ac = sum_pol,
               cero = rep(0, 5000),
               vis = as.factor(rep(paste("Visit", j), 5000)))
  
  df2 <- rbind(df2, df)
  
}

set.seed(1234)
try <- 
  apply(pp_sim[c(sample(1:ncol(pp_simNB), size = 50, replace = F)),], 1, 
        simplify = 'list', FUN = 
          function(x) {
            
            df_2 <- tibble()
            
            for (j in 1:10) {
              
              sum_vis <- vector("double", 2000)
              for (i in 1:2000) {
                sum_vis[[i]] <- sum(sample(x, size = j, T))
              }
              
              df_ <- tibble(poll_ac = sum_vis,
                            cero = rep(0, 2000),
                            vis = as.factor(rep(paste("Visit", j), 2000)))
              
              df_2 <- rbind(df_2, df_)
              
            }
            df_2
          })

names(try) <- paste('sim', 1:length(try))

for (i in seq_along(try)) {
  try[[i]]$sim <- rep(paste('sim', i), nrow(try[[i]]))
}

try <- do.call('rbind', try)
try$sim_vis <- paste(try$vis, try$sim)

df2$sim <- rep('obs', nrow(df2))
df2$sim_vis <- paste(df2$vis, df2$sim)


try <- rbind(try, df2)

colors <- grepl('obs', levels(as.factor(try$sim_vis)))
lines <- grepl('obs', unique(try$sim_vis))

try |> 
  ggplot(aes(poll_ac, y = after_stat(scaled))) +
  geom_density(aes(color = sim_vis), alpha = 0, linewidth = 0.05) +
  scale_color_manual(values = ifelse(colors == T, 'red', 'black')) +
  geom_density(data = try[try$sim == 'obs',], aes(poll_ac, y = after_stat(scaled)),
               alpha = 0, linewidth = 0.3, color = 'red') +
  geom_vline(xintercept = c(112, 274), linetype = 2) +
  geom_vline(xintercept = 192, linetype = 3) + 
  scale_x_continuous(limits = c(0, 500)) +
  facet_wrap(~ vis, scales = "free") +
  theme_bw() +
  labs(x = "SPL after sequential honeybee visits", y = NULL) +
  theme(legend.position = "none",
        panel.grid = element_blank())

try <- try[try$sim != 'obs', ]

try$sim <- as.factor(try$sim)

try <- split(try, list(try$vis, try$sim))

ci_prop_vis_ <- 
  lapply(try,
         function(x) {
           vec <- x$poll_ac
           vec2 <- sum(vec >= 112 & vec <= 274)/length(vec)
           
           tibble(prop = vec2)
         })


for (i in seq_along(ci_prop_vis_)) {
  ci_prop_vis_[[i]]$vis <- gsub('^(.*)\\.(.*)', "\\1", names(ci_prop_vis_)[i])
}

ci_prop_vis_ <- do.call('rbind', ci_prop_vis_)

ci_prop_vis_$vis <- as.factor(ci_prop_vis_$vis)

ci_prop_vis_ <- 
  lapply(split(ci_prop_vis_, ci_prop_vis_$vis), FUN = 
           function(x) {
             tibble(mean = mean(x$prop),
                    u_ic = quantile(x$prop, probs = 0.975),
                    l_ic = quantile(x$prop, probs = 0.025))
           })

for (i in seq_along(ci_prop_vis_)) ci_prop_vis_[[i]]$vis <- names(ci_prop_vis_)[i]

ci_prop_vis_ <- do.call('rbind', ci_prop_vis_)

ci_prop_vis_$vis <- as.numeric(gsub('^(.*)\\s(.)', '\\2', ci_prop_vis_$vis))

ci_prop_vis_ <- ci_prop_vis_[order(ci_prop_vis_$vis), ]

ci_prop_vis_$vis <- as.factor(ci_prop_vis_$vis)

ggplot(ci_prop_vis_) + 
  geom_errorbar(aes(ymin = l_ic, ymax = u_ic, x = vis),
                width = 0.2, color = 'tan1') + 
  geom_line(aes(vis, mean, group = 1), linetype = 2) +
  geom_point(aes(vis, mean), color = 'tan1', size = 2) +
  theme_bw() +
  labs(x = 'Number of honeybee visits',
       y = 'Prop. Optmimal SPLs') +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10.5))



# ====== Hoenybee pollen deposition ====
apisSVP <- readRDS('honeybee_svpd_data.rds')

apisSVP <- apisSVP[-which.max(apisSVP)]

tibble(x = c(0, 70)) |> 
  ggplot(aes(x)) +
  stat_function(fun = dlnorm, args = c(3, 1)) +
  labs(x = NULL, y = NULL, 
       title = expression(mu['honeybee pollen deposition'~'~'~'LogNormal(3, 1)']))

dat_apis <- list(visits = apisSVP)

plot(density(rnbinom(1e3, mu = exp(3), exp(1))))

model_apisNB <- 
  ulam(
    alist(
      visits ~ dgampois(lambda, phi),
      log(lambda) <- mu,
      mu ~ dnorm(3, 1),
      phi ~ dexp(1)
    ), data = dat_apis, chains = 3, cores = 3, warmup = 500, iter = 2e3
  )

precis(model_apisNB, depth = 2)
exp(precis(model_apisNB, depth = 2))

post_apis_NB <- extract.samples(model_apisNB)

pp_simNB <- sapply(1:100, simplify = 'array', FUN = 
                     function(i) {
                       with(post_apis_NB, 
                            {rnbinom(1e3, mu = exp(mu[i]), size = exp(phi[i]))})
                     })

model_apisPois <- 
  ulam(
    alist(
      visits ~ dpois(lambda),
      log(lambda) <- mu,
      mu ~ dnorm(3, 1)
    ), data = dat_apis, chains = 3, cores = 3, warmup = 500, iter = 2e3
  )

precis(model_apisPois, depth = 2)
exp(precis(model_apisPois, depth = 2))

post_apisPois <- extract.samples(model_apisPois)

pp_simPois <- sapply(1:100, FUN = 
                       function(i) {
                         with(post_apisPois, 
                              {rpois(1e3, exp(mu[i]))})
                       })
plot(density(rlnorm(1e3, log(17), log(2))))

model_apislog <- 
  ulam(
    alist(
      visits ~ dlnorm(mu, sigma),
      mu <- alpha,
      alpha ~ dnorm(17, 1),
      sigma ~ dexp(1)
    ), data = dat_apis, chains = 3, 
    cores = 3, warmup = 500, iter = 2e3
  )

precis(model_apislog, depth = 2)
exp(precis(model_apislog, depth = 2))

post_apislog <- extract.samples(model_apislog)

pp_simLOG <- sapply(1:100, FUN = 
                      function(i) {
                        with(post_apislog, 
                             {rlnorm(1e3, alpha[i], sigma[i])})
                      })

plot(density(apisSVP), col = 'tan1', lwd = 2, ylim = c(0, 0.08), main = '')
for (i in 1:100) lines(density(pp_simNB[, i]), lwd = 0.1)
for (i in 1:100) lines(density(pp_simPois[, i]), lwd = 0.1, col = 'red')
for (i in 1:100) lines(density(pp_simLOG[, i]), lwd = 0.1, col = 'green')
lines(density(apisSVP), col = 'tan1', lwd = 3)

# mejor el modelo gamma-possion (binomial negativa)



exp(mean(post_apis_NB$mu)) # Average pollen deposition
exp(quantile(post_apis_NB$mu, probs = c(0.025, 0.975))) # Credibility intervals

df <- tibble()
df2 <- tibble()


set.seed(1234)
for (j in 1:10) {
  
  sum_pol <- vector("double", 5000)
  for (i in 1:5000) {
    sum_pol[[i]] <- sum(sample(apisSVP, j, T))
  }
  
  df <- tibble(poll_ac = sum_pol,
               cero = rep(0, 5000),
               vis = as.factor(rep(paste("Visit", j), 5000)))
  
  df2 <- rbind(df2, df)
  
}

set.seed(1234)
try <- 
  apply(pp_sim[c(sample(1:ncol(pp_simNB), size = 50, replace = F)),], 1, 
        simplify = 'list', FUN = 
          function(x) {
            
            df_2 <- tibble()
            
            for (j in 1:10) {
              
              sum_vis <- vector("double", 2000)
              for (i in 1:2000) {
                sum_vis[[i]] <- sum(sample(x, size = j, T))
              }
              
              df_ <- tibble(poll_ac = sum_vis,
                            cero = rep(0, 2000),
                            vis = as.factor(rep(paste("Visit", j), 2000)))
              
              df_2 <- rbind(df_2, df_)
              
            }
            df_2
          })

names(try) <- paste('sim', 1:length(try))

for (i in seq_along(try)) {
  try[[i]]$sim <- rep(paste('sim', i), nrow(try[[i]]))
}

try <- do.call('rbind', try)
try$sim_vis <- paste(try$vis, try$sim)

df2$sim <- rep('obs', nrow(df2))
df2$sim_vis <- paste(df2$vis, df2$sim)


try <- rbind(try, df2)

colors <- grepl('obs', levels(as.factor(try$sim_vis)))
lines <- grepl('obs', unique(try$sim_vis))

try |> 
  ggplot(aes(poll_ac, y = after_stat(scaled))) +
  geom_density(aes(color = sim_vis), alpha = 0, linewidth = 0.05) +
  scale_color_manual(values = ifelse(colors == T, 'red', 'black')) +
  geom_density(data = try[try$sim == 'obs',], aes(poll_ac, y = after_stat(scaled)),
               alpha = 0, linewidth = 0.3, color = 'red') +
  geom_vline(xintercept = c(112, 274), linetype = 2) +
  geom_vline(xintercept = 192, linetype = 3) + 
  scale_x_continuous(limits = c(0, 500)) +
  facet_wrap(~ vis, scales = "free") +
  theme_bw() +
  labs(x = "SPL after sequential honeybee visits", y = NULL) +
  theme(legend.position = "none",
        panel.grid = element_blank())

try <- try[try$sim != 'obs', ]

try$sim <- as.factor(try$sim)

try <- split(try, list(try$vis, try$sim))

ci_prop_vis_ <- 
  lapply(try,
         function(x) {
           vec <- x$poll_ac
           vec2 <- sum(vec >= 112 & vec <= 274)/length(vec)
           
           tibble(prop = vec2)
         })


for (i in seq_along(ci_prop_vis_)) {
  ci_prop_vis_[[i]]$vis <- gsub('^(.*)\\.(.*)', "\\1", names(ci_prop_vis_)[i])
}

ci_prop_vis_ <- do.call('rbind', ci_prop_vis_)

ci_prop_vis_$vis <- as.factor(ci_prop_vis_$vis)

ci_prop_vis_ <- 
  lapply(split(ci_prop_vis_, ci_prop_vis_$vis), FUN = 
           function(x) {
             tibble(mean = mean(x$prop),
                    u_ic = quantile(x$prop, probs = 0.975),
                    l_ic = quantile(x$prop, probs = 0.025))
           })

for (i in seq_along(ci_prop_vis_)) ci_prop_vis_[[i]]$vis <- names(ci_prop_vis_)[i]

ci_prop_vis_ <- do.call('rbind', ci_prop_vis_)

ci_prop_vis_$vis <- as.numeric(gsub('^(.*)\\s(.)', '\\2', ci_prop_vis_$vis))

ci_prop_vis_ <- ci_prop_vis_[order(ci_prop_vis_$vis), ]

ci_prop_vis_$vis <- as.factor(ci_prop_vis_$vis)

ggplot(ci_prop_vis_) + 
  geom_errorbar(aes(ymin = l_ic, ymax = u_ic, x = vis),
                width = 0.2, color = 'tan1') + 
  geom_line(aes(vis, mean, group = 1), linetype = 2) +
  geom_point(aes(vis, mean), color = 'tan1', size = 2) +
  theme_bw() +
  labs(x = 'Number of honeybee visits',
       y = 'Prop. Optmimal SPLs') +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10.5))





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








