pks <- c('tidyverse', 'rethinking', 'rstan', 'magrittr', 'cmdstanr',
         'ggdag', 'dagitty', 'readxl', 'brms')
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


plot(density(visits_day), main = '', xlab = 'Number of floral visits of\n a single forager honeybee bee per day',
     lwd = 4, col = 'lightblue')

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
     xlab = 'blueberry fruit set', col = 'lightblue', lwd = 4)

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
           seed = 123) {
    
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
                   
                   t1 <- Sys.time()
                   for (day in 1:2) { # number of pollination days (blueberry days receptivity)
                     
                     message(paste('Starting day', day, 'of pollination'))
                     
                     for (trips_bee in seq_along(hive)) { # trips per bee
                       
                       flower_id <- sample.int(sum_flowers, hive[trips_bee]) # visited flowers
                       
                       for (k in seq_along(pla)) flowers[flower_id[k]] <- flowers[flower_id[k]] + 1
                       
                     }
                     
                   }
                   message(paste('Execution time:'))
                   print(Sys.time() - t1)
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
                   
                   t1 <- Sys.time()
                   for (day in 1:2) { # number of pollination days (blueberry days receptivity)
                     
                     message(paste('Starting day', day, 'of pollination'))
                     
                     for (trips_bee in seq_along(hive)) { # trips per bee
                       
                       flower_id <- sample.int(sum_flowers, hive[trips_bee]) # visited flowers
                       
                       for (k in seq_along(pla)) flowers[flower_id[k]] <- flowers[flower_id[k]] + 1
                       
                     }
                     
                   }
                   message(paste('Execution time:'))
                   print(Sys.time() - t1)
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
    
    return(plants)
    
  }


t1 <- Sys.time()
p <- simulated_visits(p_ha = p_01ha,
                      flowers_plant = total_flowers, 
                      visits_bee = visits_day, 
                      bees_hive = hives_ha(20), 
                      hive_aggregate = T)
Sys.time() - t1

sum(p$Hive1$plant1)

sum(p$Hive2$plant1) - sum(p$Hive1$plant1)

sum(p$Hive7$plant1) - sum(p$Hive6$plant1)

length(p)

par(mfrow = c(2, 2), mar = c(4.2, 4.2, 1.5, 1.5))

plot(density(p$Hive1[[1]]), ylim = c(0, 4.5), xlim = c(-0.5, 5),
     xlab = paste('Honeybee visits per flower\n (density:', 
                  1, ' hive per ha)'), 
     lwd = 0.3, main = '', col = 1)

for (j in 1:100) {
  lines(density(p$Hive1[[j]]), col = 1, lwd = 0.3) 
}


plot(density(p$Hive7[[1]]), ylim = c(0, 0.5), xlim = c(-0.5, 11),
     xlab = paste('Honeybee visits per flower\n (density:', 
                  7, ' hive per ha)'), 
     lwd = 0.3, main = '', col = 2)

for (j in 1:100) {
  lines(density(p$Hive7[[j]]), col = 2, lwd = 0.3) 
}

plot(density(p$Hive14[[1]]), ylim = c(0, 0.25), xlim = c(-0.5, 17),
     xlab = paste('Honeybee visits per flower\n (density:', 
                  14, ' hive per ha)'), 
     lwd = 0.3, main = '', col = 3)

for (j in 1:100) {
  lines(density(p$Hive14[[j]]), col = 3, lwd = 0.3) 
}

plot(density(p$Hive20[[1]]), ylim = c(0, 0.2), xlim = c(-0.5, 20),
     xlab = paste('Honeybee visits per flower\n (density:', 
                  20, ' hive per ha)'), 
     lwd = 0.3, main = '', col = 4)

for (j in 1:100) {
  lines(density(p$Hive20[[j]]), col = 4, lwd = 0.3) 
}

par(mfrow = c(1, 1))


plot(density(unlist(lapply(p$Hive20, median), use.names = F)), 
     xlim = c(0, 40), main = '', lwd = 0, ylim = c(0, 0.25), 
     xlab = 'Number of honeybee visits per flower')
for (i in seq_along(p)) {
  lines(density(unlist(lapply(p[[i]], median), use.names = F)), 
        lwd = 0.5, col = i)
}

lapply(p$output$Hive3, 
       FUN = function(x) quantile(x))

filt <- unlist(lapply(p$output$Hive3, 
                      FUN = function(x) quantile(x)[2] > 2), 
               use.names = F)

outpu <- p$output$Hive2[filt]

record <- p$cont[[7]][filt]

p$output$Hive1$plant3

plot(density(outpu$plant3))
record$plant3
length(p$cont[[4]])

vis_cult <- 
  lapply(p, FUN = 
           function(x) {
             
             t <- unlist(lapply(x, median), use.names = F)
             
             tibble(plant = paste('plant', 1:length(t), sep = ''), 
                    visits = t)
             
           })

for (i in seq_along(vis_cult)) {
  vis_cult[[i]]$n_hives <- paste(i)
}

vis_cult <- do.call('rbind', vis_cult)

vis_cult |> 
  group_by(n_hives) |> 
  transmute(mu = median(visits), 
            li = quantile(visits, 0.025),
            ls = quantile(visits, 0.975)) |> 
  unique() |> 
  ggplot(aes(fct_reorder(n_hives, ls), mu, ymin = li, ymax = ls)) +
  geom_errorbar(width = 0) +
  geom_point() +
  labs(x = 'Hive density per ha', y = 'Average floral visits per plant') +
  theme_bw() +
  theme(panel.grid = element_blank())

t <- lapply(p, FUN = 
              function(x) {
                lapply(x, quantile)
              })

t$Hive7

# ====== 5. Hoenybee pollen deposition ====
apisSVP <- readRDS('honeybee_svpd_data.rds')

apisSVP <- apisSVP[-which.max(apisSVP)] # removing value with 155 pollen deposition

plot(density(rnbinom(1e4, size = 2, mu = 15)), 
     main = expression(mu['honeybee pollen deposition'~'~'~'NegBinom(size = 2, mu = 15)']), 
     xlab = 'Single visit pollen deposition')

dat_apis <- list(visits = apisSVP, 
                 N = length(apisSVP))

cat(file = 'model_apis_SVP.stan', 
    '
    data{
      int N;
      array[N] int visits;
    }
    
    parameters{
      real mu;
      real<lower = 0> phi;
    }
    
    model {
      real lambda;
      mu ~ normal(2.5, 0.5);
      phi ~ exponential(1);
      
      lambda = exp(mu);
      visits ~ neg_binomial_2(lambda, phi);
    }
    ')

file <- paste(getwd(), '/model_apis_SVP.stan', sep = '')
fit_apis_svp <- cmdstan_model(file, compile = T)

model_apis_svp <- 
  fit_apis_svp$sample(
    data = dat_apis, 
    chains = 3, 
    iter_sampling = 2e3, 
    iter_warmup = 500, 
    parallel_chains = 3, 
    thin = 3, 
    refresh = 200, 
    seed = 123
  )

out_put_apis_svp <- model_apis_svp$summary()

model_apis_svp$summary()

par(mfrow = c(1, 3))
for (i in 1:3) trace_plot(model_apis_svp, out_put_apis_svp$variable[i], 3)
par(mfrow = c(1, 1))

post_apis_svp <- model_apis_svp$draws(format = 'df')[-1]
post_apis_svp$mu <- exp(post_apis_svp$mu)

pp_check_svp <- 
  sapply(1:500, FUN = 
           function(x) {
             
             rnbinom(1e3, mu = post_apis_svp$mu[x], 
                     size = post_apis_svp$phi[x])
             
           })

plot(density(pp_check_svp[1, ]), lwd = 0.1, main = '', 
     xlab = 'Single visit pollen deposition', ylim = c(0, 0.03))
for (i in 1:100) lines(density(pp_check_svp[i, ]), lwd = 0.1)
lines(density(dat_apis$visits), col = 'red')

set.seed(1234)
sim_visits <- 
  apply(pp_check_svp[c(sample(1:ncol(pp_check_svp), size = 50, replace = F)),], 1, 
        simplify = 'list', FUN = 
          function(x) {
            
            df_2 <- tibble()
            
            for (j in 1:10) {
              
              sum_vis <- vector("double", 2000)
              for (i in 1:2000) {
                sum_vis[[i]] <- sum(sample(x, size = j, T))
              }
              
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
plot(NULL, xlim = c(0, 20), ylim = c(0, 100))
for (i in 1:50) curve(alpha[i] + beta[i]*x, add = T, lwd = 0.3)

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

output_mod_slope <- mod_svp_slope$summary()
output_mod_slope

post_svp_slope <- mod_svp_slope$draws(format = 'df')
colnames(post_svp_slope)
#for (i in 2:3) post_svp_slope[[i]] <- exp(post_svp_slope[[i]])

par(mfrow = c(2, 2))
for (i in 1:4) trace_plot(mod_svp_slope, output_mod_slope$variable[i], 3)
par(mfrow = c(1, 1))

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

plot(NULL, lwd = 0.1, main = '', 
     xlab = 'Single visit pollen deposition', 
     xlim = c(0, 600), ylim = c(0, 0.007))
for (i in 1:50) lines(density(ppcheck_slope[i, ]), lwd = 0.1)
lines(density(sim_visits$poll_ac), col = 'red')

# here I have to fit the poisson model and then 


# ========= 6. Visit to pollen ========

pollen_fruit <- readRDS('datos_experimento.rds')
pollen_fruit <- pollen_fruit[, c("tratamiento", "fruto_diam", 
                                 "carga_poli", "carga_poli2")]

quantile(pollen_fruit$carga_poli)
mu_asintota <- round(((370-187)/2) + 187)
asymptote <- rnorm(2001, mu_asintota, 15)
slope <- post_svp_slope$beta



plot(NULL, xlim = c(0, 30), ylim = c(0, 320), xlab = 'Honeybee visits', 
     ylab = 'Single visit pollen deposition')
for (i in 1:200) curve(asymptote[i]*(1-exp(-slope[i]*x)), 
                       add = T, lwd = 0.1)
curve(mean(asymptote)*(1-exp(-mean(slope)*x)), 
      add = T, lwd = 2, col = 'red')

pollen_deposition_fun <- function(x, mu_est = T, n_posterior = 10) {
  if (mu_est) {
    a <- mean(asymptote)
    b <- mean(slope)
    
    return(round(a * (1-exp(-b*x))))
  } else {
    a <- sample(asymptote, n_posterior, T)
    b <- sample(slope, n_posterior, T)
    
    return(round(a * (1-exp(-b*x))))
  }
}

pollen_deposition <- 
  lapply(p, FUN = 
           function(x) {
             lapply(x, pollen_deposition_fun)
           })


par(mfrow = c(2, 2))

plot(density(pollen_deposition$Hive1[[1]]), ylim = c(0, 0.2), 
     xlim = c(-10, 159),
     xlab = paste('Pollen deposition per flower', 
                  '(1 hive per ha)'), 
     lwd = 0.3, main = '', col = 1)

for (j in 1:100) {
  lines(density(pollen_deposition$Hive1[[j]]), col = 1, lwd = 0.3) 
}


plot(density(pollen_deposition$Hive7[[1]]), ylim = c(0, 0.05), 
     xlim = c(-10, 300),
     xlab = paste('Pollen deposition per flower', 
                  '(7 hives per ha)'), 
     lwd = 0.3, main = '', col = 2)

for (j in 1:100) {
  lines(density(pollen_deposition$Hive7[[j]]), col = 2, lwd = 0.3) 
}


plot(density(pollen_deposition$Hive14[[1]]), ylim = c(0, 0.1), 
     xlim = c(-10, 300),
     xlab = paste('Pollen deposition per flower', 
                  '(14 hives per ha)'), 
     lwd = 0.3, main = '', col = 3)

for (j in 1:100) {
  lines(density(pollen_deposition$Hive14[[j]]), col = 3, lwd = 0.3) 
}

plot(density(pollen_deposition$Hive20[[1]]), ylim = c(0, 0.1), 
     xlim = c(-10, 300),
     xlab = paste('Pollen deposition per flower', 
                  '(20 hives per ha)'), 
     lwd = 0.4, main = '', col = 3)

for (j in 1:100) {
  lines(density(pollen_deposition$Hive20[[j]]), col = 4, lwd = 0.3) 
}

par(mfrow = c(1, 1))



# ============= 7. Pollen to fruit =====

optimal_model <- bf(fruto_diam ~ carga_poli + carga_poli2 + tratamiento)

get_prior(optimal_model, data = pollen_fruit, family = "Gaussian")

prior_optMod <- c(set_prior('normal(5, 3)', class = 'Intercept'),
                  set_prior('normal(0, 1)', class = 'b', coef = 'tratamientocs'),
                  set_prior('normal(0, 1)', class = 'b', coef = 'tratamientol'),
                  set_prior('normal(0, 1)', class = 'b', coef = 'tratamientols'),
                  set_prior('normal(0, 0.2)', class = 'b', coef = 'carga_poli'),
                  set_prior('normal(0, 0.2)', class = 'b', coef = 'carga_poli2'),
                  set_prior('exponential(1)', class = 'sigma'))

optim_model <- brm(formula = optimal_model, data = modelo, family = gaussian,
                   prior = prior_optMod, future = T, chains = 3,
                   thin = 3, iter = 25e3, warmup =  500, seed = 222)

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








