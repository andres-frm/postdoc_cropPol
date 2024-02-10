pks <- c('tidyverse', 'rethinking', 'rstan', 'magrittr', 'cmdstanr',
         'ggdag', 'dagitty')
sapply(pks, library, character.only = T)


sim_mod <- 
  dagify(Y ~ FW + nF, 
         nF ~ FloD + Pfs, 
         FW ~ PollD,
         Pfs ~ nVis,
         PollD ~ nVis + wild,
         nVis ~ apis + percF, 
         apis ~ propF + apisN + timeF,
         BHd ~ BHs,
         propF ~ BHs + BHd,
         exposure = 'BHd', 
         outcome = 'Y')

coordinates(sim_mod) <- 
  list(
    y = c(BHs = 1, propF = 0, BHd = -1, apis = 0, 
          apisN = 1, timeF = -1,  nVis = 0, PollD = -1, wild = -2, Pfs = 1,
          nF = 1, FW = -1, FloD = 2, Y = 0, percF = 1), 
    
    x = c(BHs = 1, propF = 2, BHd = 1, apis = 3, 
          apisN = 3,  timeF = 3, nVis = 4, PollD = 5, wild = 5, Pfs = 5, 
          nF = 6, FW = 6, FloD = 6, Y = 7, percF = 4)
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

hive_ha <- 15

time_visit <- rgamma(2e3, 
                     shape = (17^2)/8, 
                     rate = 17/8)
plot(density(time_visit))

prop_foragers <- rbeta(1e3, 4.5, 8)
plot(density(rbeta(1e3, 4.5, 8)))

beehive_strength <- rnbinom(hive_ha, size = 10, mu = 1e4)

time <- vector('double', length = 1)
vis <- 0

visits_bee <- vector('double', beehive_strength[1])

foraging_time <- sample(foraging_time, beehive_strength[1], T)

t1 <- Sys.time()
for (i in 1:beehive_strength[1]) {
  
  time <- 0
  vis <- 0
  
  while (time < foraging_time[[i]]) {
    time <- time + sample(time_visit, 1)
    vis <- vis + 1
  }
  
  visits_bee[[i]] <- vis
}
Sys.time() - t1


plot(density(visits_bee))

t1 <- Sys.time()
visits_bee <- 
  sapply(1:beehive_strength[1], FUN = 
           function(x) {
             
             time <- 0
             vis <- 0
             
             while (time < foraging_time[[x]]) {
               time <- time + sample(time_visit, 1)
               vis <- vis + 1
             }
             
             vis
             
           })
Sys.time() - t1

plot(density(visits_bee))

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








