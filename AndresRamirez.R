paquetes <- c("ggplot2","fitdistrplus","brms",
              "future","GGally",  "ggeffects",
              "performance", "cowplot",
              "ggdist", "arm", "DHARMa","qqplotr",
              "reshape2", "bayesplot", "parallel", 
              "doBy","ggbreak", "ggridges", "ggdist",
              "pROC", 'tidyverse', 'magrittr', 'emmeans')


for (i in 1:length(paquetes)) require(paquetes[[i]], character.only = T)

# ============== Pregunta 1 ===========

p1 <- read.csv('examen_2022/pregunta1.csv', stringsAsFactors = T)

str(p1)
summary(p1)

p1 %$% aggregate(drupes ~ site + Treatment, FUN = length)


pois <- fitdist(p1$drupes, 'pois')
nbin <- fitdist(p1$drupes, 'nbinom')
cdfcomp(list(pois, nbin), main = NULL, plotstyle = 'ggplot')
qqcomp(list(pois, nbin), main = NULL, plotstyle = 'ggplot')



p1$cm_ashF <- as.factor(as.character(p1$cm_ash))

p1 |> 
  ggplot(aes(fct_reorder(site, drupes, .desc = T), drupes)) +
  geom_boxplot(aes(fill = Treatment),  alpha = 0.5) +
  scale_fill_manual(values = c('cyan4', 'tan1')) +
  geom_jitter(aes(color = Treatment), 
              position = position_dodge(width = 0.75)) +
  scale_color_manual(values = c('cyan4', 'tan1')) +
  labs(x = NULL, y = 'Número de drupas') +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank()) +
  geom_hline(yintercept = mean(p1$drupes), linetype = 3)


p1 |> 
  ggplot(aes(fct_reorder(cm_ashF, drupes), drupes, fill = Treatment)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c('cyan4', 'tan1')) + 
  geom_jitter(aes(color = Treatment), 
              position = position_dodge(width = 0.75)) +
  scale_color_manual(values = c('black', 'black')) +
  theme_bw() +
  labs(x = 'Cantidad de ceniza en el suelo', y = 'Número de drupas') +
  theme(panel.grid = element_blank(),
        legend.position = 'top') +
  stat_summary(fun = 'median', geom = 'point', position = position_dodge(width = 0.3)) +
  stat_summary(aes(group = 1, color = Treatment), fun = 'median', geom = 'line', 
               position = position_dodge(width = 0.3), color = 'red', linetype = 2)

p1 %$% aggregate(drupes ~ site + Treatment + cm_ash, FUN = length) # justamente los
                            # sitios varian en su nivel de ceniza volc[anica]

# definir las previas

get_prior(formula = drupes ~ Treatment*cm_ashF + (1|site), 
          data = p1, family = negbinomial)

# tenemos un grupo de refernecia freely.pol:cm_ashF0 y los demas son
# diferenciales de intercepto

modP1 <- bf(drupes ~ Treatment*cm_ashF + (1|site))

priors_p1m1 <- c(set_prior('normal(log(50), log((80-40)/4))', class = 'Intercept'),
                 set_prior('normal(0, 10)', class = 'b'),
                 set_prior('normal(0, 10)', class = 'sd', group = 'site'))

previa_refPM1 <- 
  tibble(x = 0:130,
       y = dnbinom(0:130, size = 10, mu = 50)) |> 
  ggplot(aes(x, y)) +
  geom_line(color = 'cyan4') +
  labs(x = expression(mu[referencia]), y = NULL) +
  theme_bw() +
  ggtitle(expression('BiNeg('~mu~'=50,'~theta~'=10)')) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14))
  
previa_difPM1 <- 
  tibble(x = c(-40, 50)) |> 
  ggplot(aes(x)) +
  stat_function(fun = dnorm, n = 1000,
                args = list(mean = 0, sd = 10), 
                color = 'cyan4') +
  labs(x = 'Diferencial de interceptos', y = NULL) +
  theme_bw() +
  ggtitle(expression('Norm('~mu~'=0,'~'SD'~'=10)')) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14))

plot_grid(previa_refPM1, previa_difPM1)  

# Previa para el intercepto segun Saez et al 2014 and Aizen et al 2014.
plot(density(rgamma(1000, shape = 1)))
plot(density(rlnorm(1000, 0, 1)))

mod1_p1 <- brm(formula = modP1, data = p1, family = negbinomial(),
               prior = priors_p1m1, warmup = 1500, future = T,
               chains = 3, iter = 4000, thin = 3)
summary(mod1_p1)
saveRDS(mod1_p1, 'mod1_p1.rds')
mod1_p1 <- readRDS('mod1_p1.rds')

mcmc_dens_overlay(mod1_p1, regex_pars = c("^b", "sigma"))
mcmc_dens_chains(mod1_p1, regex_pars = c("^b", "sigma")) 
mcmc_trace(mod1_p1, regex_pars = c("^b", "sigma"))

mcmc_acf(mod1_p1, regex_pars = c("^b", "sigma"))

bayes_R2(mod1_p1)

posteriorM1P1 <- as_draws_df(mod1_p1) # postriores de los coeficientes

# verificar las previas en general
previa_postm1 <- 
  ggplot() + # previa intercepto y posterior
  geom_density(data = posteriorM1P1, aes(x = exp(b_Intercept)), 
               linetype = 3, color = 'cyan4') +
  geom_line(data = tibble(x = 0:130,
                          y = dnbinom(0:130, size = 10, mu = 50)),
            aes(x, y), color = 'cyan4') +
  scale_x_continuous(limits = c(-0.5, 130)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = expression(mu[ref]), y = NULL)

previa_posteriotm1_2 <- 
  previa_difPM1 + 
  geom_density(data = posteriorM1P1, aes(x = b_Treatmentsuppl.pl, y = ..scaled..), 
               linetype = 2, color = 'dimgray') +
  geom_density(data = posteriorM1P1, aes(x = b_cm_ashF10, y = ..scaled..), 
               linetype = 3, color = 'dimgray') +
  geom_density(data = posteriorM1P1, aes(x = b_cm_ashF3, y = ..scaled..), 
               linetype = 4, color = 'dimgray') +
  geom_density(data = posteriorM1P1, aes(x = `b_Treatmentsuppl.pl:cm_ashF10`, y = ..scaled..), 
               linetype = 5, color = 'dimgray') +
  geom_density(data = posteriorM1P1, aes(x = `b_Treatmentsuppl.pl:cm_ashF3`, y = ..scaled..), 
               linetype = 6, color = 'dimgray') +
  labs(x = 'Dif. medias') +
  theme(title = element_blank(), 
        axis.title.x = element_text(size = 11))

zoom_posteriorm1 <- 
  ggplot() + 
  geom_density(data = posteriorM1P1, aes(x = b_Treatmentsuppl.pl, y = ..scaled..), 
               linetype = 2, color = 'dimgray') +
  geom_density(data = posteriorM1P1, aes(x = b_cm_ashF10, y = ..scaled..), 
               linetype = 3, color = 'dimgray') +
  geom_density(data = posteriorM1P1, aes(x = b_cm_ashF3, y = ..scaled..), 
               linetype = 4, color = 'dimgray') +
  geom_density(data = posteriorM1P1, aes(x = `b_Treatmentsuppl.pl:cm_ashF10`, y = ..scaled..), 
               linetype = 5, color = 'dimgray') +
  geom_density(data = posteriorM1P1, aes(x = `b_Treatmentsuppl.pl:cm_ashF3`, y = ..scaled..), 
               linetype = 6, color = 'dimgray') +
  theme_bw() +
  labs(x = 'Dif. medias') +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank())
  
  
plot_grid(previa_postm1, previa_posteriotm1_2,
            zoom_posteriorm1, ncol = 3, 
          labels = paste('(', letters[1:3],')', sep = ''),
          label_size = 10)
  
  
post_response_scale <- posteriorM1P1

post_response_scale$free_ash0 <- exp(posteriorM1P1$b_Intercept)
post_response_scale$supp_ash0 <- exp(posteriorM1P1$b_Intercept + 
                                       posteriorM1P1$b_Treatmentsuppl.pl)
post_response_scale$free_ash3 <- exp(posteriorM1P1$b_Intercept +
                                       posteriorM1P1$b_cm_ashF3)
post_response_scale$supp_ash3 <- exp(posteriorM1P1$b_Intercept +
                                       posteriorM1P1$b_Treatmentsuppl.pl +
                                       posteriorM1P1$b_cm_ashF3 +
                                       posteriorM1P1$`b_Treatmentsuppl.pl:cm_ashF3`)
post_response_scale$free_ash10 <- exp(posteriorM1P1$b_Intercept +
                                        posteriorM1P1$b_cm_ashF10)
post_response_scale$supp_ash10 <- exp(posteriorM1P1$b_Intercept +
                                        posteriorM1P1$b_Treatmentsuppl.pl +
                                        posteriorM1P1$b_cm_ashF10 +
                                        posteriorM1P1$`b_Treatmentsuppl.pl:cm_ashF3`)

post_response_scale <- post_response_scale[, grep('_ash[0-3]', colnames(post_response_scale))]

post_response_scale <- gather(post_response_scale)

post_response_scale$ash <- ifelse(grepl('h0$', post_response_scale$key), 'ceniza 0',
                                  ifelse(grepl('h3$', post_response_scale$key), 'ceniza 3', 
                                         'ceniza 10'))
post_response_scale$poll <- ifelse(grepl('^f', post_response_scale$key), 'libre',
                                   'suplementada')

post1 <- 
  post_response_scale |> 
  ggplot(aes(value, linetype = key)) +
  geom_density(color = 'cyan4') +
  scale_x_continuous(limits = c(-5, 200)) +
  theme_bw() +
  labs(x = expression(mu[drupas])) +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.7))

post2 <- 
  post_response_scale |> 
  ggplot(aes(value, linetype = poll)) +
  geom_density(color = 'cyan4') +
  scale_x_continuous(limits = c(-5, 200)) +
  theme_bw() +
  labs(x = expression(mu[drupas])) +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.7)) +
  geom_vline(xintercept = c(63.21434, 66.60), linetype = 1:2)


post3 <- 
  post_response_scale |> 
  ggplot(aes(value, linetype = ash)) +
  geom_density(color = 'cyan4') +
  scale_x_continuous(limits = c(-5, 200)) +
  theme_bw() +
  labs(x = expression(mu[drupas])) +
    scale_linetype_manual(values = c(1:3)) +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.7))

plot_grid(post2, post3)

post_response_scale %$% aggregate(value ~ poll, FUN = median)


sum((post_response_scale[post_response_scale$poll == 'suplementada' &
                           post_response_scale$ash == 'ceniza 10' , ]$value >= 
       quantile(post_response_scale$value, probs = c(0.75)))) /
  nrow(post_response_scale[post_response_scale$poll == 'suplementada' &
                             post_response_scale$ash == 'ceniza 10', ]) -
  sum((post_response_scale[post_response_scale$poll == 'libre' &
                             post_response_scale$ash == 'ceniza 10', ]$value >= 
         quantile(post_response_scale$value, probs = c(0.75)))) /
  nrow(post_response_scale[post_response_scale$poll == 'libre' &
                             post_response_scale$ash == 'ceniza 10', ])

sum((post_response_scale[post_response_scale$ash == 'ceniza 10', ]$value >= 
       quantile(post_response_scale$value, probs = c(0.75)))) /
  nrow(post_response_scale[post_response_scale$ash == 'ceniza 10', ]) -
  sum((post_response_scale[post_response_scale$ash == 'ceniza 0', ]$value >= 
         quantile(post_response_scale$value, probs = c(0.75)))) /
  nrow(post_response_scale[post_response_scale$ash == 'ceniza 0', ])




# diapo 16 practico 2
# calcular la diferencia entre grupos y calcular una distribuci[on]
# posterior de las diferencias entre grupos

conditional_effects(mod1_p1)

ppmp1 <- pp_check(mod1_p1, ndraws = 1000)

ppmp1 <- ppmp1$data

levels(ppmp1$rep_label)

ppmp1 <- ppmp1[ppmp1$rep_label != 'italic(y)', ]

ggplot() +
  geom_density(data = ppmp1, 
               mapping = aes(x = value, color = rep_label),
               alpha = 0.5, size = 0.5) +
  geom_density(data = p1, aes(drupes), color = 'red') +
  theme_bw() +
  scale_color_manual(values = rep('black', nrow(ppmp1))) +
  labs(x = 'Número de drupas') +
  theme(legend.position = 'none',
        panel.grid = element_blank())
  

post_predm1p1 <- predict(mod1_p1, ndraws = 1000, summary = F)

post_dfp1 <- data.frame()

for (i in 1:nrow(post_predm1p1)) {
  post_dfp1 <- rbind(post_dfp1,
                     data.frame(sim = rep(paste(i), 
                                          ncol(post_predm1p1)),
                                data = post_predm1p1[i, 1:354]))
}

post_dfp1 |> 
  ggplot() +
  geom_density(aes(data, color = sim)) + 
  geom_density(data = p1, aes(drupes)) +
  theme(legend.position = 'none') # otra forma de generar la distribuci[on]
                                  # predictiva posterior


# residual analysis 

q_resP1M1 <- createDHARMa(simulatedResponse = t(post_predm1p1),
                          observedResponse = p1$drupes,
                          fittedPredictedResponse = apply(post_predm1p1, 2, median),
                          integerResponse = T)

res_m1p1 <- qnorm(residuals(q_resP1M1))

res_m1p1 <- cbind(res = res_m1p1,
                  x = p1[, 3:5],
                  fitted = fitted(mod1_p1, ndraws = 1000)[, 1],
                  pareto = loo(mod1_p1, pointwise = T)$diagnostics$pareto_k)

r1p1 <- 
  res_m1p1 |> 
  ggplot(aes(fitted, res)) +
  geom_point(color = 'cyan4', alpha = 0.4) +
  geom_hline(yintercept = 0, color = 'cyan4', linetype = 3) +
  labs(x = 'Valores ajustados',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 7),
        text = element_text(size = 7))

r2p1 <- 
  res_m1p1 |> 
  ggplot(aes(sample = res)) +
  stat_qq_line(color = 'cyan4') +
  stat_qq_band(alpha = 0.4, fill = 'cyan4') +
  stat_qq_point(color = 'tan1', alpha = 0.3) +
  labs(x = 'Cuantiles teóricos',
       y = 'Cuantiles empíricos') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 7),
        text = element_text(size = 7))

r3p1 <- 
  res_m1p1 |> 
  ggplot(aes(x.Treatment, res)) +
  geom_boxplot(color = 'cyan4', fill = 'cyan4', alpha = 0.4) +
  geom_jitter(color = 'cyan4', width = 0.1, alpha = 0.5) +
  geom_hline(yintercept = 0, color = 'black', linetype = 3) +
  scale_x_discrete(labels = c(paste("Polinización", "libre", sep = "\n"),
                              paste('Polinización', "suplementada", sep = '\n'))) +
  labs(x = 'Porcentaje de arbustos',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 7),
        text = element_text(size = 7),
        axis.title.x = element_blank())

r4p1 <- 
  res_m1p1 |> 
  ggplot(aes(x.cm_ashF, res)) +
  geom_boxplot(color = 'cyan4', fill = 'cyan4', alpha = 0.4) +
  geom_jitter(color = 'cyan4', width = 0.1, alpha = 0.5) +
  geom_hline(yintercept = 0, color = 'black', linetype = 3) +
  #scale_x_discrete(labels = c(paste("Polinización", "libre", sep = "\n"),
  #                            paste('Polinización', "suplementada", sep = '\n'))) +
  labs(x = 'Cantidad de ceniza volcánica en el suelo (cm)',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 7),
        text = element_text(size = 7))



r5p1 <- 
  res_m1p1 |> # pareto
  ggplot(aes(1:nrow(res_m1p1), pareto)) +
  geom_point(color = 'cyan4', alpha = 0.4) +
  geom_hline(yintercept = 0.5, color = 'cyan4', linetype = 3) +
  geom_hline(yintercept = 0.7, color = 'cyan4', linetype = 3) +
  geom_hline(yintercept = 1, color = 'cyan4', linetype = 3) +
  labs(x = 'Dato',
       y = 'Pareto') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        title = element_text(size = 7),
        text = element_text(size = 7))

plot_grid(r1p1, r2p1, r3p1,
          r4p1, r5p1, nrow = 2, align = 'hv', 
          labels = letters[1:5])

# entender grupo de referencia y difereniales


# ========= Pregunta 2 ==============


p2 <- read.csv('examen_2022/pregunta2.csv', header = T, stringsAsFactors = T)

str(p2)

summary(p2)

# dado descarto variables que no se incluyen en el analisis

p2 <- p2[, !grepl('area', colnames(p2))]

levels(p2$Sp.name)

p2 %$% aggregate(Biomass ~ Genus + Species, FUN = length) # n por sp genero

# veamos como se relaciona nuestra variable respuesta Biomasa con las variables 
# predictoras por genero

genus <- split(p2, p2$Genus)

genus <- lapply(genus, FUN = 
                  function(x) {
                    ggcorr(x[, -c(1:3)], label = T, 
                           label_size = 2,
                           size = 2.5, hjust = 0.75) +
                      theme(legend.position = 'none',
                            title = element_text(size = 8))
                  })


for (i in seq_along(genus)) {
  genus[[i]] <- genus[[i]] + ggtitle(paste(names(genus)[[i]]))
}

plot_grid(genus$Dipterocarpus, genus$Dryobalanops,
          genus$Hopea, genus$Parashorea,
          genus$Shorea) 

# En todos los casos las correlaciones con entre biomasa y los predictores son 
# positivas

# ajuste a la distribuci[on de probabilidad
norm <- fitdist(p2$Biomass, 'norm')
lnorm <- fitdist(p2$Biomass, 'lnorm')
gamma <- fitdist(p2$Biomass, 'gamma')
cdf <- cdfcomp(list(norm, lnorm, gamma), plotstyle = 'ggplot')
qq <- qqcomp(list(norm, lnorm, gamma), plotstyle = 'ggplot') # mejor gamma
gridExtra::grid.arrange(cdf, qq, ncol = 1)

# antes de ajustar el modelo estandarizo los predictores para comparra sus 
# coeficientes en el modelo


# previas

plot(seq(0, 3, 0.01), dlnorm(seq(0, 3, 0.01), 0, 1))

theta <- 
  data.frame(x = c(0, 4)) |> 
  ggplot(aes(x)) +
  stat_function(fun = dcauchy, n = 1000, args = list(location=0, scale=3), 
                color = 'cyan4') +
  labs(x = expression(theta), y = 'Dens. Prob.') +
  ggtitle('lnorm(mean = 0, SD = 1)') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        title = element_text(size = 7))

beta0 <- 
  data.frame(x = c(-3, 3)) |> 
  ggplot(aes(x = x)) +
  stat_function(fun = dnorm, n = 1000, args = list(0, 1), 
                color = 'cyan4') +
  labs(x = expression(beta[0]), y = '') +
  ggtitle('norm(mean = 0, SD = 1)') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        title = element_text(size = 7))

beta_i <- 
  data.frame(x = c(-4, 4)) |> 
  ggplot(aes(x = x)) +
  stat_function(fun = dnorm, n = 1000, args = list(1, 1), 
                color = 'cyan4') +
  geom_area(stat = 'function',
            fun = dnorm, args = list(1, 1),
            color = 'cyan4',
            xlim = c(qnorm(0.000001, 1, 1),
                     qnorm(pnorm(0, 1, 1), 1, 1)), 
            fill = 'cyan4', alpha = 0.5) +
  labs(x = expression(beta[i]), y = '') +
  ggtitle('norm(mean = 1, SD = 1)') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        title = element_text(size = 7))

plot_grid(theta, beta0, beta_i, ncol = 3)

p2s <- p2

p2s[] <- lapply(p2, 
                function(x) {
                  if (!is.factor(x)) as.vector(scale(x)) else(x)
                })

summary(p2s)

p2s$Biomass <- p2$Biomass

# definimos el modelo y las previas

mod_p2 <- bf(Biomass ~ Height + 
               leaves + Diam + (1|Species/Genus))

get_prior(formula = mod_p2, 
                    data = p2s, family = lognormal())

prio_p2 <- c(set_prior('normal(1, 1)', class = 'b'),
               set_prior('normal(0, 1)', class = 'Intercept'),
               set_prior("lognormal(0, 1)", class = 'sigma'))

m1_p2 <- brm(formula = mod_p2, data = p2s, family = lognormal(), prior = prio_p2,
           future = T, warmup = 1000, iter = 4000, chains = 3, thin = 3)
summary(m1_p2)

saveRDS(m1_p2, 'model_p2.rds')

m1_p2 <- readRDS('model_p2.rds')

bayes_R2(m1_p2)
cond_effp2 <- conditional_effects(m1_p2, type = 'response')
cond_effp2 <- lapply(cond_effp2, as_tibble)

f1 <- 
  ggplot(data = cond_effp2$Height,
       aes(Height, estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__),
              fill = 'tan1', alpha = 0.3) +
  geom_line(color = 'cyan4') +
  theme_bw() +
  labs(y = 'Biomasa') +
  theme(panel.grid = element_blank())


f2 <- 
  ggplot(data = cond_effp2$leaves,
       aes(leaves, estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__),
              fill = 'tan1', alpha = 0.3) +
  geom_line(color = 'cyan4') +
  theme_bw() +
  labs(y = '') +
  theme(panel.grid = element_blank())

f3 <- 
  ggplot(data = cond_effp2$Diam,
       aes(Diam, estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__),
              fill = 'tan1', alpha = 0.3) +
  geom_line(color = 'cyan4') +
  theme_bw() +
  labs(y = '') +
  theme(panel.grid = element_blank())

plot_grid(f1, f2, f3, ncol = 3, align = 'hv', labels = letters[1:3])


posterior_p2 <- as_tibble(as_draws_df(m1_p2, add_chains = T))

theta <- tibble(x = seq(0, 4, 0.01),
                dens = dlnorm(seq(0, 4, 0.01), 0, 1))

f1 <- ggplot() +
  geom_density(data = posterior_p2,
               aes(sigma, y = ..scaled..),
               color = 'cyan4', fill = 'cyan4', alpha = 0.5) +
  geom_line(data = theta,
            aes(x, dens), color = 'red') +
  labs(x = expression(theta), y = 'Dens. Prob.') +
  theme_bw() +
  theme(panel.grid = element_blank())

f4 <- ggplot() +
  geom_density(data = posterior_p2,
               aes(sigma, y = ..scaled..),
               color = 'cyan4', fill = 'cyan4', alpha = 0.5) +
  labs(x = expression(theta), y = 'Dens. Prob.') +
  theme_bw() +
  theme(panel.grid = element_blank())


b0 <- tibble(x = seq(-3, 3, 0.01),
                dens = dnorm(seq(-3, 3, 0.01), 0, 1))

f2 <- ggplot() +
  geom_density(data = posterior_p2,
               aes(b_Intercept, y = ..scaled..),
               color = 'cyan4', fill = 'cyan4', alpha = 0.5) +
  geom_line(data = b0,
            aes(x, dens), color = 'red') +
  labs(x = expression(beta[0]), y = '') +
  theme_bw() +
  theme(panel.grid = element_blank())

f5 <- ggplot() +
  geom_density(data = posterior_p2,
               aes(exp(b_Intercept), y = ..scaled..),
               color = 'cyan4', fill = 'cyan4', alpha = 0.5) +
  labs(x = expression(beta[0]), y = '') +
  theme_bw() +
  theme(panel.grid = element_blank())

b1 <- tibble(x = seq(-1, 3, 0.01),
             dens = dnorm(seq(-1, 3, 0.01), 1, 1))

f3 <- ggplot() +
  geom_density(data = posterior_p2,
               aes(b_Height, y = ..scaled..),
               color = 'cyan4', fill = 'cyan4', alpha = 0.5) +
  geom_density(data = posterior_p2,
               aes(b_leaves, y = ..scaled..),
               color = 'tan1', fill = 'tan1', alpha = 0.5) +
  geom_density(data = posterior_p2,
               aes(b_Diam, y = ..scaled..),
               color = 'gray3', fill = 'gray3', alpha = 0.5) +
  geom_line(data = b1,
            aes(x, dens), color = 'red') +
  labs(x = expression(beta[i]), y = '') +
  theme_bw() +
  theme(panel.grid = element_blank())

f6 <- ggplot() +
  geom_density(data = posterior_p2,
               aes(exp(b_Height), y = ..scaled..),
               color = 'cyan4', fill = 'cyan4', alpha = 0.5) +
  geom_density(data = posterior_p2,
               aes(exp(b_leaves), y = ..scaled..),
               color = 'tan1', fill = 'tan1', alpha = 0.5) +
  geom_density(data = posterior_p2,
               aes(exp(b_Diam), y = ..scaled..),
               color = 'gray3', fill = 'gray3', alpha = 0.5) +
  labs(x = expression(beta[i]), y = '') +
  theme_bw() +
  theme(panel.grid = element_blank())


plot_grid(f1, f2, f3,
          f4, f5, f6, ncol = 3, align = 'hv', labels = letters[1:6])

ppcheck_p2 <- predict(m1_p2, ndraws = 1000, summary = F)

post_df <- data.frame()

for (i in 1:nrow(ppcheck_p2)) {
  post_df <- rbind(post_df,
                   data.frame(sim = rep(paste(i), 290),
                              data = ppcheck_p2[i, 1:290]))
}


ggplot() +
  geom_density(data = post_df, 
               aes(data, color = sim), size = 0.1) +
  scale_color_manual(values = rep('cyan4', length(unique(post_df$sim)))) +
  geom_density(data = p2s, 
               aes(Biomass), color = 'tan1', size = 1) +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = 'Biomasa', y = 'Dens. Prob.')
  

q_resM1 <- createDHARMa(simulatedResponse = t(ppcheck_p2),
                        observedResponse = p2s$Biomass,
                        fittedPredictedResponse = apply(ppcheck_p2, 2, median),
                        integerResponse = T)
res_m1 <- data.frame(res = qnorm(residuals(q_resM1)))

res_m1 <- cbind(res_m1,
                p2s[, c("Height", "Diam", "leaves")],
                fitted = fitted(m1_p2, ndraws = 1000)[, 1], 
                pareto = loo(m1_p2, pointwise = T)$diagnostics$pareto_k)

f2 <- 
  res_m1 |> # qqplot
  ggplot(aes(sample = res)) +
  stat_qq_line(color = 'cyan4') +
  stat_qq_band(alpha = 0.4, fill = 'cyan4') +
  stat_qq_point(color = 'tan1', alpha = 0.3) +
  labs(x = 'Cuantiles teóricos',
       y = 'Cuantiles empíricos') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 7),
        text = element_text(size = 7))


f1 <- 
  res_m1 |> # fitted vs res
  ggplot(aes(fitted, res)) +
  geom_point(color = 'tan1', alpha = 0.3) +
  geom_hline(yintercept = 0, color = 'cyan4', linetype = 3) +
  labs(x = 'Valores ajustados',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 7),
        text = element_text(size = 7))

f3 <- 
  res_m1 |> # residuals vs variables explicativas
  ggplot(aes(Height, res)) +
  geom_point(color = 'tan1', alpha = 0.3) +
  geom_hline(yintercept = 0, color = 'cyan4', linetype = 3) +
  labs(x = 'Altura planta',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 7),
        text = element_text(size = 7))

f4 <- 
  res_m1 |> # residuals vs variables explicativas
  ggplot(aes(Diam, res)) +
  geom_point(color = 'tan1', alpha = 0.3) +
  geom_hline(yintercept = 0, color = 'cyan4', linetype = 3) +
  labs(x = 'Diámetro hoja',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        title = element_text(size = 7),
        text = element_text(size = 7))

f5 <- 
  res_m1 |> # residuals vs variables explicativas
  ggplot(aes(leaves, res)) +
  geom_point(color = 'tan1', alpha = 0.3) +
  geom_hline(yintercept = 0, color = 'cyan4', linetype = 3) +
  labs(x = 'Número de hojas',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        title = element_text(size = 7),
        text = element_text(size = 7))


f6 <- 
  res_m1 |> # pareto
  ggplot(aes(1:nrow(res_m1), pareto)) +
  geom_point(color = 'tan1', alpha = 0.3) +
  geom_hline(yintercept = 0.5, color = 'cyan4', linetype = 3) +
  geom_hline(yintercept = 0.7, color = 'cyan4', linetype = 3) +
  geom_hline(yintercept = 1, color = 'cyan4', linetype = 3) +
  labs(x = 'Dato',
       y = 'Pareto') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        title = element_text(size = 7),
        text = element_text(size = 7))


plot_grid(f1, f2, f3,
          f4, f5, f6, nrow = 2,
          align = 'hv', labels = letters[1:6])



# ======== Pregunta 3 ===============

p3 <- read.csv('examen_2022/pregunta3.csv', header = T, stringsAsFactors = T)

str(p3)
summary(p3)

f1 <- 
  p3 |> 
  ggplot(aes(PERSHRUB, RODENTSP)) +
  geom_smooth(method = 'glm', method.args = list(family = 'binomial'),
              se = F, color = 'tan1') +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'Arbustos área', y = 'Presencia de roedores')


f2 <- 
  p3 |> 
  ggplot(aes(DISTX, RODENTSP)) +
  geom_smooth(method = 'glm', method.args = list(family = 'binomial'),
              se = F, color = 'tan1') +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = 'Cañón más cercano (m)')


f3 <- 
  p3 |> 
  ggplot(aes(AGE, RODENTSP)) +
  geom_smooth(method = 'glm', method.args = list(family = 'binomial'),
              se = F, color = 'tan1') +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = 'Años desde fragmentación', y = 'Presencia de roedores')

f4 <- 
  p3 |> 
  ggplot(aes(AGE, PERSHRUB)) +
  geom_smooth(method = 'lm',
              se = F, color = 'tan1') +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = 'Años desde fragmentación')

plot_grid(f1, f2, f3, f4, ncol = 2, align = 'hv')


GGally::ggcorr(p3[, 1:3], label = T)

p3s <- as_tibble(scale(p3[, 1:3], center = T, scale = T))
p3s$RODENTSP <- p3$RODENTSP

modelo1 <- bf(RODENTSP ~ PERSHRUB + DISTX)

get_prior(modelo1, data = p3s, family = 'Binomial')

pb1 <- 
  tibble(x = c(-4, 4)) |> 
  ggplot(aes(x)) +
  stat_function(fun = dnorm, n = 1000, args = list(0.7, 1), 
                color = 'cyan4') +
  geom_area(stat = 'function',
            fun = dnorm, args = list(0.7, 1),
            fill = 'cyan4', alpha = 0.7,
            xlim = c(qnorm(0.00001, 0.7, 1),
                     qnorm(pnorm(0, 0.7, 1), 0.7, 1))) +
  geom_area(stat = 'function',
            fun = dnorm, args = list(0.7, 1),
            fill = 'tan1', alpha = 0.5, 
            xlim = c(qnorm(pnorm(0, 0.7, 1), 0.7, 1),
                     qnorm(0.00001, 0.7, 1, lower.tail = F))) +
  geom_label(label = paste('P(x <= 0)', round(pnorm(0, 0.7, 1), 2), sep = '\n'),
             x = -3, y = 0.2, size = 3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = expression(beta[1]), y = 'Dens. Prob.')

pb0 <- 
  tibble(x = c(-4, 4)) |> 
  ggplot(aes(x)) +
  stat_function(fun = dnorm, n = 1000, args = list(0, 1), 
                color = 'cyan4') +
  geom_area(stat = 'function',
            fun = dnorm, args = list(0, 1),
            fill = 'cyan4', alpha = 0.7,
            xlim = c(qnorm(0.00001, 0, 1),
                     qnorm(pnorm(0, 0, 1), 0, 1))) +
  geom_area(stat = 'function',
            fun = dnorm, args = list(0, 1),
            fill = 'tan1', alpha = 0.5, 
            xlim = c(qnorm(pnorm(0, 0, 1), 0, 1),
                     qnorm(0.00001, 0, 1, lower.tail = F))) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = expression(beta[0]), y = 'Dens. Prob.')

plot_grid(pb0, pb1, ncol = 2, align = 'hv')

plot(density(rnorm(100, -0.5, 1)))
priorM1 <- c(set_prior('normal(0.7, 1)', class = 'b', coef = 'PERSHRUB'),
             set_prior('normal(0, 1)', class = 'b', coef = 'DISTX'),
             set_prior('normal(0, 1)', class = 'Intercept'))

m1_ <- brm(modelo1, data = p3s, family = bernoulli, prior = priorM1,
          future = T, warmup = 1000, iter = 3000, chains = 3, thin = 3)
summary(m1_)
saveRDS(m1_, 'modelo1P3.rds')

m1_ <- readRDS('modelo1P3.rds')

priorM2 <- c(set_prior('normal(0, 3)', class = 'b', coef = 'PERSHRUB'),
             set_prior('normal(0, 3)', class = 'b', coef = 'DISTX'),
             set_prior('normal(0, 3)', class = 'Intercept')) # previas no informativas

m2_ <-  brm(modelo1, data = p3s, family = bernoulli, prior = priorM2,
           future = T, warmup = 1000, iter = 3000, chains = 3, thin = 3)
summary(m2_)
saveRDS(m2_, 'modelo2P3.rds')
m2_ <- readRDS('modelo2P3.rds')

r2_m1 <- bayes_R2(m1_, summary = F)
r2_m2 <- bayes_R2(m2_, summary = F)

rp3 <- rbind(tibble(r2 = r2_m1,
                    mod = rep('Pre. Inf.', nrow(r2_m1))),
             tibble(r2 = r2_m2,
                    mod = rep('Pre. No Inf.', nrow(r2_m1))))

rp3 |> 
  ggplot(aes(r2, color = mod, fill = mod)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = c('tan1', 'cyan4')) +
  scale_fill_manual(values = c('tan1', 'cyan4')) +
  theme_bw() +
  labs(x = expression(r^2), y = 'Dens. Prob.') +
  theme(panel.grid = element_blank(),
        legend.position = c(0.35, 0.8), 
        legend.title = element_blank())


tapply(rp3$r2, rp3$mod, function(x) data.frame(mean = mean(x),
                                           sd = sd(x),
                                           q1 = quantile(x, probs = c(0.25))[[1]],
                                           q3 = quantile(x, probs = c(0.75))[[1]]))


variables(m1_)

mcmc_dens_overlay(m1_, regex_pars = c('^b_')) +
               geom_density()

mcmc_trace(m1_)

mcmc_acf(m1_, regex_pars = c('^b_')) +
  geom_line()

mcmc_areas(m1_, regex_pars = c('^b_')) +
  geom_line()

mcmc_areas_data(m1_)

ppcM1 <- as_tibble(as_draws_df(m1_, add_chain = T)) # distribuci[on posterior de todos los par[ametros]]

ppcM1$mod <- rep('Pre. Inf.', nrow(ppcM1))

logit_rev <- function(x) exp(x)/(1 + exp(x))

for (i in 1:3) ppcM1[[i]] <- logit_rev(ppcM1[[i]])

b1 <- dnorm(seq(from = -4, to = 4, by = 0.001), 0.7, 1)
b1 <- tibble(x = seq(from = -4, to = 4, by = 0.001),
             dens = b1)

previa_postP3 <- 
  ggplot(data = data.frame(x = c(-3, 4)), aes(x)) +
  geom_density(data = ppcM1, aes(b_PERSHRUB), color = 'cyan4',
               fill = 'cyan4', alpha = 0.5) +
  stat_function(fun = dnorm, n = 1000, args = list(0.7, 1), 
                color = 'red', alpha = 0.5) +
  theme_bw() +
  labs(x = expression(beta[1]), y = 'Dens. Prob.') +
  theme(panel.grid = element_blank())
  
sum(logit_rev(ppcM1$b_PERSHRUB) <= 0)/nrow(ppcM1)

previa_postP3_ <- 
  ggplot(data = data.frame(x = c(-3, 4)), aes(x)) +
  geom_density(data = ppcM1, aes(b_DISTX), color = 'cyan4',
               fill = 'cyan4', alpha = 0.5) +
  stat_function(fun = dnorm, n = 1000, args = list(0, 1), 
                color = 'red', alpha = 0.5) +
  theme_bw() +
  labs(x = expression(beta[0]), y = 'Dens. Prob.') +
  theme(panel.grid = element_blank())

sum(ppcM1$b_DISTX <= 0)/nrow(ppcM1)



tibble(x = c(-4, 4)) |> 
  ggplot(aes(x)) +
  stat_function(fun = dnorm, n = 1000, args = list(0.7, 1), 
                color = 'cyan4') +
  geom_area(stat = 'function',
            fun = dnorm, args = list(0.7, 1),
            fill = 'cyan4', alpha = 0.7,
            xlim = c(qnorm(0.00001, 0.7, 1),
                     qnorm(0.00001, 0.7, 1, ))) +
  geom_area(stat = 'function',
            fun = dnorm, args = list(0.7, 1),
            fill = 'tan1', alpha = 0.5, 
            xlim = c(qnorm(pnorm(0, 0.7, 1), 0.7, 1),
                     qnorm(0.00001, 0.7, 1, lower.tail = F))) +
  geom_label(label = paste('P(x <= 0)', round(pnorm(0, 0.7, 1), 2), sep = '\n'),
             x = -3, y = 0.2, size = 3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = expression(beta[1]), y = 'Dens. Prob.')


tibble(x = c(-4, 4)) |> 
  ggplot(aes(x)) +
  stat_function(fun = dnorm, n = 1000, args = list(0, 1), 
                color = 'cyan4') +
  geom_area(stat = 'function',
            fun = dnorm, args = list(0, 1),
            fill = 'cyan4', alpha = 0.7,
            xlim = c(qnorm(0.00001, 0, 1),
                     qnorm(pnorm(0, 0, 1), 0, 1))) +
  geom_area(stat = 'function',
            fun = dnorm, args = list(0, 1),
            fill = 'tan1', alpha = 0.5, 
            xlim = c(qnorm(pnorm(0, 0, 1), 0, 1),
                     qnorm(0.00001, 0, 1, lower.tail = F))) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = expression(beta[0]), y = 'Dens. Prob.')



# residuales 

res <- as_tibble(residuals(m1_, type = 'pearson', ndraws = 1000, summary = T))
fitted_ <- fitted(m1_, ndraws = 1000, summary = T, scale = 'linear')

plot(fitted_[,1], res$Estimate); abline(h = 0) # es porque es una binomial, tranqui
qqnorm(res$Estimate); qqline(res)

res |> 
  ggplot(aes(sample = Estimate)) +
  stat_qq_line() +
  stat_qq_point() +
  stat_qq_band(alpha = 0.5)

# efectos condicionales 

cond_eff <- conditional_effects(m2_)
cond_eff$PERSHRUB
cond_eff$DISTX



# distribucion predictiva posterior
post_predM1 <- predict(m1_, ndraws = 1000, summary = F)

post_df <- data.frame()

for (i in 1:nrow(post_predM1)) {
  post_df <- rbind(post_df,
                   data.frame(sim = rep(paste(i), 25),
                              data = post_predM1[i, 1:25]))
}

post_predM2 <- predict(m2_, ndraws = 1000, summary = F)

post_df2 <- data.frame()

for (i in 1:nrow(post_predM2)) {
  post_df2 <- rbind(post_df2,
                    data.frame(sim = rep(paste(i), 25),
                               data = post_predM2[i, 1:25]))
}

f1 <- 
  post_df |> 
  ggplot(aes(data, color = sim)) +
  geom_density(size = 0.1) +
  geom_density(data = p3, aes(RODENTSP), color = 'red') +
  scale_color_manual(values = rep('black', 1000)) +
  theme_bw() +
  labs(x = "Observacion roedor", y = 'Density',
       title = 'Previas informativas') +
  theme(legend.position = 'none',
        title = element_text(size = 7)) 
  

f2 <- 
  post_df2 |> 
  ggplot(aes(data, color = sim)) +
  geom_density(size = 0.1) +
  geom_density(data = p3, aes(RODENTSP), color = 'red') +
  scale_color_manual(values = rep('black', 1000)) +
  theme_bw() +
  labs(x = "Observacion roedor", title = 'Previas no informativas') +
  theme(legend.position = 'none',
        title = element_text(size = 7),
        axis.title.y = element_blank())

plot_grid(f1, f2, ncol = 2, align = 'hv') # sensibilidad previas

presences <- function(x) mean(x == 1)
pp_check(m1_, ndraws = 1000, type = 'stat', 
         stat = 'presences') +
  geom_histogram() +
  theme(axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        legend.position = "none")


q_resM1 <- createDHARMa(simulatedResponse = t(post_predM1),
             observedResponse = p3s$RODENTSP,
             fittedPredictedResponse = apply(post_predM1, 2, median),
             integerResponse = T)
res_m1 <- data.frame(res = qnorm(residuals(q_resM1)))

res_m1 <- cbind(res_m1,
                p3s[, 1:3],
                fitted = fitted(m1_, ndraws = 1000)[, 1], 
                pareto = loo(m1_, pointwise = T)$diagnostics$pareto_k)

f2 <- 
  res_m1 |> # qqplot
  ggplot(aes(sample = res)) +
  stat_qq_line(color = 'cyan4') +
  stat_qq_band(alpha = 0.4, fill = 'cyan4') +
  stat_qq_point(color = 'tan1') +
  labs(x = 'Cuantiles teóricos',
       y = 'Cuantiles empíricos') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 7),
        text = element_text(size = 7))


f1 <- 
  res_m1 |> # fitted vs res
  ggplot(aes(fitted, res)) +
  geom_point(color = 'tan1') +
  geom_hline(yintercept = 0, color = 'cyan4', linetype = 3) +
  labs(x = 'Valores ajustados',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 7),
        text = element_text(size = 7))

f3 <- 
  res_m1 |> # residuals vs variables explicativas
  ggplot(aes(PERSHRUB, res)) +
  geom_point(color = 'tan1') +
  geom_hline(yintercept = 0, color = 'cyan4', linetype = 3) +
  labs(x = 'Porcentaje de arbustos',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 7),
        text = element_text(size = 7))

f4 <- 
  res_m1 |> # residuals vs variables explicativas
  ggplot(aes(DISTX, res)) +
  geom_point(color = 'tan1') +
  geom_hline(yintercept = 0, color = 'cyan4', linetype = 3) +
  labs(x = 'Distancia a caño más cercano',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        title = element_text(size = 7),
        text = element_text(size = 7))


f5 <- 
  res_m1 |> # pareto
  ggplot(aes(1:nrow(res_m1), pareto)) +
  geom_point(color = 'tan1') +
  geom_hline(yintercept = 0.5, color = 'cyan4', linetype = 3) +
  geom_hline(yintercept = 0.7, color = 'cyan4', linetype = 3) +
  geom_hline(yintercept = 1, color = 'cyan4', linetype = 3) +
  labs(x = 'Dato',
       y = 'Pareto') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        title = element_text(size = 7),
        text = element_text(size = 7))


plot_grid(f1, f2, f3,
          f4, f5, nrow = 2,
          align = 'hv')


# efectos condicionales

temp_plot <- ggpredict(m1_, terms = c("PERSHRUB"))
temp_plot

ggpredict(m1_, terms = c("PERSHRUB")) |> 
  ggplot(aes(x, predicted)) +
  geom_line(color = 'cyan4') +
  geom_line(aes(x, conf.low), linetype = 3, color = 'cyan4') +
  geom_line(aes(x, conf.high), linetype = 3, color = 'cyan4') +
  geom_jitter(data = p3s, 
             aes(PERSHRUB, RODENTSP), height = 0.02, color = 'tan1') +
  labs(y = 'Probabilidad de capturar roedores',
       x = 'porcentaje de bosque') +
  theme_bw() +
  theme(panel.grid = element_blank())

ggpredict(m1_, terms = c("DISTX")) |> 
  ggplot(aes(x, predicted)) +
  geom_line(color = 'cyan4') +
  geom_line(aes(x, conf.low), linetype = 3, color = 'cyan4') +
  geom_line(aes(x, conf.high), linetype = 3, color = 'cyan4') +
  geom_jitter(data = p3s, 
              aes(PERSHRUB, RODENTSP), height = 0.02, color = 'tan1') +
  labs(y = 'Distancia al ca;o m[as cercano',
       x = 'porcentaje de bosque') +
  theme_bw() +
  theme(panel.grid = element_blank())















