library(tidyverse)
library(magrittr)
library(lattice)
library(readxl)
library(brms)

VR <- as_tibble(read.csv('censosPRI_SJ.csv', sep = ';', dec = ','))
fruit_quality <- as_tibble(read.csv('fruit_qualityPRI_SJ.csv', sep = ';', dec = ','))

VR$tv <- VR$a.mellifera.fv. / VR$no.flo

summary(VR$hora)

VR$hora2 <- VR$hora^2

densityplot(VR$tv)

acti_mod <- lm(tv ~ hora + hora2, data = VR)
summary(acti_mod)

VR$tv_sim <- -0.292878 + (VR$hora*0.060952) + (VR$hora2*-0.002081)

a <- coef(acti_mod)[[3]]
b <- coef(acti_mod)[[2]]

x_value <- -b/(2*a) # cantidad de polen que maximiza el tamano del fruto
y_value <- (a*(x_value)^2) + (b*x_value) + coef(acti_mod)[1] # as?nto

asintota <- data.frame(x = x_value, 
                       y = y_value)

VR |> 
  ggplot() +
  geom_point(aes(hora, tv, color = finca), alpha = 0.6) +
  geom_line(aes(hora, tv_sim)) +
  geom_point(data = asintota, aes(x, y), color = 'red') +
  theme_bw() +
  labs(x = 'Hora', y = 'Tasa de visita Apis') +
  theme(panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman')) +
  geom_vline(xintercept = c(12.3, 15.3), linetype = 3)

VR <- VR[, c("fecha", "finca", "variedad", "lote",
             "planta", "no.flo", "a.mellifera.fv.")]

for (i in 1:(length(VR)-2)) VR[[i]] <- as.factor(VR[[i]])

VR$tv <- VR$a.mellifera.fv. / VR$no.flo

VR <- 
  VR |> 
  group_by(finca, variedad, lote, planta) |> 
  transmute(tv = mean(tv, na.rm = T))

fruit_quality <- 
  fruit_quality |> 
  group_by(finca, variedad, lote, planta) |> 
  transmute(diamF = mean(diamF, na.rm = T)) |> 
  unique()

for (i in 1:(length(fruit_quality) - 1)) fruit_quality[[i]] <- as.factor(fruit_quality[[i]])

VR$finca <- factor(VR$finca, labels = c('citro', 'sta.lu',
                                        'TDA', 'TDB'))
VR$lote <- factor(VR$lote, labels = paste('l', 1:6, sep = ''))

fruit_quality$lote <- factor(fruit_quality$lote, 
                             labels = paste('l', 1:6, sep = ''))

fruit_quality$planta <- factor(fruit_quality$planta, labels = 
                                 gsub('^(.*)(l)(.*)', '\\1\\3', 
                                      levels(fruit_quality$planta)))

VR$planta <- factor(VR$planta, labels = gsub('^(.*)(l)(.*)', '\\1\\3', 
                                             levels(VR$planta)))

lapply(VR[, -length(VR)], levels)

lapply(fruit_quality[, -length(fruit_quality)], levels)

VR$tot_vis <- VR$tv * 12 * 6  * 2
  
df <- 
  full_join(VR, fruit_quality, 
            by = c('finca', 'variedad', 'lote', 'planta'))

df |> 
  ggplot(aes(tot_vis, diamF, color = variedad)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  facet_wrap(~variedad)

summary(df)

apply(df, 2, function(x) sum(is.na(x)))

df <- na.omit(df)

df$cat_poll <- 
  ifelse(df$tot_vis < 7, 'a 0-7', 
       ifelse(df$tot_vis > 7 & df$tot_vis < 12, 'b 7-12',
              ifelse(df$tot_vis > 12 & df$tot_vis < 18, 'c 12-18',
                     ifelse(df$tot_vis > 18 & df$tot_vis < 27, 'd 18-27', 
                            ifelse(df$tot_vis > 27 & df$tot_vis < 33, 'e 27-33',
                                   ifelse(df$tot_vis > 33 & df$tot_vis < 38, 'f 33-38', 'g >38'))))))

df %$% aggregate(diamF ~ cat_poll + variedad, FUN = mean)
df %$% aggregate(diamF ~ cat_poll + variedad, FUN = length)

df$cat_poll <- as.factor(df$cat_poll)

df$cat_poll <- factor(df$cat_poll, labels = gsub('^([a-z])(.)(.*)$', 
                                                 '\\3', levels(df$cat_poll)))

df |> 
  ggplot(aes(cat_poll, diamF)) +
  stat_summary(fun = 'mean', geom = 'point', 
               position = position_dodge(width = 0.3)) +
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar',
               position = position_dodge(width = 0.3), 
               width = 0.1) +
  stat_summary(aes(group = 1), fun = 'mean', geom = 'line', 
               position = position_dodge(width = 0.3), linetype = 3) +
  #facet_wrap(~variedad, scales = 'free_y') +
  theme_bw() +
  labs(x = 'Range of honeybee visits', y = 'Fruit diameter (mm)') +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(family = 'Times New Roman', angle = 45,
                                   hjust = 1, vjust = 1))
quantile(df[df$cat_poll == '>38',]$tot_vis)



# ====== bayes models =====

paquetes <- c("ggplot2","fitdistrplus","brms",
              "future","GGally",  "ggeffects",
              "performance", "cowplot",
              "ggdist", "arm", "DHARMa","qqplotr",
              "reshape2", "bayesplot", "parallel", 
              "doBy","ggbreak", "ggridges", "ggdist",
              "pROC", 'tidyverse', 'magrittr', 'readxl', 
              'magrittr', 'tidyverse', 'lattice', 'qqplotr', 'DHARMa',
              'ggfortify') 


for (i in 1:length(paquetes)) require(paquetes[[i]], character.only = T) 


df |> 
  ggplot(aes(diamF)) +
  geom_histogram(aes(y = after_stat(density)), 
                 color = 'black', fill = 'white') +
  geom_density() +
  facet_wrap(~variedad) +
  theme_bw()

df$tot_vis1 <- as.vector(scale(df$tot_vis, center = T, scale = T))

modelo <- bf(diamF ~ tot_vis1:variedad + (1|lote/finca))
modelo <- bf(diamF ~ tot_vis:variedad + (1|lote/finca))

get_prior(modelo, data = df, fammily = gaussian)

#plot(density(rnorm(1000, mean = 15, sd = 1)))

prior_modelo <- c(set_prior('normal(0, 2)', class = 'b'),
                  set_prior('normal(14, 2)', class = 'Intercept'),
                  set_prior('normal(0, 1)', class = 'sd', coef = 'Intercept', 
                            group = 'lote:finca'),
                  set_prior('normal(0, 1)', class = 'sd', coef = 'Intercept', 
                            group = 'lote'))

mod_fru <- brm(formula = modelo, data = df, family = gaussian,
               prior = prior_modelo,
               warmup = 1500, future = T, chains = 3, iter = 10000, thin = 3)

summary(mod_fru)

mcmc_dens_overlay(mod_fru, regex_pars = c("^b"))

mcmc_trace(mod_fru, regex_pars = '^b')

mcmc_acf(mod_fru, regex_pars = '^b')

post_fru <- as_draws_df(mod_fru)[, 2:3]

post_fru <- gather(post_fru)

post_fru$key <- as.factor(post_fru$key)

post_fru$key <- factor(post_fru$key, labels = c('Primadonna',
                                                'San Joaquin'))

post_fru |> 
  filter(post_fru$key != 'Intercept') |> 
  ggplot() +
  geom_density(aes(value)) +
  facet_wrap(~key, scales = 'free_x') +
  labs(x = expression(beta[1]), y = "Density")
  
pp_image <- pp_check(mod_fru, ndraws = 1000)

pp_image <- pp_image$data

pp_image <- pp_image[pp_image$rep_label != 'italic(y)',]

ggplot() +
  geom_density(data = pp_image, 
               mapping = aes(x = value, color = rep_label),
               alpha = 0.2, size = 0.1) +
  geom_density(data = df, aes(diamF), color = 'red') +
  theme_bw() +
  scale_color_manual(values = rep('black', nrow(pp_image))) +
  labs(x = 'Fruir size (mm)') +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman'))


post_pred_image <- predict(mod_fru, ndraws = 1000, summary = F)

q_resP1M1 <- createDHARMa(simulatedResponse = t(post_pred_image),
                          observedResponse = na.omit(df$diamF),
                          fittedPredictedResponse = apply(post_pred_image, 2, median),
                          integerResponse = T)

res_m1p1 <- qnorm(residuals(q_resP1M1))

res_m1p1 <- cbind(res = res_m1p1,
                  x = na.omit(df),
                  fitted = fitted(mod_fru, ndraws = 1000)[, 1],
                  pareto = loo(mod_fru, pointwise = T)$diagnostics$pareto_k)

apply(df, 2, function(x) sum(is.na(x)))

res_m1p1 |> 
  ggplot(aes(fitted, res)) +
  geom_point(color = 'cyan4', alpha = 0.4) +
  geom_hline(yintercept = 0, color = 'cyan4', linetype = 3) +
  labs(x = 'Valores ajustados',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 11),
        text = element_text(size = 10, family = 'Times New Roman'))


res_m1p1 |> 
  ggplot(aes(sample = res)) +
  stat_qq_line(color = 'cyan4') +
  stat_qq_band(alpha = 0.4, fill = 'cyan4') +
  stat_qq_point(color = 'tan1') +
  labs(x = 'Cuantiles te?ricos',
       y = 'Cuantiles emp?ricos') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 11),
        text = element_text(size = 10, family = 'Times New Roman'))


res_m1p1 |> 
  ggplot(aes(variedad, res)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color = 'cyan4', linetype = 3) +
  labs(x = 'Valores ajustados',
       y = 'Residuales') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        title = element_text(size = 11),
        text = element_text(size = 10, family = 'Times New Roman'))


# ===== modelo probabilidad entre intervalos =====

df
modelo <- bf(diamF ~ cat_poll:variedad + (1|lote/finca))

get_prior(modelo, data = df, fammily = gaussian)

#plot(density(rnorm(1000, mean = 15, sd = 1)))

prior_modelo <- c(set_prior('normal(14, 2)', class = 'b'),
                  set_prior('normal(14, 2)', class = 'Intercept'),
                  set_prior('normal(0, 1)', class = 'sd', coef = 'Intercept', 
                            group = 'lote:finca'),
                  set_prior('normal(0, 1)', class = 'sd', coef = 'Intercept', 
                            group = 'lote'))

mod_fru <- brm(formula = modelo, data = df, family = gaussian,
               prior = prior_modelo,
               warmup = 1500, future = T, chains = 3, iter = 15000, thin = 3)

summary(mod_fru)
saveRDS(mod_fru, 'mod_catPOLL.rds')
mod_fru <- readRDS('mod_catPOLL.rds')


post_fru <- as_draws_df(mod_fru)[,2:15]

post_fru <- gather(post_fru)

post_fru$variedad <- gsub('^(.*)(:)(variedad)(.*)$', '\\4', post_fru$key)
post_fru$cat_fru <- gsub('(b_cat_poll)(.*)(:)(.*)', '\\2', post_fru$key)

post_fru$cat_fru <- as.factor(post_fru$cat_fru)
post_fru$cat_fru <- factor(post_fru$cat_fru, labels = c(">38", "0-7", "12-18", 
                                      "18-27", "27-33", "33-38", "7-12" ))

post_fru$cat_fru <- factor(post_fru$cat_fru, 
                           levels = c("0-7", "7-12", "12-18", 
                                      "18-27", "27-33", 
                                      "33-38", ">38"))

levels(post_fru$cat_fru)

post_fru$variedad <- as.factor(post_fru$variedad)

levels(post_fru$variedad)

post_fru$variedad <- factor(post_fru$variedad, labels = c('Primadonna',
                                                'San Joaquin'))

post_fru |> 
  ggplot() +
  geom_density(aes(value, color = cat_fru, 
                   fill = cat_fru), alpha = 0.2) +
  facet_wrap(~variedad, scales = 'free_x') +
  labs(x = expression(beta[1]), y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank())

pp_image <- pp_check(mod_fru, ndraws = 1000)

pp_image <- pp_image$data

pp_image <- pp_image[pp_image$rep_label != 'italic(y)',]

ggplot() +
  geom_density(data = pp_image, 
               mapping = aes(x = value, color = rep_label),
               alpha = 0.2, size = 0.1) +
  geom_density(data = df, aes(diamF), color = 'red') +
  theme_bw() +
  scale_color_manual(values = rep('black', nrow(pp_image))) +
  labs(x = 'Fruir size (mm)') +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman'))

library(ggridges)

post_fru |> 
  filter(variedad == 'Primadonna') |> 
  ggplot(aes(value, cat_fru, fill = factor(after_stat(quantile)))) +
  #geom_density_ridges_gradient() +
  #facet_wrap(~variedad, scales = 'free_x') +
  #scale_fill_viridis_c(name = 'Fruis size') +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = F,
    quantiles = c(0.025, 0.5, 0.975), quantile_lines = TRUE,
    alpha = 0.6
  ) +
  scale_fill_manual(values = c('red', 'gray', 'gray', 'red')) +
  labs(x = expression(mu['Fruit diameter']), y = 'Total honeybee visits') +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        text = element_text(family = 'Times New Roman'))
  

post_fru |> 
  group_by(variedad, cat_fru) |> 
  transmute(CI1 = quantile(value, probs = 0.025),
            CI2 = quantile(value, probs = 0.975),
            mu = median(value)) |> 
  unique() |> 
  ggplot() +
    geom_errorbar(aes(ymin = CI1, 
                    ymax = CI2, 
                    x = cat_fru), width = 0.2) +
  stat_summary(aes(group = 1, y = mu, x = cat_fru), fun = 'mean', 
               geom = 'line', position = 
                 position_dodge(), linetype = 3) +
  geom_point(aes(cat_fru, mu), color = 'red') +
  facet_wrap(~variedad) +
  labs(y = expression(mu['Fruit diameter']), x = 'Total honeybee visits') +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        text = element_text(family = 'Times New Roman'))

  

prim <- 
  split(post_fru[post_fru$variedad == "Primadonna",], 
      list(post_fru$variedad, post_fru$cat_fru))

prim <- prim[unlist(lapply(prim, FUN = nrow), use.names = F) != 0]

perc_probs <- function(x1, x2, menor_que = T) {
  if (menor_que == T) {
    return((sum(x1 <= mean(x2)) / length(x1)))
  } else {
    return(sum(x1 >= mean(x2)) / length(x1))
  }
}

prim <- tibble(vis = levels(post_fru$cat_fru),
       probs = sapply(prim, FUN = 
                        function(x) {
                          perc_probs(x$value,
                                     prim$`Primadonna.0-7`$value, menor_que = T)
                        }, USE.NAMES = F)) 

prim$vis <- as.factor(prim$vis)
prim$vis <- factor(prim$vis, levels = levels(post_fru$cat_fru))

prim |> 
  ggplot(aes(vis, probs)) +
  geom_line(aes(group = 1), linetype = 3) +
  geom_point() +
  labs(y = 'Probability of fruit decrease', x = 'Total honeybee visits') +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        text = element_text(family = 'Times New Roman'))


sanJ <- 
  split(post_fru[post_fru$variedad == "San Joaquin",], 
        list(post_fru$variedad, post_fru$cat_fru))

sanJ <- sanJ[unlist(lapply(sanJ, FUN = nrow), use.names = F) != 0]


sanJ <- tibble(vis = levels(post_fru$cat_fru),
               probs = sapply(sanJ, FUN = 
                                function(x) {
                                  perc_probs(x$value,
                                             sanJ$`Primadonna.0-7`$value, menor_que = T)
                                }, USE.NAMES = F)) 

sanJ$vis <- as.factor(sanJ$vis)
sanJ$vis <- factor(sanJ$vis, levels = levels(post_fru$cat_fru))

sanJ |> 
  ggplot(aes(vis, probs)) +
  geom_line(aes(group = 1), linetype = 3) +
  geom_point() +
  labs(y = 'Probability of fruit decrease', x = 'Total honeybee visits') +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        text = element_text(family = 'Times New Roman'))



