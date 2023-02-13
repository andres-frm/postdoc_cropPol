library(tidyverse)
library(bbmle)

df <- iris[iris$Species == 'virginica',]

plot(density(df$Sepal.Length))


library(fitdistrplus)

norm <- fitdist(df$Petal.Length, 'norm')
lnorm <- fitdist(df$Petal.Length, 'lnorm')
gamma <- fitdist(df$Petal.Length, 'gamma')
cdfcomp(list(norm, lnorm, gamma))
qqcomp(list(norm, lnorm, gamma))
gofstat(list(norm, lnorm, gamma))

x <- rnorm(50, mean = 3, sd = 2)

# mle for the mean 
LL <- function(mu, sigma) {
  -sum(log(dnorm(df$Petal.Width, mean = mu, sd = sigma)))
}

m_LL <- mle2(minuslogl = LL, 
            start = list(mu = 1, sigma = 1))

summary(m_LL)
confint(m_LL)
tibble(mu = mean(df$Petal.Width), sigma = sd(df$Petal.Width))

# mle for relation among variables

LL <- function(b0, b1, mu, sigma) {
  
  n <- df$Petal.Width - (b0 + df$Sepal.Width * b1)
  
  -sum(log(dnorm(n, mean = mu, sd = sigma)))
}

m_LL <- mle2(minuslogl = LL, 
             start = list(mu = 0, sigma = 1,
                          b0 = 1,
                          b1 = 0.2))

summary(m_LL)
summary(lm(df$Petal.Width ~ df$Sepal.Width))
confint(m_LL)

# mle mu paramether
s <- seq(0.1, 20, by = 0.01)

LL_res <- vector('double', length(s))

for (i in 1:length(s)) {
  LL_res[[i]] <- sum(log(dnorm(df$Sepal.Length, mean = s[[i]], sd = 0.636)))
}

q_ll <- tibble(LL = LL_res,
               shape = s,
               shape2 = s^2)

q_ll <- lm(LL ~ shape + shape2, data = q_ll)

q_ll <- summary(q_ll)

vertex <- -q_ll$coefficients[2,1]/(2 * q_ll$coefficients[3,1])

tibble(LL = LL_res,
       shape = s) |> 
  ggplot(aes(s, LL)) +
  geom_line() +
  geom_vline(xintercept = vertex, linetype = 3) +
  geom_label(label = round(vertex, 2), x = vertex + 1, y = -9000) +
  labs(x = expression(mu), y = "log likelihood")


# mtcars data frame 
GGally::ggcorr(mtcars, label = T)

plot(mtcars$disp, mtcars$hp)

plot(density(mtcars$disp))



LL <- function(mu, sigma) {
  -sum(log(dnorm(mtcars$disp, mean = mu, sd = sigma)))
}

m_mle <- mle2(minuslogl = LL,
              start = list(mu = 200, sigma = 50))
summary(m_mle)

# reg

LL <- function(b0, b1, mu, sigma) {
  y <- mtcars$disp - (b0 + b1 * mtcars$hp)
  
  - sum(log(dnorm(y, mean = mu, sd = sigma)))
  
}

n_mle <- mle2(minuslogl = LL, 
              start = list(b0 = 20,
                           b1 = 0.5,
                           mu = 200,
                           sigma = 100))

summary(n_mle)
summary(lm(disp ~ hp, data = mtcars))

df <- na.omit(dplyr::starwars)

df %$% plot(height, mass)

df <- df[df$mass != min(df$mass), ]


LL <- function(b0, b1, mu, sigma) {
  y <- df$mass - (b0 + b1 * df$height)
  -sum(log(dnorm(y, mean = mu, sd = sigma)))
}

m_mle <- mle2(minuslogl = LL,
              start = list(b0 = 10,
                           b1 = 0.5,
                           mu = 80,
                           sigma = 100))

summary(m_mle)
summary(df %$% lm(mass ~ height))









