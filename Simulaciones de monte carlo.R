library(statsr)
library(ggplot2)
data(tapwater)
str(tapwater)
summary(tapwater)
var(tapwater$tthm)
library(statsr)
library(ggthemes)
library(tidyverse)
library(rjags)
install.packages("rtools")
install.packages("statsr")
library(statsr)
install.packages("BAS")
library(BAS)
library(ggplot2)
library(dplyr)
install.packages("BayesFactor")
library(BayesFactor)
library(knitr)
library(rjags)
install.packages("coda")
library(coda)
install.packages("latex2exp")
library(latex2exp)
install.packages("foreign")
library(foreign)
install.packages("BHH2")
library(BHH2)
install.packages("scales")
library(scales)
install.packages("logspline")
library(logspline)
install.packages("cowplot")
library(cowplot)
library(ggthemes)
remotes::install_github("statswithr/statsr@BayesFactor")



# prior hyperparameters
m_0 = 35; n_0 = 25; s2_0 = 156.25; v_0 = n_0 - 1
# sample summaries
Y = tapwater$tthm
ybar = mean(Y)
s2 = var(Y)
n = length(Y)
# posterior hyperparamters
n_n = n_0 + n
m_n = (n*ybar + n_0*m_0)/n_n
v_n = v_0 + n
s2_n = ((n-1)*s2 + v_0*s2_0 + n_0*n*(m_0 - ybar)^2/n_n)/v_n

bayes_inference(tthm, data=tapwater, prior="NG",
                hypothesis_prior =c(mu_0 = m_0, n_0=n_0, s_0 = sqrt(s2_0), v_0 = v_0),
                statistic="mean", type= c("ci"), method="theoretical",
                cred_level = 0.95, show_res=TRUE, show_summ=TRUE, show_plot=FALSE)
?bayes_inference

# Previa de hiper-parámetros
m_0 = 35; n_0 = 25; s2_0 = 156.25; v_0 = n_0 - 1
# Datos
data(tapwater); Y = tapwater$tthm
ybar = mean(Y); s2 = var(Y); n = length(Y)
# Hiper-parámetros posteriores
n_n = n_0 + n
m_n = (n*ybar + n_0*m_0)/n_n
v_n = v_0 + n
s2_n = ((n-1)*s2 + v_0*s2_0 + n_0*n*(m_0 - ybar)^2/n_n)/v_n

#semilla para que la simulación sea consisstente a lo largo del tiempo
set.seed(42)

#Generación de la simulación Monte Carlo
phi = rgamma(1000, shape = v_n/2, rate=s2_n*v_n/2)
?sort()
?geom_density()
?geom_line()
df = data.frame(phi = sort(phi))
df = mutate(df,
            density = dgamma(phi,
                             shape = v_n/2,
                             rate=s2_n*v_n/2))

#Plot de la simulación
ggplot(data=df, aes(x=phi)) +
  geom_histogram(aes(x=phi, y=..density..), bins = 50) +
  geom_density(aes(phi, ..density..), color="black") +
  geom_line(aes(x=phi, y=density), color="orange") +
  xlab(expression(phi)) + theme_tufte()


#Parámetro aproximado ye intervalo de credivilidad
mean(phi)
quantile(phi, c(0.025, 0.975))

#Usar el parámetro phi para obtener sigma
sigma= 1/sqrt(phi)
tapwater<- cbind(tapwater, sigma=sigma)
  
mean(sigma)
quantile(sigma, c(0.025, 0.975))

#Distribución de sigma
df2=data_frame(sigma= sort(sigma))

ggplot(data=df, aes(x=sigma)) +
  geom_histogram(aes(x=sigma, y=..density..), bins = 50) +
  geom_density(aes(sigma, ..density..), color="black") +
  xlab(expression(sigma)) + theme_tufte()

#Simulación para obtener parámetros previos 
m_0 = (60+10)/2; s2_0 = ((60-10)/4)^2;
n_0 = 2; v_0 = n_0 - 1
set.seed(1234)
S = 10000
phi = rgamma(S, v_0/2, s2_0*v_0/2)
sigma = 1/sqrt(phi)
mu = rnorm(S, mean=m_0, sd=sigma/(sqrt(n_0)))
Y = rnorm(S, mu, sigma)
quantile(Y, c(0.025,0.975))

#Nueva n_0
m_0 = (60+10)/2; s2_0 = ((60-10)/4)^2;
n_0 = 25; v_0 = n_0 - 1
set.seed(1234)
phi = rgamma(10000, v_0/2, s2_0*v_0/2)
sigma = 1/sqrt(phi)
mu = rnorm(10000, mean=m_0, sd=sigma/(sqrt(n_0)))
y = rnorm(10000, mu, sigma)
quantile(y, c(0.025,0.975))

#Obtenemos el estimado de Y por debajo de cero
sum(y < 0)/length(y)

#Uso de monte carlo para predicción 
set.seed(1234)
phi = rgamma(10000, v_n/2, s2_n*v_n/2)
sigma = 1/sqrt(phi)
post_mu = rnorm(10000, mean=m_n, sd=sigma/(sqrt(n_n)))
pred_y = rnorm(10000,post_mu, sigma)
quantile(pred_y, c(.025, .975))

#Ahora lo que hacemos es generar la predicción de que pase cierto límite
sum(pred_y > 80)/length(pred_y)


phi = rgamma(10000, (n-1)/2, s2*(n-1)/2)
sigma = 1/sqrt(phi)
post_mu = rnorm(10000, mean=ybar, sd=sigma/(sqrt(n)))
pred_y = rnorm(10000,post_mu, sigma)
quantile(pred_y, c(.025, .975))

# Markov Chain Monte Carlo (MCMC) pseudo código (no sirve realmente)
sigma[1] = 1; n_0[1]=1; mu[1]=m_0
#draw from full conditional distributions
for (i in 2:S) {
  mu[i] = p_mu(sigma[i-1], n_0[i-1], m_0, r, data)
  sigma[i] = p_sigma(mu[i], n_0[i-1], m_0, r, data)
  n_0[i] = p_n_0(mu[i], sigma[i], m_0, r, data)
}

#Ahora sí el bueno pero que hace todo tras bambalinas
bayes_inference(y=tthm, data=tapwater, statistic="mean",
                mu_0 = 35, rscale=1, prior="JZS",
                type="ci", method="sim")

#Ejemplo con la distribución de Jeffrey-Zellener-Siow comparando dos medias no independientes
library(statsr)
?bayes_inference
bayes_inference(difference, data=zinc, statistic="mean", type="ht",
                prior="JZS", mu_0=0, method="theo", alt="twosided")

#Comparación de dos medias independientes
library(statsr)
data(nc)
bayes_inference(y=gained, x=mature, data=nc, type= "ht",
                statistic="mean", alternative="twosided", null=0,
                prior="JZS", r=1, method="theo", show_sum= FALSE)


#Otra prueba usando el JZS pero con r=0.05
library(statsr)
data(nc)
out=bayes_inference(y=weight, x=habit, data =nc, type= "ht",
                    statistic = "mean", alternative = "twosided", null = 0,
                    prior="JZS", r=0.05, method="sim", show_summ = F)

#Intervalos de credibilidad
out.ci = bayes_inference(y=weight, x=habit, data=nc, type='ci',
                         statistic='mean', prior='JZS', mu_0=0,
                         r=.5, method='sim', verbose=FALSE)

print(out.ci$summary, digits=2)


library(ggplot2)
out = bayes_inference(y=weight, x=habit, data=nc,type='ht',
                      statistic='mean', alternative='twosided',
                      prior='JZS', null=0, r=.5, method='theo',
                      show_summ=FALSE, show_res=FALSE, show_plot=TRUE)