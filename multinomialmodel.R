library(readr)
library(tidyverse)
library(stats)

# Read csv
# Database avalaible at:
# https://www.pns.icict.fiocruz.br/wp-content/uploads/2023/11/pns2019.zip
pnsdata <- read_delim("PNS_2019\\PNS_2019.csv", delim=';')  

#Selecting Columns
pnsdata <- pnsdata[complete.cases(pnsdata[, c('C008', 'I00102','J012')]), ]

# Define column 'sinistro'
pnsdata$sinistro <- ifelse(pnsdata$J012 >= 3, 1, 0)

# Define column 'cobertura_plano'
pnsdata$cobertura_plano <- (-1) * (pnsdata$I00102 - 2)

# Creating age bins
bins <- c(0, 18, 23, 28, 33, 38, 43, 48, 53, 59, Inf)
labels <- c('0-18', '19-23', '24-28', '29-33', '34-38', '39-43', '44-48', '49-53', '54-59', '59+')

pnsdata$age_group <- cut(pnsdata$C008, breaks = bins, labels = labels, right = FALSE)

########################
### CANONICAL MODEL
# Probit model for 'cobertura_plano' (coverage)
modelocobertura <- glm(cobertura_plano ~ age_group, data = pnsdata, family = binomial(link = "probit"))

# Probit model for 'sinistro' (risk of loss)
modelosinistro <- glm(sinistro ~ age_group, data = pnsdata, family = binomial(link = "probit"))

# Extracting residuals
residuals1 <- modelocobertura$residuals  
residuals2 <- modelosinistro$residuals

numerator <- sum(residuals1 * residuals2)^2
denominator <- sum((residuals1^2) * (residuals2^2))
W_statistic <- numerator / denominator

W_statistic

########################
########################
### MULTINOMIAL MODEL (KIM ET. AL 2009)
## ORDERED PROBIT FOR 'cobertura' (coverage)
library(MASS)

vetor_temp <- ifelse(pnsdata$I010010 == 3, 2, pnsdata$cobertura_plano)
pnsdata$cobertura_plano_ordered <- replace(vetor_temp, is.na(vetor_temp), 0)
pnsdata$cobertura_plano_ordered


probitordenado <- polr(as.ordered(cobertura_plano_ordered)~age_group, data = pnsdata,method='probit')
summary(probitordenado)

#NOW CALCULATE RESIDUALS
X_matrix <- model.matrix(modelocobertura)
XtimesB <- X_matrix %*% modelocobertura$coefficients#[-1]
xis <- XtimesB - probitordenado$zeta[1]

#Define parameters
mu <- 0
sigma <- 1

residuoespecial1 <- (dnorm(xis, mean = mu, sd = sigma)*(pnsdata$cobertura_plano-(pnorm(xis, mean = mu, sd = sigma)))) / ((pnorm(xis, mean = mu, sd = sigma) * (1 - pnorm(xis, mean = mu, sd = sigma))))
xis <- XtimesB - probitordenado$zeta[2]

pnsdata$cobertura_plano2 <- pnsdata$cobertura_plano_ordered-1
pnsdata$cobertura_plano2 <- ifelse(pnsdata$cobertura_plano2 == -1, 0, pnsdata$cobertura_plano2)
pnsdata$cobertura_plano2
residuoespecial2 <- (dnorm(xis, mean = mu, sd = sigma)*(pnsdata$cobertura_plano2-(pnorm(xis, mean = mu, sd = sigma)))) / ((pnorm(xis, mean = mu, sd = sigma) * (1 - pnorm(xis, mean = mu, sd = sigma))))


########################
## Sinistro/Risk of Loss (Poisson QMLE or Negative Binomial)
## Poisson QMLE
model_poisson <- glm(J012 ~ age_group +residuoespecial1 + residuoespecial2, 
                     data = pnsdata, 
                     family = quasipoisson,maxit=1000)

summary(model_poisson)

## Negative Binomial

nb<-glm.nb(formula = J012 ~ c(age_group)+residuoespecial1+residuoespecial2, data = bla, init.theta = 40, link = "log",maxit=3000)

summary(nb)







