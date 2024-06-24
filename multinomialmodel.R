library(readr)
library(tidyverse)
library(stats)

# Ler csv
# Base de dados disponível em:
# https://www.pns.icict.fiocruz.br/wp-content/uploads/2023/11/pns2019.zip
pnsdata <- read_delim("PNS_2019\\PNS_2019.csv", delim=';')  


pnsdata <- pnsdata[complete.cases(pnsdata[, c('C008', 'I00102','J012')]), ]

# Definir coluna 'sinistro'
pnsdata$sinistro <- ifelse(pnsdata$J012 >= 3, 1, 0)

# Definir coluna 'cobertura_plano'
pnsdata$cobertura_plano <- (-1) * (pnsdata$I00102 - 2)


# Criar faixas de idade
bins <- c(0, 18, 23, 28, 33, 38, 43, 48, 53, 59, Inf)
labels <- c('0-18', '19-23', '24-28', '29-33', '34-38', '39-43', '44-48', '49-53', '54-59', '59+')

pnsdata$age_group <- cut(pnsdata$C008, breaks = bins, labels = labels, right = FALSE)

# Modelo probit para 'cobertura_plano'
modelocobertura <- glm(cobertura_plano ~ age_group, data = pnsdata, family = binomial(link = "probit"))

# Modelo probit para 'sinistro'
modelosinistro <- glm(sinistro ~ age_group, data = pnsdata, family = binomial(link = "probit"))


# Extraindo resíduos
residuals1 <- modelocobertura$residuals  
residuals2 <- modelosinistro$residuals


numerator <- sum(residuals1 * residuals2)^2
denominator <- sum((residuals1^2) * (residuals2^2))
W_statistic <- numerator / denominator

W_statistic


## Probit Ordenado Modelo Multinomial
library(MASS)

vetor_temp <- ifelse(pnsdata$I010010 == 3, 2, pnsdata$cobertura_plano)
pnsdata$cobertura_plano_ordered <- replace(vetor_temp, is.na(vetor_temp), 0)
pnsdata$cobertura_plano_ordered


probitordenado <- polr(as.ordered(cobertura_plano_ordered)~age_group, data = pnsdata,method='probit')
summary(probitordenado)

#calculando residuos
X_matrix <- model.matrix(modelocobertura)

XtimesB <- X_matrix %*% modelocobertura$coefficients#[-1]

xis <- XtimesB - probitordenado$zeta[1]

#definindo parametros 
mu <- 0
sigma <- 1


residuoespecial1 <- (dnorm(xis, mean = mu, sd = sigma)*(pnsdata$cobertura_plano-(pnorm(xis, mean = mu, sd = sigma)))) / ((pnorm(xis, mean = mu, sd = sigma) * (1 - pnorm(xis, mean = mu, sd = sigma))))

xis <- XtimesB - probitordenado$zeta[2]

pnsdata$cobertura_plano2 <- pnsdata$cobertura_plano_ordered-1
pnsdata$cobertura_plano2 <- ifelse(pnsdata$cobertura_plano2 == -1, 0, pnsdata$cobertura_plano2)
pnsdata$cobertura_plano2
residuoespecial2 <- (dnorm(xis, mean = mu, sd = sigma)*(pnsdata$cobertura_plano2-(pnorm(xis, mean = mu, sd = sigma)))) / ((pnorm(xis, mean = mu, sd = sigma) * (1 - pnorm(xis, mean = mu, sd = sigma))))


## Poisson QMLE
model_poisson <- glm(J012 ~ age_group +residuoespecial1 + residuoespecial2, 
                     data = pnsdata, 
                     family = quasipoisson,maxit=1000)

summary(model_poisson)

## Negative Binomial

nb<-glm.nb(formula = J012 ~ c(age_group)+residuoespecial1+residuoespecial2, data = bla, init.theta = 40, link = "log",maxit=3000)

summary(nb)







