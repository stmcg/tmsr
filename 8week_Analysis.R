library('MASS')
library(readr) 
library(magic)
library(metafor)


get_cov <- function(mean_val, sd_val, n, rho = 0.5, n_MonteCarlo = 1e6, seed = 1234){
  set.seed(seed)
  nTime <- length(mean_val)
  yi <- (mean_val[2:nTime] - mean_val[1:(nTime-1)]) / mean_val[1:(nTime-1)] 
  
  Sigma <- matrix(NA, nrow = nTime, ncol = nTime)
  for (i in 1:nTime){
    for (j in 1:nTime){
      if (i == j){
        Sigma[i, j] <- sd_val[i]^2 / n[i]
      } else {
        Sigma[i, j] <- rho * sd_val[i] * sd_val[j] / sqrt(n[i] * n[j])
      }
    }
  }
  monte_temp <- mvrnorm(n_MonteCarlo, mu = mean_val, Sigma = Sigma)
  monte_foldchanges <- (monte_temp[,2:ncol(monte_temp)] - monte_temp[,1:(ncol(monte_temp)-1)]) / 
    monte_temp[,1:(ncol(monte_temp)-1)]
  monte_foldchanges <- as.matrix(monte_foldchanges)
  
  # Checking for potential issues with numerical stability
  for (i in ncol(monte_foldchanges)){
    large_ind <- abs(monte_foldchanges[,i] - yi[i]) > 8 * sd(monte_foldchanges[,i])
    monte_foldchanges[large_ind,i] <- NA
    temp <- 100 * sum(large_ind) / n_MonteCarlo
    print(paste0('Removed ', sum(large_ind), ' observations (', temp, '%) extremely large values in Monte Carlo integration'))
  }
  monte_foldchanges <- as.matrix(monte_foldchanges[complete.cases(monte_foldchanges),])
  V <- cov(monte_foldchanges)
  return(list(yi = yi, V = V))
}


# Including Ferrian
crp <- readxl:: read_xlsx("./TMSR_FinalData_20210803.xlsx", sheet = "CRP")

crp$N = as.numeric(crp$N)
crp$Weeks = as.numeric(crp$Weeks)

crp <- crp[crp$Weeks %in% c(0,8) & ! (crp$Author %in% c('Franscisco', 'Moraes')),] #reducing data set to correct studies
crp <- crp[complete.cases(crp), ]

all_res <- vector(mode = "list", length = length(unique(crp$Author)))

i <- 1

for (author in unique(crp$Author)){
  study_crp <- crp[crp$Author == author,]
  all_res[[i]] <- get_cov(mean_val = study_crp$Mean, 
                          sd_val = study_crp$SD, 
                          n = study_crp$N, 
                          rho = 0.5)
  i <- i + 1
}

length(all_res)

yi <- c(unlist(sapply(all_res, '[[', 'yi')))
V <- c(unlist(sapply(all_res, '[[', 'V')))


crp_temp <- crp[crp$Weeks != 0,]

metacrp <- rma.uni(yi = yi, vi = V, data = crp_temp)
summary(metacrp)




# Excluding Ferrian
crp <- readxl:: read_xlsx("./TMSR_FinalData_20210803.xlsx", sheet = "CRP")

crp$N = as.numeric(crp$N)
crp$Weeks = as.numeric(crp$Weeks)

crp <- crp[crp$Weeks %in% c(0,8) & ! (crp$Author %in% c('Franscisco', 'Moraes', 'Ferrian- SR', 'Ferrian-FR')),] #reducing data set to correct studies
crp <- crp[complete.cases(crp), ]

all_res <- vector(mode = "list", length = length(unique(crp$Author)))

i <- 1

for (author in unique(crp$Author)){
  study_crp <- crp[crp$Author == author,]
  all_res[[i]] <- get_cov(mean_val = study_crp$Mean, 
                          sd_val = study_crp$SD, 
                          n = study_crp$N, 
                          rho = 0.5)
  i <- i + 1
}

length(all_res)

yi <- c(unlist(sapply(all_res, '[[', 'yi')))
V <- c(unlist(sapply(all_res, '[[', 'V')))


crp_temp <- crp[crp$Weeks != 0,]

metacrp <- rma.uni(yi = yi, vi = V, data = crp_temp)
summary(metacrp)




il6 <- readxl:: read_xlsx("./TMSR_FinalData_20210803.xlsx", sheet = "IL6")

il6$N = as.numeric(il6$N)
il6$Weeks = as.numeric(il6$Weeks)

il6 <- il6[il6$Weeks %in% c(0,8) & ! (il6$Author %in% c('Djoba Siawaya')),] #reducing data set to correct studies
il6 <- il6[complete.cases(il6), ]

all_res3 <- vector(mode = "list", length = length(unique(il6$Author)))

i <- 1

for (author in unique(il6$Author)){
  study_il6 <- il6[il6$Author == author,]
  all_res3[[i]] <- get_cov(mean_val = study_il6$Mean, 
                           sd_val = study_il6$SD, 
                           n = study_il6$N, 
                           rho = 0.5)
  i <- i + 1
}

length(all_res3)

yi3 <- c(unlist(sapply(all_res3, '[[', 'yi')))
V3 <- c(unlist(sapply(all_res3, '[[', 'V')))



il6_temp <- il6[il6$Weeks != 0,]

metail6 <- rma.uni(yi = yi3, vi = V3, data = il6_temp)

summary(metail6)



ip10 <- readxl:: read_xlsx("./TMSR_FinalData_20210803.xlsx", sheet = "IP10")

ip10$Weeks = as.numeric(ip10$Weeks)

ip10$n = as.numeric(ip10$n)
ip10$sd_value = as.numeric(ip10$sd_value)

ip10 <- ip10[ip10$Weeks %in% c(0,8) & ! (ip10$Author %in% c('Djoba', 'Hong- FR', 'Hong-SR', 'Francisco')),] #reducing data set to correct studies
ip10 <- ip10[complete.cases(ip10), ]

all_res2 <- vector(mode = "list", length = length(unique(ip10$Author)))

i <- 1

for (author in unique(ip10$Author)) {
  study_ip10 <- ip10[ip10$Author == author,]
  all_res2[[i]] <- get_cov(mean_val = study_ip10$mean_value, 
                           sd_val = study_ip10$sd_value, 
                           n = study_ip10$n, 
                           rho = 0.5)
  i <- i+ 1
  
}

length(all_res2)

yi2 <- c(unlist(sapply(all_res2, '[[', 'yi')))
V2 <- c(unlist(sapply(all_res2, '[[', 'V')))

ip10_temp <- ip10[ip10$Weeks != 0,]

metaip10 <- rma.uni(yi = yi2, vi = V2, data = ip10_temp)

summary(metaip10)




tnf <- readxl:: read_xlsx("./TMSR_FinalData_20210803.xlsx", sheet = "TNFa")

tnf$Weeks = as.numeric(tnf$Weeks)

tnf$n = as.numeric(tnf$n)

tnf <- tnf[tnf$Weeks %in% c(0,8) & ! (tnf$Author %in% c('Djoba Siawaya')),] #reducing data set to correct studies
tnf <- tnf[complete.cases(tnf), ]

all_res4 <- vector(mode = "list", length = length(unique(tnf$Author)))

i <- 1

for (author in unique(tnf$Author)) {
  study_tnf <- tnf[tnf$Author == author,]
  all_res4[[i]] <- get_cov(mean_val = study_tnf$mean_value, 
                           sd_val = study_tnf$sd_value, 
                           n = study_tnf$n, 
                           rho = 0.5)
  i <- i+ 1
  
}

length(all_res4)

yi4 <- c(unlist(sapply(all_res4, '[[', 'yi')))
V4 <- c(unlist(sapply(all_res4, '[[', 'V')))


tnf_temp <- tnf[tnf$Weeks != 0,]

metatnf <- rma.uni(yi = yi4, vi = V4, data = ip10_temp)

summary(metatnf)

