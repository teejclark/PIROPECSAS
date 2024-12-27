#### Bootstrap Analysis

# calculate conversion factor
hist(densities$conv_fact)
mean(densities$conv_fact)
median(densities$conv_fact)

(conversion.factor <- median(densities$conv_fact))
(conversion.factor.ci <- quantile(densities$conv_fact, c(0.05, 0.95)))


reps <- 1000
#conversion.factor.i <- rnorm(reps, mean = 0.9601281, sd = 0.4267296)
growth.rates.i <- rep(NA,reps)
  

# for loop
for (i in 1:reps){
  
conversion.factor.i <- rnorm(1, mean = 0.9601281, sd = 0.4267296)

# PIROP data
PIROP_estimates <- left_join(PIROP_abund, PIROP_yearly_effort, by = join_by(year == Region.Label)) %>% 
  mutate(mean_unadj = (tot_birds/Effort),
         mean = (tot_birds/Effort)/conversion.factor.i) %>%
  # If using gam to predict equivalent ECSAS linear density for PIROP lin dens uncomment next line
  # mean = predict(m1, newdata = data.frame(PIROP_lindens = tot_birds/Effort))) %>%
  select(year, mean, tot_birds, Effort, mean_unadj)

PIROP_estimates <- PIROP_estimates %>%
  mutate(cv = head(gannet_PIROP_300$Abundance_CV, - 1),
         se = mean * cv,
         LCI = lognormal.ci(mean, cv)$lcl,
         UCI = lognormal.ci(mean, cv)$ucl
  )

overall_estimates <- data.frame(year = seq(1969,2021,1),
                                PIROP_mean = c(PIROP_estimates$mean, rep(NA,length(seq(1984,2021)))),
                                PIROP_LCI = c(PIROP_estimates$LCI, rep(NA,length(seq(1984,2021)))),
                                PIROP_UCI = c(PIROP_estimates$UCI, rep(NA,length(seq(1984,2021)))),
                                ECSAS_mean = c(rep(NA,length(seq(1969,2005))), ECSAS_estimates$mean),
                                ECSAS_LCI = c(rep(NA,length(seq(1969,2005))), ECSAS_estimates$LCI),
                                ECSAS_UCI = c(rep(NA,length(seq(1969,2005))), ECSAS_estimates$UCI))


# calculate growth rate
fit <- data.frame(year = overall_estimates$year,
                  mean = c(overall_estimates$PIROP_mean[1:37],
                           overall_estimates$ECSAS_mean[38:53]),
                  Program = c(rep("PIROP", 37), rep("ECSAS", 16)))
fitlm <- lm(fit$mean ~ fit$year)
predlm <- predict.lm(fitlm, newdata = data.frame(fit$year),
                     type = "response")

(predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]
mean((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)])
quantile((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)], c(0.05,0.95))

# Same as above - just simpler to understand
numerator <- rep(coef(fitlm)[2], length(predlm) - 1)
denominator <- head(predlm, length(predlm) - 1)
growth.rates <- numerator/denominator
growth.rates

(med <- round(median(growth.rates), 4))
mean(growth.rates)
(quant <- round(quantile(growth.rates, c(0.05,0.95)), 4))

growth.rates.i[i] <- round(median(growth.rates), 4)

}

#######################################################################################################################################

# Randomization Analyses

densities$conv_fact

conversion.factor <- median(densities$conv_fact)

# function to calculate median

calculate_median <- function(vector) {
  sorted_vector <- sort(vector)
  n <- length(sorted_vector)
  if (n %% 2 == 0) {
    return((sorted_vector[n/2] + sorted_vector[n/2 + 1]) / 2)
  } else {
    return(sorted_vector[(n + 1) / 2])
  }
}

# Function to calculate medians for different combinations of lengths
calculate_median_combinations <- function(vector) {
  subsets <- list()
  n <- length(vector)
  for (i in 1:n) {
    for (j in 1:(n - i + 1)) {
      subset <- vector[j:(j + i - 1)]
      subsets <- c(subsets, list(subset))
    }
  }
  
  results <- data.frame(length = integer(), median = numeric())
  for (k in 1:length(subsets)) {
    subset <- subsets[[k]]
    med <- calculate_median(subset)
    len <- length(subset)
    results <- rbind(results, data.frame(length = len, median = med))
  }
  
  return(results)
}

# test
vector <- c(7, 3, 1, 4, 6, 2, 5)
medians <- calculate_median_combinations(vector)
print(medians)


random.medians <- calculate_median_combinations(densities$conv_fact)


reps <- nrow(random.medians)
conversion.factor.i <- random.medians$median
growth.rates.i <- rep(NA,reps)


# for loop
for (i in 1:reps){
  
  # PIROP data
  PIROP_estimates <- left_join(PIROP_abund, PIROP_yearly_effort, by = join_by(year == Region.Label)) %>% 
    mutate(mean_unadj = (tot_birds/Effort),
           mean = (tot_birds/Effort)/conversion.factor.i[i]) %>%
    # If using gam to predict equivalent ECSAS linear density for PIROP lin dens uncomment next line
    # mean = predict(m1, newdata = data.frame(PIROP_lindens = tot_birds/Effort))) %>%
    select(year, mean, tot_birds, Effort, mean_unadj)
  
  PIROP_estimates <- PIROP_estimates %>%
    mutate(cv = head(gannet_PIROP_300$Abundance_CV, - 1),
           se = mean * cv,
           LCI = lognormal.ci(mean, cv)$lcl,
           UCI = lognormal.ci(mean, cv)$ucl
    )
  
  overall_estimates <- data.frame(year = seq(1969,2021,1),
                                  PIROP_mean = c(PIROP_estimates$mean, rep(NA,length(seq(1984,2021)))),
                                  PIROP_LCI = c(PIROP_estimates$LCI, rep(NA,length(seq(1984,2021)))),
                                  PIROP_UCI = c(PIROP_estimates$UCI, rep(NA,length(seq(1984,2021)))),
                                  ECSAS_mean = c(rep(NA,length(seq(1969,2005))), ECSAS_estimates$mean),
                                  ECSAS_LCI = c(rep(NA,length(seq(1969,2005))), ECSAS_estimates$LCI),
                                  ECSAS_UCI = c(rep(NA,length(seq(1969,2005))), ECSAS_estimates$UCI))
  
  
  # calculate growth rate
  fit <- data.frame(year = overall_estimates$year,
                    mean = c(overall_estimates$PIROP_mean[1:37],
                             overall_estimates$ECSAS_mean[38:53]),
                    Program = c(rep("PIROP", 37), rep("ECSAS", 16)))
  fitlm <- lm(fit$mean ~ fit$year)
  predlm <- predict.lm(fitlm, newdata = data.frame(fit$year),
                       type = "response")
  
  (predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]
  mean((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)])
  quantile((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)], c(0.05,0.95))
  
  # Same as above - just simpler to understand
  numerator <- rep(coef(fitlm)[2], length(predlm) - 1)
  denominator <- head(predlm, length(predlm) - 1)
  growth.rates <- numerator/denominator
  growth.rates
  
  (med <- round(median(growth.rates), 4))
  mean(growth.rates)
  (quant <- round(quantile(growth.rates, c(0.05,0.95)), 4))
  
  growth.rates.i[i] <- round(median(growth.rates), 4)
  
}

random.medians$growth.rate <- growth.rates.i

# plot relationship between # of values and growth rate
plot(random.medians$length, random.medians$growth.rate)

ggplot(random.medians)+
  geom_point(aes(length, growth.rate), size = 3)+
  geom_hline(yintercept = random.medians$growth.rate[reps], linetype = "dashed")






