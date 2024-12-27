#### Calculating trends in Gannets

library(tidyverse)
library(ggplot2)
library(lubridate)
library(broom)

###############################################################################################################################

# load the data from Atlantic Canada
# NOTE = ignore Whitehorse, Machias Seal Island, and green
atlantic <- read.csv("./Data/RData/NOGAColonialSeabirdData_ALLAtlanticRegion.csv",
                     fileEncoding="latin1") 

# extract year, stupid conversions for pre-70s data
atlantic$year <- year(as.Date(atlantic$Date, "%d-%b-%y"))
atlantic$year[atlantic$year > 2020] <- atlantic$year[atlantic$year > 2020] - 100

# plot each time-series by year for each area (6 colonies)
atlantic %>%
filter(Colony.Id != "Machias Seal Island, NB" & Colony.Id != "Whitehorse Island, NB" & Colony.Id != "Green Island 1, NS") %>%
ggplot(aes(year, Colony.size, color = Colony.Id)) +
  geom_point()+
  geom_line(lwd=1.5)

# load data from Quebec
quebec <- read.csv("./Data/Rdata/quebec gannet censuses.csv")[,1:4]

# plot each time-series by year
quebec %>%
  ggplot(aes(year, census, color = colony.id)) +
  geom_point()+
  geom_line(lwd=2)

quebec_blah <- quebec %>%
  filter(colony.id == "Anticosti Island" | colony.id == "Bonaventure Island" | colony.id == "Bird Rocks")
quebec_blah <- data.frame(colony.id = quebec_blah$colony.id,
                          size = quebec_blah$census,
                          year = quebec_blah$year)

# combine datasets
total_gannets <- data.frame(colony.id = atlantic$Colony.Id, 
                            size = atlantic$Colony.size, 
                            year = atlantic$year) %>%
  filter(colony.id != "Machias Seal Island, NB" & colony.id != "Whitehorse Island, NB" & colony.id != "Green Island 1, NS") %>%
  add_row(quebec_blah) %>%
  complete(colony.id, year = 1934:2020, # add missing vals
           fill = list(size = NA))

total_gannets$colony.id[total_gannets$colony.id == "Baccalieu Island, NF"] <- "Baccalieu Island"
total_gannets$colony.id[total_gannets$colony.id == "Cape St. Marys, NF"] <- "Cape St. Marys"
total_gannets$colony.id[total_gannets$colony.id == "Funk Island, NF"] <- "Funk Island"

# graph by colony
total_gannets %>%
  filter(year >= 1972) %>%
  ggplot(aes(year, pop_size))+
  geom_point(size = 2)+
  geom_line(size = 1.5)

#### FIGURE 5 GRAPH
na.omit(total_gannets) %>%
  filter(year >= 1968) %>%
  ggplot(aes(year, size))+
  geom_point(size=2, aes(color = colony.id))+
  geom_line(lwd=1.5, aes(color = colony.id))+
  facet_wrap(vars(colony.id), scales = "free_y")+
  #geom_smooth(method="lm")+
  scale_color_manual(values = c("#ff595e","#ff7d00","#ffca3a","#8ac926","#1982c4","#6a4c93"))+
  labs(x="Year",
       y = "Population Size")+
  theme_classic()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"),
    legend.position = "none",
    strip.text.x = element_text(size = 20),
    plot.margin = margin(10,10,10,10))
  
# creating new figure 4 graph (combined with fig. 3 for comparison)
# creating new data.frame
total_gannets_alt <- data.frame(
  year = c(1972, 1984, 1994, 1999, 2004, 2009, 2010, 2011, 2012, 2013, 2018),
  pop = c(32732, 40102, 56993, 77735, 101846, 116825, 100992, 104788, 104032, 108322, 113524)
)

all_years <- data.frame(year = 1972:2018)

result <- all_years %>% 
  left_join(total_gannets_alt, by = "year")

b <- na.omit(result) %>%
  ggplot(aes(year, pop)) +
  geom_smooth(method = "lm", color = "black")+
  geom_point(size = 3, color = "#3A4F41")+
  geom_line(lwd = 1.5, color = "#3A4F41")+
  labs(x = "Year",
       y = "Population Size")+
  scale_x_continuous(limits = c(1970,2020))+
  scale_y_continuous(labels = scales::comma)+
  theme_classic()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 18))
b

ggsave("figure5_alt.png", width = 15, height = 10, dpi = 600)

##################################################################################################################################

#### Calculate Rate of Increase

# fit linear model to data

# Anticosti Island
anticosti <- total_gannets %>%
  filter(colony.id == "Anticosti Island") %>%
  filter(between(year, 1972, 2018))
ggplot(anticosti, aes(year, size))+
  geom_point()+
  geom_line()+
  geom_smooth(method = "lm")
fitlm <- lm(anticosti$size ~ anticosti$year)
predlm <- predict.lm(fitlm, newdata = data.frame(anticosti$year),
                     type = "response")
anti.mean <- mean((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 0.1% change, insignificant
hist((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)])
anti.med <- median((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 0.19

# Baccalieu Island
baccalieu <- total_gannets %>%
  filter(colony.id == "Baccalieu Island") %>%
  filter(between(year, 1972, 2018))
ggplot(baccalieu, aes(year, size))+
  geom_point()+
  geom_line()+
  geom_smooth(method = "lm")
fitlm <- lm(baccalieu$size ~ baccalieu$year)
predlm <- predict.lm(fitlm, newdata = data.frame(baccalieu$year),
                     type = "response")
bacc.mean <- mean((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 7.97% change, significant
hist((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)])
bacc.med <- median((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 4.0

# Bird Rocks
birdrocks <- total_gannets %>%
  filter(colony.id == "Bird Rocks") %>%
  filter(between(year, 1972, 2018))
ggplot(birdrocks, aes(year, size))+
  geom_point()+
  geom_line()+
  geom_smooth(method = "lm")
fitlm <- lm(birdrocks$size ~ birdrocks$year)
predlm <- predict.lm(fitlm, newdata = data.frame(birdrocks$year),
                     type = "response")
brocks.mean <- mean((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 6.64% change, significant
hist((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)])
brocks.med <- median((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 3.70

# Bonaventure Island
bonaventure <- total_gannets %>%
  filter(colony.id == "Bonaventure Island") %>%
  filter(between(year, 1972, 2018))
ggplot(bonaventure, aes(year, size))+
  geom_point()+
  geom_line()+
  geom_smooth(method = "lm")
fitlm <- lm(bonaventure$size ~ bonaventure$year)
predlm <- predict.lm(fitlm, newdata = data.frame(bonaventure$year),
                     type = "response")
bon.mean <- mean((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 2.95% change, significant
hist((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)])
bon.med <- median((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 2.50

# Cape St. Marys
stmarys <- total_gannets %>%
  filter(colony.id == "Cape St. Marys") %>%
  filter(between(year, 1972, 2018))
ggplot(stmarys, aes(year, size))+
  geom_point()+
  geom_line()+
  geom_smooth(method = "lm")
fitlm <- lm(stmarys$size ~ stmarys$year)
predlm <- predict.lm(fitlm, newdata = data.frame(stmarys$year),
                     type = "response")
csm.mean <- mean((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 3.35% change, significant
hist((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)])
csm.med <- median((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 2.72

# Funk Island
funk <- total_gannets %>%
  filter(colony.id == "Funk Island") %>%
  filter(between(year, 1972, 2018))
ggplot(funk, aes(year, size))+
  geom_point()+
  geom_line()+
  geom_smooth(method = "lm")
fitlm <- lm(funk$size ~ funk$year)
predlm <- predict.lm(fitlm, newdata = data.frame(funk$year),
                     type = "response")
funk.mean <- mean((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 2.70% change, significant
hist((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)])
funk.med <- median((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]) # 2.34

meds <- c(funk.med, csm.med, bon.med, brocks.med, bacc.med, anti.med)
means <- c(funk.mean, csm.mean, bon.mean, brocks.mean, bacc.mean, anti.mean)
mean(means)

# total average, across all islands?
# mean(c(2.70,3.35,2.95,6.64,7.97, 0.01)) # 3.93%
weight <- c(10964,14598,52417,26234,3488,96)
sum(2.70*.101,3.35*.135,2.95*.486,6.64*.243,7.97*.032,0.01*0.0008)

#95%CIs with means
mean(means)
round(quantile(means, probs = c(0.025, 0.975)), 4)

weighted.mean(means, weight)
Hmisc::wtd.quantile(means, weight, probs = c(0.025, 0.975))

#95%CIs with medians
round(median(meds), 4)
round(quantile(meds, probs = c(0.025, 0.975)), 4)

# weighted
round(matrixStats::weightedMedian(meds,weight), 4)
round(Hmisc::wtd.quantile(meds, weight, probs = c(0.025, 0.975)), 4)

