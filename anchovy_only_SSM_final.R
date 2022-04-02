# CALCULATE MISSING ANCHOVY DATA FOR ALL STATIONS FOR 2018-2019 ----

load("Mars_Anchovy/CentralBay_MARSS_2022.RData")

#### Set Working Directory ----
#setwd("")

#### Load Libraries ----
library(tidyverse)
library(MARSS)
library(broom)
library(forecast)
library(lubridate)
library(ggplot2)

#### Load and Modify Data ----

### Load Data
CentralBayFish <- read.csv("Mars_Anchovy/input/CentralBayFish.csv") #Have one date and ordinal date for each survey

### Vector modification
CentralBayFish$Year <- as.integer(CentralBayFish$Year)
CentralBayFish$fYear <- as.factor(CentralBayFish$Year)

CentralBayFish$Date <- mdy(CentralBayFish$Date) #Convert Date column to Date vector

CentralBayFish$Survey <- as.factor(CentralBayFish$Survey)
CentralBayFish$fStation <- as.factor(CentralBayFish$Station)

CentralBayFish$NORANC <- as.numeric(CentralBayFish$NORANC)


str(CentralBayFish)

### Subset data ----
subset <- CentralBayFish %>%
  filter(Year>1987) %>%
  select(Year, fYear, Survey, Time, RepDate, Date, fStation, NORANC)

# Add 1 to fish count so can take log (can't take log of 0)
subset$NORANC <- subset$NORANC+1

#### Reshape into Matrix ----
anc <- reshape2::acast(subset, fStation ~ Year + Survey, value.var = "NORANC")

# Take log of matrix
anc_log <- log(anc)


#### MARSS Modeling: NORANC ----

# State process:
# x[t] = x[t-1] + u + w[t] where w[t]~ N(0,q)
# x[0]=intercept

# Observation process:
#y[i,t] = x[t] + a[i] + v[i] where v[i]~N(0,r)

### Set up seasonal effects

# Months x Years
TT <- dim(anc_log)[2]

# Number of months per year
period <- 12

# Seasonal effect as Fourier series
# Define the 2 x T seasonal covariate matrix as a combination of 1 cosine and 1 sine wave
cos.t <- cos(2 * pi * seq(TT)/period)
sin.t <- sin(2 * pi * seq(TT)/period)
c.Four <- rbind(cos.t, sin.t)

# Z matrix: Each station as one "population"
Z <- diag(11)

#Z <- as.matrix(rep(1, 11))

# Set up MARSS model parameters
mod.list = list(
  B = "identity",
  U = "unequal",  # Allow each station to have own trend
  C = "unconstrained", # Select unconstrained when inputing a custom covariate matrix
  c = c.Four,
  Q = "equalvarcov",
  Z = Z,
  A = "scaling",
  R = "unconstrained",
  x0 = "unequal",
  D = "zero",
  d = "zero",
  tinitx =0)

m5 <- MARSS(anc_log, model=mod.list, control=list(maxit=10000))

#save.image("D:/Telework/R_Directory/Alcatraz/MARSS/CentralBay_MARSS_2022.RData")

summary(m1$model) # information on the structure of the MARSS model that was fit
library(forecast)

x_num <- as.numeric(unlist(m5$fitted.ytT[2,]))

fr <- forecast(x_num)
ggplot2::autoplot(fr)

AIC(m1,m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m14, m16)

m5_CIs <- MARSSparamCIs(m5)

plot(m5, plot.type="model.residuals")

ggplot2::autoplot(m5)

# Tidy data
anc_est <- broom::tidy(m5, type="fitted.ytT", conf.int=TRUE, conf.level = 0.85) # Estimated population size
#Available types:
# "parameters"
# "xtT" Expected value of X (states) conditioned on all data
# ytT expected values of the y (left side of the y equation with the error terms)
# fitted.ytT expected values of the y (left side of the y equation withOUT the error terms)



# Plot mean LOGGED estimates
ggplot(data = anc_est) + 
  geom_line(aes(t, estimate)) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.2, fill="blue")+
  labs(x ="Time", y ="Estimate (logged)")
#+facet_wrap(~.rownames)

# Exponentiate to get "real" values
anc_exp <- anc_est
anc_exp$y <- exp(anc_exp$y)
anc_exp$estimate <- exp(anc_exp$estimate)
anc_exp$conf.low <- exp(anc_exp$conf.low)
anc_exp$conf.high <- exp(anc_exp$conf.high)
anc_exp$Year <- rep(1988:2019, each = 12, times = 1)
anc_exp$Month <- rep(1:12, 32)
anc_exp$Day <- 15

ggplot(data = anc_exp, aes(x = t, y = estimate)) + 
  geom_line(aes(t, estimate)) +
  #geom_point(data = fitted, aes(t, .fitted)) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.4, fill="red")+
  ylim(0, 6000)+
  facet_wrap(~.rownames)+
  theme_classic()

# Get summed total estimate for Central Bay
summed <- anc_exp %>% 
  group_by(t) %>% 
  mutate(mo_total = sum(estimate, na.rm = TRUE))%>% 
  mutate(mo_totalk = mo_total/1000) %>%
  ungroup(t) %>% 
  group_by(Year) %>% 
  mutate(yr_total = sum(estimate, na.rm = TRUE)) %>% 
  mutate(yr_totalk = yr_total/1000)

summed_orig <- CentralBayFish %>% 
  group_by(Year, Survey) %>% 
  mutate(mo_total = sum(NORANC, na.rm = TRUE))%>% 
  mutate(mo_totalk = mo_total/1000)%>%
  ungroup(Survey) %>% 
  group_by(Year) %>% 
  mutate(yr_total = sum(NORANC, na.rm = TRUE)) %>% 
  mutate(yr_totalk = yr_total/1000)


summed_orig <- summed_orig[,c(1:14,85,163:166)]
summed_orig$RepDate <- mdy(summed_orig$RepDate) #Convert Date column to Date vector

#write.csv(summed, "output/anc_predictions_Mar22.csv", row.names = FALSE)


#### Plot anchovy data ----
total_anc_m5 <- ggplot(summed, aes(x = Year, y = yr_totalk))+
  geom_line(col='blue', size=1)+
  geom_point(size = 1.5)+
  geom_line(summed_orig, mapping = aes(x = Year, y = yr_totalk), color = "red")+
  labs(x = "Year", y = "Annual Central Bay Anchovy Population (in thousands)") +
  theme_classic() +
  theme(axis.line = element_line(size = 1),   # These must be added after overall plot theme
        axis.title.x = element_text(size=12, face = "bold"), 
        axis.title.y = element_text(size=12, face = "bold"),
        axis.text.x = element_text(size=11, face="bold", color ="black"),
        axis.text.y = element_text(size=11, face="bold", color ="black"),
        axis.ticks = element_line(size=1))+
  scale_y_continuous(limits=c(0, 160), breaks = seq(0, 150, by =25), expand = c(0,0))+
  scale_x_continuous(limits=c(1988,2020.5), breaks =seq(1988,2020, by =2), expand = c(0,0))
total_anc_m5


library('Cairo')
png(filename="total_anchovy.png",
    type="cairo",units="px", width=2500, height=1500,  pointsize=12, res=300)
print(total_anc )
dev.off()


#### Interpolate Fish ----

### NORTHERN ANCHOVY

# Convert month and year to date
anc_int <- summed %>%
  mutate(date = make_date(Year,Month,Day))

spline_anc <- splinefun(x = anc_int$date, y = anc_int$mo_totalk)
Date <- seq.Date(ymd(min(anc_int$date)), ymd(max(anc_int$date)), by = 1)
SplineFit_anc <- spline_anc(Date)

plot(anc_int$date, anc_int$mo_totalk, ylim = c(0, 50),  xlab = "Date", ylab="Anchovy Population Estimate (in thousands)")
lines(Date, SplineFit_anc, col = "blue")

anc_splined <- as.data.frame(cbind(Date, SplineFit_anc))
anc_splined$Date <- as_date(anc_splined$Date)
colnames(anc_splined)[2] <-"IntTotalk" # Change column name to Interpolated Total
anc_splined$IntTotal <- anc_splined$IntTotalk*1000

# Replace values <0 with 1
anc_splined$IntTotal <- ifelse(anc_splined$IntTotal<0, 1, anc_splined$IntTotal)

# Take integral
anc_splined$Integral <- 1:length(anc_splined[,1])
for (i in 2:length(anc_splined[,1])) {
  # A3 - A2
  anc_splined$Integral[i] <- ((anc_splined[i, 1] - anc_splined[i - 1, 1]) *
  # B3 + B2
    (anc_splined[i, 3] + anc_splined[i - 1, 3])) / 2
}



# Original Data
spline_anc_orig <- splinefun(x = summed_orig$RepDate, y = summed_orig$mo_totalk)
Date <- seq.Date(ymd(min(summed_orig$RepDate)), ymd(max(summed_orig$RepDate)), by = 1)
SplineFit_anc_orig <- spline_anc_orig(Date)

plot(summed_orig$RepDate, summed_orig$mo_totalk, ylim = c(0, 350),  xlab = "Date", ylab="Anchovy Population Estimate (in thousands)")
lines(Date, SplineFit_anc_orig, col = "blue")

anc_splined_orig <- as.data.frame(cbind(Date, SplineFit_anc_orig))
anc_splined_orig$Date <- as_date(anc_splined_orig$Date)
colnames(anc_splined_orig)[2] <-"IntTotalk" # Change column name to Interpolated Total
anc_splined_orig$IntTotal <- anc_splined_orig$IntTotalk*1000

# Replace values <0 with 1
anc_splined_orig$IntTotal <- ifelse(anc_splined_orig$IntTotal<0, 1, anc_splined_orig$IntTotal)

# Take integral
anc_splined_orig$Integral <- 1:length(anc_splined_orig[,1])
for (i in 2:length(anc_splined_orig[,1])) {
  # A3 - A2
  anc_splined_orig$Integral[i] <- ((anc_splined_orig[i, 1] - anc_splined_orig[i - 1, 1]) *
                                # B3 + B2
                                (anc_splined_orig[i, 3] + anc_splined_orig[i - 1, 3])) / 2
}



# Export 
#write.csv(anc_splined,"D:/Telework/R_Directory/Alcatraz/FishInt/output/splined_anchovy_integral2.csv", row.names = FALSE)


#### Plot splined anchovy data ----
library(cowplot)# for draw_image
library(magick)
library(patchwork)

splined_anc_plot <- ggplot(anc_splined, aes(x = Date, y = IntTotalk))+
  geom_line(col='blue', size=1)+
  #geom_point(size = 1.5)+
  geom_line(anc_splined_orig, mapping = aes(x = Date, y = IntTotalk), color = "red")+
  labs(x = "Year", y = "Daily Central Bay Anchovy Population (in thousands)") +
  theme_classic() +
  theme(axis.line = element_line(size = 1),   # These must be added after overall plot theme
        axis.title.x = element_text(size=12, face = "bold"), 
        axis.title.y = element_text(size=12, face = "bold"),
        axis.text.x = element_text(size=11, face="bold", color ="black"),
        axis.text.y = element_text(size=11, face="bold", color ="black"),
        axis.ticks = element_line(size=1))+
  scale_y_continuous(limits=c(0, 40), breaks = seq(0, 40, by =5), expand = c(0,0))+
  scale_x_date(date_breaks = "1 year", date_labels= "%Y", date_minor_breaks = "1 year", limits = as.Date(c('1988-01-01','2020-06-15'), by = "2 years"), expand = c(0,0))
splined_anc_plot

splined_anc_plot2 <- ggdraw(splined_anc_plot) +
  draw_image("D:/Telework/R_Directory/Alcatraz/MARSS/NorAnchovy.png", scale = 0.3, halign = 0.95, valign = 0.9)
splined_anc_plot2

library('Cairo')
png(filename="Model5_ANC_ytT.png",
    type="cairo",units="px", width=3800, height=2000,  pointsize=12, res=300)
print(splined_anc_plot2)
dev.off()
