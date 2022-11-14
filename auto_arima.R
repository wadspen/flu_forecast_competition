library(lubridate)
library(forecast)
library(tidyr)
library(dplyr)
library(ggplot2)


raw_ILI <- read.csv('https://raw.githubusercontent.com/cdcepi/Flusight-forecast-data/master/data-truth/truth-Incident%20Hospitalizations.csv')

quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
levels <- 100*(quantiles[23:13]-quantiles[1:11])



# forecast_date <- ceiling_date(today(),unit='week',2)
forecast_date <- ymd('2022-11-14')
forecast_plots <- list()
all_forecasts <- data.frame()
m <- 1
locs <- unique(raw_ILI$location)
loc_names <- unique(raw_ILI$location_name)
for (i in 1:length(locs)) {
  mod <- auto.arima(ts(raw_ILI$value[raw_ILI$location==locs[i]]))
  forc <- forecast(mod,level=levels)
  # forecast_plots[[m]] <- autoplot(forecast(mod,4))
  ggp <- autoplot(forecast(mod,4)) + ggtitle(loc_names[i]) + 
          ylab('hopitalizations') +
          theme_bw()
  forecast_plots[[i]] <- ggp
  m <- m + 1
  low <- forc$lower[1:4,]
  low <- low[,order(colnames(low))]
  colnames(low) <- quantiles[11:1]
  upp <- forc$upper[1:4,]
  colnames(upp) <- quantiles[13:23]

  forecasts <- cbind(low,'0.5'=forc$mean[1:4],upp)
  forecasts <- data.frame(t(forecasts))
  forecasts$quantile <- rownames(forecasts)
  targets <- c('1 wk ahead inc flu hosp','2 wk ahead inc flu hosp',
    '3 wk ahead inc flu hosp','4 wk ahead inc flu hosp')
  
  colnames(forecasts)[1:4] <- targets
  l_forecasts <- forecasts %>% 
    pivot_longer(1:4,names_to='target') %>% 
    mutate(type='quantile') %>% 
    select(target,type,quantile,value)

  points <- l_forecasts %>% 
    group_by(target) %>% 
    summarise(
      value=median(value)
    ) %>% 
    mutate(type='point',quantile=NA) %>% 
    select(target,type,quantile,value)

  point_forecasts <- rbind(l_forecasts,points) %>% 
    arrange(target,type) %>% 
    mutate(location=i)
  
  all_forecasts <- rbind(all_forecasts,point_forecasts)
}

all_forecasts <- all_forecasts %>% 
  mutate(forecast_date,
         target_end_date=forecast_date + 7*as.numeric(substr(target,1,1))-2) %>% 
  select(forecast_date,target,target_end_date,location,type,quantile,value) 

all_forecasts <- all_forecasts %>% 
  mutate(value = ifelse(value >=0,value,0))

file_name <- paste(getwd(),'/ISU_NiemiLab/',
                   forecast_date,'-ISU_NiemiLab.csv',sep='')
# write.csv(all_forecasts,file_name,row.names = FALSE)

library(gridExtra)
pdf("plots.pdf", onefile = TRUE)
for (i in seq(length(forecast_plots))) {
  do.call("grid.arrange", args=list(top='Flu Forecast'),forecast_plots[[i]])  
}

do.call("grid.arrange", c(forecast_plots[1:12], nrow = 3)) 
do.call("grid.arrange", c(forecast_plots[13:24], nrow = 3))
do.call("grid.arrange", c(forecast_plots[25:36], nrow = 3))
do.call("grid.arrange", c(forecast_plots[37:48], nrow = 3))
do.call("grid.arrange", c(forecast_plots[49:length(forecast_plots)], nrow = 3))

