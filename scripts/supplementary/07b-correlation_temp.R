# Computes the correlation between diversification and dispersal, wit
# temperature.

# load libraries ---------------------------------------------------------
library(divDyn)
library(magrittr)
library(reshape2)
library(forecast)
library(nlme)


# Load data ---------------------------------------------------------------

load(file.path("data", "wrangled_data.RData")) # sp at species level
load(file.path("output", "turnover_rates.RData"))
load(file.path("output", "mig_stage.RData"))

# Correlation with temperature --------------------------------------------

source(file.path("code", "supplementary", "07a-temp_conversion_paleo.R"))
paleotemp <- df
colnames(paleotemp) <- c("age", "temperature")

#bin according to stage
paleotemp$stg <- NA

for (i in 1:nrow(stages)){
  temp <- stages[i,]
  
  paleotemp$stg[which(paleotemp$age <= temp$bottom & paleotemp$age > temp$top)] <- temp$stg
}

paleotemp <- tapply(paleotemp$temperature, paleotemp$stg, mean, na.rm=TRUE)
paleotemp <- data.frame(paleotemp) %>% 
  mutate(stg=as.numeric(rownames(.))) %>% 
  left_join(int.tbl %>% dplyr::select(stg,mid), by="stg")


acf(paleotemp$paleotemp)

orig.rate <- orig.rate %>% 
  filter(mid <67 & method =="sqs") %>% 
  select(-method) %>% 
  reshape2::melt(id.vars=c("mid", "env", "group")) %>% 
  group_by(mid, env, group) %>% 
  summarise(mean=mean(value, na.rm=TRUE))

ext.rate <- ext.rate %>% 
  filter(mid <67 & method =="sqs") %>% 
  select(-method) %>% 
  reshape2::melt(id.vars=c("mid", "env", "group")) %>% 
  group_by(mid, env, group) %>% 
  summarise(mean=mean(value, na.rm=TRUE))

mig.tot <- mig.tot %>% left_join(stages[,c("stg", "mid")], by=c("bin"="stg")) %>% 
  filter(mid <67 & method =="sqs") %>% 
  select(-method, -env,-bin) %>% 
  reshape2::melt(id.vars=c("mid", "group"), variable.name="env") %>% 
  group_by(mid, env, group) %>% 
  summarise(mean=mean(value, na.rm=TRUE))

# * all -------------------------------------------------------------------
rate <- rbind(cbind(orig.rate , rate="origination"),
              cbind(ext.rate,  rate="extinction"),
              cbind(mig.tot, rate="migration"))

corr <- list()

for(g in names(group_names)){
  for(e in c("t", "nt")){
    for(t in c("origination", "extinction", "migration")){
      temp <- rate %>% filter(group==g & env==e & rate==t) %>% 
        left_join(paleotemp, by="mid") %>% 
        filter(!(is.na(mean) | is.na(paleotemp))) 
      
      #linear regression with raw data
      mod <- lm(mean~paleotemp, data=temp)
      
      #find best ARIMA model from residuals
      AR <- forecast::auto.arima(mod$residuals)
      
      #extract parameters
      Phi <- AR$model$phi
      Phi<- ifelse(length(Phi) > 0, Phi, 0)  
      p <- length(Phi)
      
      Theta <- AR$model$theta
      Theta<- ifelse(length(Theta) > 0, Phi, 0)  
      q<- length(Theta)
      
      #fit model considering the autocorrelation struction from ARIMA
      mod.cor <- nlme::gls(formula(mod), data=temp, correlation = nlme::corARMA(c(Phi, Theta), form = as.formula(~1), 
                                                                                p=p, q=q),
                           method="ML")
      
      
      #save data
      corr[[paste(g,e, t)]]  <- c(group=g, env=e, rate=t, Phi = coef(mod.cor$modelStruct$corStruct,unconstrained=FALSE)[1],
                                  Intercept=summary(mod.cor)$tTable[1,"Value"],
                                  summary(mod.cor)$tTable[2,c("Value", "p-value")])
    }
  }
}

corr <- corr %>%  bind_rows()
corr$Intercept <- as.numeric(corr$Intercept)
corr$Value <- as.numeric(corr$Value)
corr$`p-value`<- as.numeric(corr$`p-value`)
corr$Phi.Phi1 <- as.numeric(corr$Phi.Phi1)
alpha <- 0.05

corr$significant <- ifelse(corr$`p-value` < alpha, 1,0)

corr <- corr %>%
  mutate(group = recode(group, !!!as.list(group_names)),
         env = recode(env, !!!as.list(env_names)),
         rate = factor(rate, levels=c("origination", "extinction", "migration"))) %>% 
  arrange(rate) %>% 
  set_names(c("Fossil Group", "Climate Zone", "Rate", "Phi","Intercept", "Coefficient", "p-value","is.significant")) %T>% 
  write.csv(file.path("output", "temperature_correlation.csv"), row.names = FALSE)
