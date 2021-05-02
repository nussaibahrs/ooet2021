# This script calculates the preferences of each marine plaknton group_by.

# Load libraries and functions ----------------------------------------------------------
library(tidyverse)
library(divDyn)
library(foreach)
library(doParallel)

source(file.path("scripts", "utils","functions.R"))
source(file.path("scripts","utils", "plot_parameters.R"))

# Load data ---------------------------------------------------------------
load(file.path("data", "wrangled_data.RData")) # sp at genus level

dat[dat$fossil.group == "FALSE",]$fossil.group <- "F"
dat <- dat %>% filter(FAD != LAD)

dat<- dat %>% mutate(per = case_when(mid > 50 ~ "G",
                                     mid < 50 & mid > 32 ~ "T1",
                                     mid < 32 & mid > 14 ~"T2",
                                     mid < 14 ~ "I")) 


dat <- dat %>% bind_rows(dat %>% mutate(fossil.group = "all"))

dat <- dat %>% bind_rows(dat %>% mutate(per = "all"))

dat$per <- factor(dat$per, levels=c("all","G", "T1", "T2", "I"))

# Bootstrap odds ratio -------------------------------------------

it=1000
chunk <- 20
n <- it/chunk


geog.ev <- levels(dat$per)
prefs <- list()

#progress bar
pb <- txtProgressBar(min = 0, max = chunk, style = 3)


for(f in unique(dat$fossil.group)){
  df <- dat %>%  filter(fossil.group==f)
  
  for(i in 1:chunk){
    # update progress bar
    setTxtProgressBar(pb, i)
    
    #create replicates of subsampled data based on number of iterations
    sub.df <- replicate(n, bootstr(df, "sqs"), simplify = FALSE)
    
    finalCoefs <- lapply(sub.df, function (x){
      temp <- df[x,]
      
      temp_ori <- temp[,c("ori", "env1", "sp", "bin")]
      temp_ori <- unique(temp_ori)
      
      temp_ext <- temp[,c("ext", "env1", "sp", "bin")]
      temp_ext <- unique(temp_ext)
      
      bins <- sort(unique(c(temp_ori$bin, temp_ext$bin)))
      bins <- bins[bins > 80]
      
      prefs <- list()
      
      for (t in bins){
        tryCatch({
          mod1 <- glm(ori~env1, data=temp_ori[temp_ori$bin==t,], family = binomial())
          mod2 <- glm(ext~env1, data=temp_ext[temp_ext$bin==t,], family = binomial())
          
          
          prefs[[t]] <- c(orig=coef(mod1)[2], ext=coef(mod2)[2], bin=t)
        }, error=function(e)prefs[[t]] <- c(orig=NA, ext=NA, bin=t))
      }
      
      return(do.call(rbind,prefs))
    })
    
    prefs[[paste(f,i)]] <- setNames(do.call(rbind.data.frame, finalCoefs), 
                                    c("orig1t", "ext1t", "bin")) %>%  
      mutate(group=f)
    
  }
}

prefs2 <- do.call(rbind.data.frame, prefs)
prefs2[,c("orig.env1t", "ext.env1t")]  <- apply(prefs2[,c("orig1t", "ext1t")] ,2,as.numeric)

save(prefs2, file=file.path("output", "prefs_sp_stage.RData"))


# Plot --------------------------------------------------------------------
load(file.path("output", "prefs_sp_stage.RData"))

p1 <- prefs2 %>% 
  group_by(group, bin) %>% 
  do(boot.summary2(.$orig.env1t)) %>% 
  left_join(int.tbl, by=c("bin"="stg")) %>% 
  ggplot(aes(x=mid, y=mean, ymin=lwr, ymax=upr, col=group, fill=group)) +
  geom_hline(yintercept = 0, linetype="dashed", col="darkgrey")+
  geom_ribbon(alpha=0.3, col=NA)+
  geom_line()+
  geom_point() +
  scale_color_manual(values=clr)+
  scale_fill_manual(values=clr)+
  scale_y_continuous(trans="reverse")+
  scale_x_continuous(trans="reverse", limits=c(66,0))+
  labs(x="Age(Ma", y="Origination Preferences")+
  facet_wrap(~group, labeller = labeller(group=group_names), nrow=2,
             scales="free_y") +
  theme_hce +
  theme(legend.position = "none")

y1 <- 16
y2 <- 20
p1 <- p1 + annotate("rect", xmin=syscol$bottom, xmax=syscol$top, ymin=y1, ymax=y2, col="white", fill="grey95") +
  annotate("text", x=syscol$mid[2], y=(y1+y2)/2, 
           label=syscol$system[2], hjust=0.5, fontface="bold", 
           col="grey20", size=3)+
  annotate("text", x=syscol$mid[3], y=(y1+y2)/2, 
           label=syscol$system[3], hjust=0.5, fontface="bold", 
           col="grey20", size=3)

p2 <- prefs2 %>% 
  group_by(group, bin) %>% 
  do(boot.summary2(.$ext.env1t)) %>% 
  left_join(int.tbl, by=c("bin"="stg")) %>% 
  ggplot(aes(x=mid, y=mean, ymin=lwr, ymax=upr, col=group, fill=group)) +
  geom_hline(yintercept = 0, linetype="dashed", col="darkgrey")+
  geom_ribbon(alpha=0.3, col=NA)+
  geom_line()+
  geom_point() +
  scale_color_manual(values=clr)+
  scale_fill_manual(values=clr)+
  scale_y_continuous(trans="reverse")+
  scale_x_continuous(trans="reverse", limits=c(66,0))+
  labs(x="Age(Ma", y="Extinction Preferences")+
  facet_wrap(~group, labeller = labeller(group=group_names), nrow=2,
             scales="free_y") +
  theme_hce +
  theme(legend.position = "none")

y1 <- 16
y2 <- 20

p2 <-p2 + annotate("rect", xmin=syscol$bottom, xmax=syscol$top, ymin=y1, ymax=y2, col="white", fill="grey95") +
  annotate("text", x=syscol$mid[2], y=(y1+y2)/2, 
           label=syscol$system[2], hjust=0.5, fontface="bold", 
           col="grey20", size=3)+
  annotate("text", x=syscol$mid[3], y=(y1+y2)/2, 
           label=syscol$system[3], hjust=0.5, fontface="bold", 
           col="grey20", size=3)

svg(file.path("figs", "Supplement", "Fig_S_preferences_stage.svg"), 
    w=10, h=10)
p1+p2 + 
  plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")") +
  plot_layout(nrow=2)
dev.off()

