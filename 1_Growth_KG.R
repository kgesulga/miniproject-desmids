rm(list=ls())


#######

# packages

library(nlsMicrobio)
library(data.table)
library(plyr)
library(tools)
library(zoo)
library(lattice)
library(ggplot2)
library(MuMIn)
library(progress)
#install.packages("devtools")
library(devtools)
library(cowplot)
# install.packages("ghit")
# library(ghit)
# install_github('padpadpadpad/TeamPhytoplankton')
# install_github('padpadpadpad/nlsLoop')
library(TeamPhytoplankton) #ask Elisa for source code if it won't load
library(nlsLoop)#ask Elisa for source code if it won't load

# dir
setwd("C:/Users/User/Desktop/Desmids/4_analyses")

# read in count data 
tt<-read.csv("data/dat_growthrate_edit2_mumax.csv") #%>%  # KG, sep=";") # LL on any machine where english isn't operating language, specify separator!
  # KG - filter days with massive outliers 
  # filter(!Day %in% c(6,15,27))
# KG - check if filtered days correctly
unique(tt$Day)

#examine the file - structure (str) and header (head)
str(tt)
head(tt)
tt$Day<-as.numeric(tt$Day)
tt$BioRep<-as.factor(tt$BioRep)
tt$Nutrient<-as.factor(tt$Nutrient)
tt$community<-as.factor(tt$community)
tt$Temp<-as.factor(tt$Temp)
#tt$Day<-tt$Day-2
tt$PFU<-as.numeric(tt$PFU)
#tt<-subset(tt,Day<11&Day>0)



tt <- within(tt, id.01 <- as.character(factor(community):factor(Nutrient):factor(Temp):factor(BioRep)))# create the identifier 
#we need to make a column with the counts as log values
tt$logcounts<-log10(tt$PFU)

#we also get rid of all NAs
tt<-tt[complete.cases(tt), ]

#check if it did the one job it had
head(tt) # YES 

tts<- tt[,c('Day','logcounts','id.01')]# we only keep the column for day, log counts, and the identifier
colnames(tts) <-c('t','LOG10N','id.01') # rename columns so they fit the model parameters as R likes them
#colnames(tts) <-c('t', 'gm', 'id.01')
id.01 <- unique(tts$id.01) # find all region/ident combinations that are unique (so we pool by bioreplicate)

head(tts)


#now make a first quick plot, inspect the data..

ggplot(tt, aes(x = Day, y = logcounts, colour=Temp, shape=Nutrient,alpha=BioRep)) +
  geom_point()+
  geom_smooth(se=F, size=0.5)+
  #scale_colour_manual("Temperature",values=c("#1B9E77","#666666", "#7570B3", "#E7298A" ,"#66A61E"))+
  theme_classic(base_size = 14, base_family = "Helvetica")+
  facet_wrap(Nutrient~community)
    # i don't care that alpha for factors is not advised. Shuddup R. 


#this is the function for the most simple logistic model
growth <- function(r, K, t, T0.biom){
  gm <- K /(1 + ((K - T0.biom)/min(T0.biom))*exp(-r*t))
  return(gm)
}


#get time at which there is µmax - note that you may need to add lag phase to that later on. 

atmu<-function (r, K, T0.biom = 2) 
{
    return(t<-((T0.biom/(1-T0.biom/K)))/r)
}


##### run a model with propeR lag ####

res_logis_lag <- nlsLoop(tts, #this will throw hissy fits until it finds a good fit! Error messages are normal as not all of the curves it tries out are going to fit
                     model = gompertzm,
                     tries = 100, #run this 1000 times 
                     id_col = 'id.01', #id can be collapsed at will. right now fits to each biorep, isolate, region, and assayT. The more you collapse the id (you do that above), the fastr your model will run
                     param_bds = c(2,10,0.5,35,0.01,3,0.5,8), #parameter space -this defines values that initial density, K, µmax and lag are allowed to take. see parameter space below for details
                     na.action = na.omit,
                     lower = c(LOG10N0=2,LOG10Nmax=0.5,mumax=0.01,lag=0.5), #lower parameter space for starting density LOG10NO (i.e. starting cells as log10), carrying capacity here called LOG10Nmax,mumax for µmax, and lag for the lag phase
                     upper = c(LOG10N0=10,LOG10Nmax=35,mumax=3,lag=8)) #upper parameter space for starting density LOG10NO (i.e. starting cells as log10), carrying capacity here called LOG10Nmax,mumax for µmax, and lag for the lag phase

#res_logis_lag <- nlsLoop(tts, #this will throw hissy fits until it finds a good fit! Error messages are normal as not all of the curves it tries out are going to fit
                         #model = gm ~ growth(K, r, t, T0.biom),
                         #tries = 50, #run this 1000 times 
                         #id_col = 'id.01', #id can be collapsed at will. right now fits to each biorep, isolate, region, and assayT. The more you collapse the id (you do that above), the fastr your model will run
                         #param_bds = c(1,20,0.1,2.1,0.5,8), #parameter space -this defines values that initial density, K, µmax and lag are allowed to take. see parameter space below for details
                         #na.action = na.omit,
                         #lower = c(K=1,r=0.1,T0.biom=0.001), #lower parameter space for starting density LOG10NO (i.e. starting cells as log10), carrying capacity here called LOG10Nmax,mumax for µmax, and lag for the lag phase
                         #upper = c(K=20,r=2.1,T0.biom=12)) #upper parameter space for starting density LOG10NO (i.e. starting cells as log10), carrying capacity here called LOG10Nmax,mumax for µmax, and lag for the lag phase


res_logis_lag$params  #these are the parameters 
# KG check these to explore if any sample has reached an upper or lower boundary
# KG if so, adjust parameters

growth_params_lag<-res_logis_lag $params 
growth_params_lag$atmum<-atmu(r=growth_params_lag$mumax, K=growth_params_lag$LOG10Nmax,T0.biom=3)+growth_params_lag$lag
head(res_logis_lag$predictions) #these are the predictions that we will use to fit the curve 
growth_preds_lag<-res_logis_lag$predictions 

#now we want to get our temp and region and isolate back into dataframe 
growth_params_lag$community<- as.vector(t(sapply(as.character(growth_params_lag$id.01), function(y) strsplit(y,split=":")[[1]][1])))
growth_params_lag$Nutrient<- as.vector(t(sapply(as.character(growth_params_lag$id.01), function(y) strsplit(y,split=":")[[1]][2])))
growth_params_lag$Temp<- as.vector(t(sapply(as.character(growth_params_lag$id.01), function(y) strsplit(y,split=":")[[1]][3])))
growth_params_lag$BioRep<- as.vector(t(sapply(as.character(growth_params_lag$id.01), function(y) strsplit(y,split=":")[[1]][4])))

#check that this worked
head(growth_params_lag) #YAY

#now we also get the details back for the predictions 
growth_preds_lag$community<- as.vector(t(sapply(as.character(growth_preds_lag$id.01), function(y) strsplit(y,split=":")[[1]][1])))
growth_preds_lag$Nutrient<- as.vector(t(sapply(as.character(growth_preds_lag$id.01), function(y) strsplit(y,split=":")[[1]][2])))
growth_preds_lag$Temp<- as.vector(t(sapply(as.character(growth_preds_lag$id.01), function(y) strsplit(y,split=":")[[1]][3])))
growth_preds_lag$BioRep<- as.vector(t(sapply(as.character(growth_preds_lag$id.01), function(y) strsplit(y,split=":")[[1]][4])))


#check that this worked
head(growth_preds_lag)  # YAY!! 
colnames(growth_preds_lag)[3] <- 'pred_rate'


#######plot with fitted curves with proper lag#####
qplot(as.numeric(Day), logcounts, data=tt, colour=Temp, alpha=BioRep)+
  #scale_color_manual("Assay Temperature",values=c("#1B9E77","#666666", "#7570B3", "#E7298A" ,"#66A61E"))+
  theme_classic(base_size = 14, base_family = "Helvetica")+
  facet_wrap(Nutrient~community)+
  geom_line(data = growth_preds_lag, aes(x = t, y = pred_rate, alpha=as.factor(BioRep),colour = as.factor(Temp)))+
  theme(legend.position='top')  # YES MUCH BETTER
#scale_x_discrete(name="Day of experiment", breaks=seq(0, 10,2)) # you will need to change 20 for the number of experimental days 


#write your new data as a csv#
write.csv(growth_params_lag,"Growth_community_KG_edit2_mumax.csv")

####Plot Growth rates
head(growth_params_lag)

p3 <- 
  ggplot(growth_params_lag, aes(x = Temp ,y = mumax, color=Temp, shape=Nutrient)) +
  geom_point(size=2)+
  stat_summary(fun=mean,
               geom="point",
               #position=dodge.posn,
               show.legend = T, 
               cex=2, color="black")+
  stat_summary(fun.data=mean_se,
               geom="errorbar",
               #position=dodge.posn,
               show.legend = F,
               width=0.2,
               cex=0.5, color="black")+
  #scale_color_manual(values=c("#1B9E77","#666666", "#7570B3", "#E7298A" ,"#66A61E"))+
  #scale_linetype_manual(values=c(1,2,3))+
  #ylim(c(0,4))+
  labs(title="maximum growth rate", cex=1)+
  facet_wrap(~community*Nutrient)


png("plots/0_mumax/0_final_edit2_mumax_i3.png", width = 10, height = 6, units = "in", res = 300)
# Print the list of plots
print(p3)
# Close the PDF device
dev.off()


## KG: This doesn't work -- haven't had the time to fuss with it
ggplot(growth_params_lag, aes(x = Temp ,y = LOG10Nmax, color=Temp)) +
  geom_point(size=2)+
  stat_summary(fun=mean,
               geom="point",
               #position=dodge.posn,
               show.legend = T,
               cex=2, color="black")+
    stat_summary(fun.data=mean_se,
               geom="errorbar",
               #position=dodge.posn,
               show.legend = F,
               width=0.2,
               cex=0.5, color="black")+
    #scale_color_manual(values=c("#1B9E77","#666666", "#7570B3", "#E7298A" ,"#66A61E"))+
  #scale_linetype_manual(values=c(1,2,3))+
  ylim(c(6.3,7))+
  labs(title="Carrying capacity", cex=1)+
  facet_wrap(~community*Nutrient)
