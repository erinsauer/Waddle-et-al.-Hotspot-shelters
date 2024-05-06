####################################################################
######## This code was developed by Erin L. Sauer      #############
######## for: Waddle et al. "Hotspot shelters          #############
######## stimulate frog resistance to chytridiomycosis #############
####################################################################

library(tidyverse)
library(mgcv)
library(itsadug)

#load all the tpref data
TT1<- read.csv("After infection tpref_T1.csv")
TT2<- read.csv("After infection tpref_T2.csv")
CT1<- read.csv("Repeated not infected tpref_T1.csv")
CT2<- read.csv("Repeated not infected tpref_T2.csv")
BL<-read.csv("basline tpref.csv")
Meso <- read.delim(pipe("pbpaste"))

#convert each file from wide to long and merge them
colnames(TT1)
tsubs <- grep("tsub", names(TT1)) #remove sub temps for now
TT1<- TT1 %>% 
  select(-tsubs) %>% #remove sub temps for now
  add_column(DF = "TT1") %>% #add column with the dataframe name
  pivot_longer(cols = starts_with("tbody"), #convert wide to long
                     names_to="timepoint",
                     names_prefix="tbody",
                     values_to=c("tbody"),
                     values_drop_na =TRUE)
colnames(CT1)
tsubs <- grep("tsub", names(CT1))
CT1 <- CT1 %>% 
  select(-tsubs) %>%
  add_column(DF = "CT1") %>%
  pivot_longer(cols = starts_with("tbody"),
               names_to="timepoint",
               names_prefix="tbody",
               values_to=c("tbody"),
               values_drop_na =TRUE)

T1 <- rbind(TT1, CT1) #merge the rows
#make a unique value for each individual
T1$ID <- paste0(T1$Gradient..,T1$DF)
#modeling 
str(T1)
T1$timepoint <- as.numeric(T1$timepoint) 
T1$EXP.group <- as.factor(T1$EXP.group)
m1 <- gamm(tbody ~ EXP.group + 
             s(timepoint , k=3), # smoothing function
      family=gaussian, 
      correlation=corCAR1(form=~timepoint|ID), data = T1)

#results
plot(m1$gam, pages=1) #this plot lets us visually assess the fit
summary(m1$gam)
anova(m1$gam)

#plot showing the two treatments
plot_smooth(m1$gam, view="timepoint", plot_all="EXP.group", rug=F,  
            main="", cex.lab=2, cex.axis=1.5,
            ylab="Tpref", xlab="Dayes since exposure")
########
#Meso data
str(Meso)
library(tidyverse)
Meso$date <- dmy(Meso$date)
Meso$day <- yday(Meso$date)

Meso$Vaccinated.Naïve.2 <- as.factor(Meso$Vaccinated.Naïve.2)
Meso$control.or.chytrid. <- as.factor(Meso$control.or.chytrid.)
Meso$shaded.or.not. <- as.factor(Meso$shaded.or.not.)
par(mar = c(5, 5, 2, 5))
###body temp
Meso <- subset(Meso, Vaccinated.Naïve.2 == "naïve" | Vaccinated.Naïve.2 == "vaccinated")

m1 <- gamm(T.body ~ Vaccinated.Naïve.2 + control.or.chytrid. + shaded.or.not. +
             time.started +
             s(day , k=5), # smoothing function
           family=gaussian, 
           correlation=corCAR1(form=~day|PIT..), data = Meso)
plot(m1$gam, pages=1) #this plot lets us visually assess the fit
summary(m1$gam)
anova(m1$gam)

#plot showing the two treatments
plot_smooth(m1$gam, view="day", plot_all="shaded.or.not.", rug=F,  
            main="", cex.lab=2, cex.axis=1.5,
            ylab="Tpref", xlab="Day")
plot_smooth(m1$gam, view="day", plot_all="Vaccinated.Naïve.2", rug=F,  
            main="", cex.lab=2, cex.axis=1.5,
            ylab="Tpref", xlab="Day")
plot_smooth(m1$gam, view="day", plot_all="control.or.chytrid.", rug=F,  
            main="", cex.lab=2, cex.axis=1.5,
            ylab="Tpref", xlab="Day")
###load
Meso2 <- subset(Meso, control.or.chytrid. == "chytrid")
Meso2 <- subset(Meso2, Vaccinated.Naïve.2 == "naïve" | Vaccinated.Naïve.2 == "vaccinated")
#
m2 <- gamm(log(ITS.copy+1) ~ Vaccinated.Naïve.2 * shaded.or.not. + 
             s(day , k=5), # smoothing function
           family=gaussian, 
           correlation=corCAR1(form=~day|PIT..), data = Meso2)
plot(m2$gam, pages=1) #this plot lets us visually assess the fit
summary(m2$gam)
anova(m2$gam)

#plot showing the two treatments
plot_smooth(m2$gam, view="day", plot_all="shaded.or.not.", rug=F,  
            main="", cex.lab=2, cex.axis=1.5,
            ylab="load", xlab="Day")
#vax by shaded interaction plots
##greenhouse treatment only
plot_smooth(m2$gam, view="day", plot_all="Vaccinated.Naïve.2", 
            cond=list(shaded.or.not.="greenhouse"), rug=F,  
            main="", cex.lab=2, cex.axis=1.5,
            ylab="load", xlab="Day")
##shaded treatment only
plot_smooth(m2$gam, view="day", plot_all="Vaccinated.Naïve.2",
            cond=list(shaded.or.not.="shaded"), rug=F,  
            main="", cex.lab=2, cex.axis=1.5,
            ylab="load", xlab="Day")


