#This document represents all the code that was written for intermediate steps or methods that didn't end up making it 
#into the final analysis -- this has EVERYTHING in it that was used, but as such, it's not organized optimally

################################ 10/19/18 ##########################################

#Running some mixed models on the data according to directions from MK yesterday. 

#glmm models
#ambd models require random effect to be a factor, so transform collection into factor

mainlice$collection <- as.factor(mainlice$collection)
str(mainlice)

calmixmod <- glmer(all.cal ~ spp - 1 + (1|collection), 
                   data = mainlice, family = 'poisson')
lepmixmod <- glmer(all.leps ~ spp - 1 + (1|collection), 
                   data = mainlice, family = 'poisson')



calAMBDmod_0 <- glmmadmb(all.cal ~ 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

calAMBDmod_1 <- glmmadmb(all.cal ~ spp - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)


lepsAMBDmod_0 <- glmmadmb(all.leps ~ 1 + (1|collection),
                          data = mainlice, family = 'nbinom', zeroInflation = FALSE)
lepsAMBDmod_1 <- glmmadmb(all.leps ~ spp - 1 + (1|collection),
                          data = mainlice, family = 'nbinom', zeroInflation = FALSE)

anova(lepsAMBDmod_0, lepsAMBDmod_1)
anova(calAMBDmod_0,calAMBDmod_1)
AIC(lepsAMBDmod_0)
AIC(lepsAMBDmod_1)
AIC(calAMBDmod_0)
AIC(calAMBDmod_1)


summary(calmixmod)
summary(lepmixmod)
summary(calAMBDmod_1)
summary(lepsAMBDmod_1)
summary(calAMBDmod_0)
summary(lepsAMBDmod_0)


AICtab(calAMBDmod_0, calAMBDmod_1)
AICtab(lepsAMBDmod_0, lepsAMBDmod_1)



## Again, use the 'emmeans' package to get our estimates here

#turn the model into class 'emmGrid'
lepmod.spyr3.l <- as.list(ref_grid(lepmod.spyr3))
lepmod.spyr3.emm <- as.emmGrid(lepmod.spyr3.l)
#check model is in correct class
class(lepmod.spyr3.emm)
#use confint() to take the model of class 'emmGrid', and calculate confidence intervals
#IMPORTANT: ensure 'adjust = 'non'' because if 'adjust = 'tukey'', then the estimates
#will be wrong, use level = 0.95 for 95% confidence intervals
lepsCI3s <- confint(lepmod.spyr3.emm, adjust = 'none',level = 0.95)
lepsCI3s
#repeat for caligus
calmod.spyr3.l <- as.list(ref_grid(calmod.spyr3))
calmod.spyr3.emm <- as.emmGrid(calmod.spyr3.l)
class(calmod.spyr3.emm)
calCI3s <- confint(calmod.spyr3.emm, adjust = 'none', level = 0.95)
calCI3s

calmod.spyr3.l1 <- as.list(ref_grid(calmod.spyr31))
calmod.spyr3.emm1 <- as.emmGrid(calmod.spyr3.l1)
class(calmod.spyr3.emm1)
calCI3s1 <- confint(calmod.spyr3.emm1, adjust = 'none', level = 0.95)
calCI3s1







############################################## 10/20/18 #################################################

## make a couple of data frames to build the plots of lice abundance estimates
coeffs <- c(exp(-1.4080),exp(-0.8103), exp(-0.9425),  exp(-3.101), exp(-2.273), exp(-4.555))
ymin <- c(exp(-1.4080 - (2*0.1192)), exp(-0.8103 - (2*0.1011)), exp(-0.9425 - (2*0.0917)), 
          exp(-3.101 - (2*0.279)), exp(-2.273 - (2*0.242)), exp(-4.555 - (2*0.339)))
ymax <- c(exp(-1.4080 + (2*0.1192)), exp(-0.8103 + (2*0.1011)), exp(-0.9425 + (2*0.0917)), 
          exp(-3.101 + (2*0.279)), exp(-2.273 + (2*0.242)), exp(-4.555 + (2*0.339)))
spp <- c('CU', 'PI', 'SO', 'CU', 'PI', 'SO')
Lice <- c('C. clemensi','C. clemensi', 'C. clemensi', 'L. salmonis', 'L. salmonis', 'L. salmonis')
numest <- data.frame(coeffs, ymin, ymax, licesp, spp)


## make the plot 


numberoflicecomp <- ggplot(data = numest, aes(x = spp, y = coeffs, shape = Lice, colour = Lice)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, width = 0.6), 
                position = position_dodge(width = 0.8), colour = 'Black') +
  geom_line() +
  geom_point(size = 3, position = position_dodge(width = 0.8)) + 
  theme_classic() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = 'Species', y = 'Average Number of Motile Lice')


numberoflicecomp


####################################################### 11/03/18 ##########################################

## Make new dataframe only including the fish that have have lep and cal copepodite data 
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

coplice <- completeFun(mainlice, c("cal.cop", "lep.cop"))


##################################################### 11/05/18 ###########################################

##join up datasets, do up models with body size/environmental data and year (see email instructions) and then
## also do a chalimus vs copepodite model with salmon as fixed effect 

forklength <- read.csv("fork_length.csv")
sienedata <- read.csv("Seine_data.csv")
surveydata <- read.csv("Survey_data.csv")


summary(forklength$length.lab)
summary(forklength$length.field)

#since there's no measurements in the field length one, I'll just use the lab length one 


mainlice <- cbind(mainlice, forklength$length.lab)
mainlice <- dplyr::rename(mainlice, fork.length = "forklength$length.lab")
names(mainlice)



head(mainlice)

## Do up new models, first one will be chalimus versus copepodite model 

mainlicechal <- mainlice[!is.na(mainlice$chal.a),]

mainlicechal <- mainlicechal %>% 
  rowwise() %>% 
  mutate(all.chal = sum(c(chal.a, chal.b), na.rm = TRUE))
mainlicechal <- mainlicechal %>% 
  rowwise() %>% 
  mutate(all.cops = sum(c(all.leps, all.cal)))
mainlicechal$collection <- as.factor(mainlicechal$collection)

chalimusmod_0 <- glmmadmb(all.chal ~ 1 + (1|collection), 
                          data = mainlicechal, family = 'nbinom', zeroInflation = FALSE)
chalimusmod_1 <- glmmadmb(all.chal ~ spp - 1 + (1|collection), 
                          data = mainlicechal, family = 'nbinom', zeroInflation = FALSE)
copepmod_0 <- glmmadmb(all.cops ~ 1 + (1|collection), 
                       data = mainlicechal, family = 'nbinom', zeroInflation = FALSE)
copepmod_1 <- glmmadmb(all.cops ~ spp - 1 + (1|collection), 
                       data = mainlicechal, family = 'nbinom', zeroInflation = FALSE)

summary(chalimusmod_0)
summary(chalimusmod_1)
summary(copepmod_0)
summary(copepmod_1)

library(bbmle)
AICtab(chalimusmod_0, chalimusmod_1)
AICtab(copepmod_0, copepmod_1)


coeffs1 <- c(exp(0.408),exp(-2.799), exp(1.400),  exp(-1.494), exp(-0.569), exp(-0.895))
ymin1 <- c(exp(0.408 - (2*0.340)), exp(-2.799 - (2*1.033)), exp(1.400 - (2*0-.111)), 
           exp(-1.494 - (2*0.593)), exp(-0.569 - (2*0.458)), exp(-0.895 - (2*0.182)))
ymax1 <- c(exp(0.408 + (2*0.340)), exp(-2.799 + (2*1.033)), exp(1.400 + (2*0-.111)), 
           exp(-1.494 + (2*0.593)), exp(-0.569 + (2*0.458)), exp(-0.895 + (2*0.182)))
spp1 <- c('CU', 'PI', 'SO', 'CU', 'PI', 'SO')
Lice <- c('chalimus','chalimus', 'chalimus', 'copepodite', 'copepodite', 'copepodite')
numest2 <- data.frame(coeffs1, ymin1, ymax1, Lice, spp1)


## make the plot 


chalvscopp <- ggplot(data = numest2, aes(x = spp1, y = coeffs1, shape = Lice, colour = Lice)) +
  geom_errorbar(aes(ymin = ymin1, ymax = ymax1, width = 0.6), 
                position = position_dodge(width = 0.8), colour = 'Black') +
  geom_line() +
  geom_point(size = 3, position = position_dodge(width = 0.8)) + 
  theme_classic() +
  scale_y_continuous(limits = c(0,5)) +
  labs(x = 'Species', y = 'Average Number of Copepodite and Chalimus Lice')
chalvscopp



#this requires a dataframe with no NAs in the fork.length colum 
mainlicefork <- mainlice[!is.na(mainlice$fork.length),]
summary(mainlicefork$fork.length)

mainlicefork$collection <- as.factor(mainlicefork$collection)

calAMBDmod_2n <- glmmadmb(all.cal ~ 1 + (1|collection), 
                          data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)
calAMBDmod_2 <- glmmadmb(all.cal ~ spp + year + fork.length - 1 + (1|collection), 
                         data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)
lepsAMBDmod_2n <- glmmadmb(all.leps ~  1 + (1|collection),
                           data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)
lepsAMBDmod_2 <- glmmadmb(all.leps ~ spp + year + fork.length - 1 + (1|collection),
                          data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

summary(calAMBDmod_2)
summary(lepsAMBDmod_2)
summary(calAMBDmod_2n)
summary(lepsAMBDmod_2n)


library(bbmle)
AICtab(calAMBDmod_2n, calAMBDmod_2)
AICtab(lepsAMBDmod_2n, lepsAMBDmod_2)


################################# Full Model Set 
#fixed effects will be spp, fork.length, year, site.region
calmodnull <- glmmadmb(all.cal ~ 1 + (1|collection), 
                       data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)
#1
calmod.sp <- glmmadmb(all.cal ~ spp - 1 + (1|collection), 
                      data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)
#2
calmod.fl <- glmmadmb(all.cal ~ fork.length - 1 + (1|collection), 
                      data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)
#3
calmod.yr <- glmmadmb(all.cal ~ year - 1 + (1|collection), 
                      data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)
#4
calmod.sr <- glmmadmb(all.cal ~ site.region - 1 + (1|collection), 
                      data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

calmod.spfl <- glmmadmb(all.cal ~ spp + fork.length - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

calmod.spyr <- glmmadmb(all.cal ~ spp + year - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

calmod.spsr <- glmmadmb(all.cal ~ spp + site.region - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

calmod.flyr <- glmmadmb(all.cal ~ year + fork.length - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

calmod.flsr <- glmmadmb(all.cal ~ fork.length + site.region - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

calmod.yrsr <- glmmadmb(all.cal ~ year + site.region - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

calmod.spflyr <- glmmadmb(all.cal ~ spp + fork.length + year - 1 + (1|collection), 
                          data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

calmod.spflsr <- glmmadmb(all.cal ~ spp + fork.length + site.region - 1 + (1|collection), 
                          data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

calmod.spyrsr <- glmmadmb(all.cal ~ spp + year + site.region - 1 + (1|collection), 
                          data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

calmod.flyrsr <- glmmadmb(all.cal ~ fork.length + year + site.region - 1 + (1|collection), 
                          data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

calmod.spflyrsr <- glmmadmb(all.cal ~ spp + fork.length + year + site.region - 1 + (1|collection), 
                            data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

####leps
lepmodnull <- glmmadmb(all.leps ~ 1 + (1|collection), 
                       data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)
#1
lepmod.sp <- glmmadmb(all.leps ~ spp - 1 + (1|collection), 
                      data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)
#2
lepmod.fl <- glmmadmb(all.leps ~ fork.length - 1 + (1|collection), 
                      data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)
#3
lepmod.yr <- glmmadmb(all.leps ~ year - 1 + (1|collection), 
                      data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)
#4
lepmod.sr <- glmmadmb(all.leps ~ site.region - 1 + (1|collection), 
                      data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

lepmod.spfl <- glmmadmb(all.leps ~ spp + fork.length - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

lepmod.spyr <- glmmadmb(all.leps ~ spp + year - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

lepmod.spsr <- glmmadmb(all.leps ~ spp + site.region - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

lepmod.flyr <- glmmadmb(all.leps ~ fork.length + year - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

lepmod.flsr <- glmmadmb(all.leps ~ fork.length + site.region - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

lepmod.yrsr <- glmmadmb(all.leps ~ year + site.region - 1 + (1|collection), 
                        data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

lepmod.spflyr <- glmmadmb(all.leps ~ spp + fork.length + year - 1 + (1|collection), 
                          data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

lepmod.spflsr <- glmmadmb(all.leps ~ spp + fork.length + site.region - 1 + (1|collection), 
                          data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

lepmod.spyrsr <- glmmadmb(all.leps ~ spp + year + site.region - 1 + (1|collection), 
                          data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

lepmod.flyrsr <- glmmadmb(all.leps ~ year + fork.length + site.region - 1 + (1|collection), 
                          data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

lepmod.spflyrsr <- glmmadmb(all.leps ~ spp + fork.length + year + site.region - 1 + (1|collection), 
                            data = mainlicefork, family = 'nbinom', zeroInflation = FALSE)

###summaries
summary(lepmodnull)
summary(lepmod.sp); AICtab(lepmodnull, lepmod.sp)
summary(lepmod.fl); AICtab(lepmodnull, lepmod.fl)
summary(lepmod.yr); AICtab(lepmodnull, lepmod.yr)
summary(lepmod.sr); AICtab(lepmodnull, lepmod.sr)
summary(lepmod.spfl); AICtab(lepmodnull, lepmod.spfl)
summary(lepmod.spyr); AICtab(lepmodnull, lepmod.spyr)
summary(lepmod.spsr); AICtab(lepmodnull, lepmod.spsr)
summary(lepmod.flyr); AICtab(lepmodnull, lepmod.flyr)
summary(lepmod.flsr); AICtab(lepmodnull, lepmod.flsr)
summary(lepmod.yrsr); AICtab(lepmodnull, lepmod.yrsr)
summary(lepmod.spflyr); AICtab(lepmodnull, lepmod.spflyr)
summary(lepmod.spflsr); AICtab(lepmodnull, lepmod.spflsr)
summary(lepmod.spyrsr); AICtab(lepmodnull, lepmod.spyrsr)
summary(lepmod.flyrsr); AICtab(lepmodnull, lepmod.flyrsr)
summary(lepmod.spflyrsr); AICtab(lepmodnull, lepmod.spflyrsr)

library(bbmle)



summary(calmodnull)
summary(calmod.sp); AICtab(calmodnull, calmod.sp)
summary(calmod.fl); AICtab(calmodnull, calmod.fl)
summary(calmod.yr); AICtab(calmodnull, calmod.yr)
summary(calmod.sr); AICtab(calmodnull, calmod.sr)
summary(calmod.spfl); AICtab(calmodnull, calmod.spfl)
summary(calmod.spyr); AICtab(calmodnull, calmod.spyr)
summary(calmod.spsr); AICtab(calmodnull, calmod.spsr)
summary(calmod.flyr); AICtab(calmodnull, calmod.flyr)
summary(calmod.flsr); AICtab(calmodnull, calmod.flsr)
summary(calmod.yrsr); AICtab(calmodnull, calmod.yrsr)
summary(calmod.spflyr); AICtab(calmodnull, calmod.spflyr)
summary(calmod.spflsr); AICtab(calmodnull, calmod.spflsr)
summary(calmod.spyrsr); AICtab(calmodnull, calmod.spyrsr)
summary(calmod.flyrsr); AICtab(calmodnull, calmod.flyrsr)
summary(calmod.spflyrsr); AICtab(calmodnull, calmod.spflyrsr)


## Full model but with fork length taken out (so using mainlice dataset in full again then)

#### Good thing to note is that all the years are fairly well represented in the dataset, 
#### wasn't as much the case in the one with fork lengths


mainlice$collection <- as.factor(mainlice$collection)
mainlice$year <- as.factor(mainlice$year)
## models
calmodnull1 <- glmmadmb(all.cal ~ 1 + (1|collection), 
                        data = mainlice, family = 'nbinom', zeroInflation = FALSE)
#1
calmod.sp1 <- glmmadmb(all.cal ~ spp - 1 + (1|collection), 
                       data = mainlice, family = 'nbinom', zeroInflation = FALSE)
#2
calmod.yr1 <- glmmadmb(all.cal ~ year - 1 + (1|collection), 
                       data = mainlice, family = 'nbinom', zeroInflation = FALSE)
#3
calmod.sr1 <- glmmadmb(all.cal ~ site.region - 1 + (1|collection), 
                       data = mainlice, family = 'nbinom', zeroInflation = FALSE)

calmod.spyr1 <- glmmadmb(all.cal ~ spp + year - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

calmod.spsr1 <- glmmadmb(all.cal ~ spp + site.region - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

calmod.yrsr1 <- glmmadmb(all.cal ~ year + site.region - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

calmod.spyrsr1 <- glmmadmb(all.cal ~ spp + year + site.region - 1 + (1|collection), 
                           data = mainlice, family = 'nbinom', zeroInflation = FALSE)

####leps
lepmodnull1 <- glmmadmb(all.leps ~ 1 + (1|collection), 
                        data = mainlice, family = 'nbinom', zeroInflation = FALSE)
#1
lepmod.sp1 <- glmmadmb(all.leps ~ spp - 1 + (1|collection), 
                       data = mainlice, family = 'nbinom', zeroInflation = FALSE)
#2
lepmod.yr1 <- glmmadmb(all.leps ~ year - 1 + (1|collection), 
                       data = mainlice, family = 'nbinom', zeroInflation = FALSE)
#3
lepmod.sr1 <- glmmadmb(all.leps ~ site.region - 1 + (1|collection), 
                       data = mainlice, family = 'nbinom', zeroInflation = FALSE)

lepmod.spyr1 <- glmmadmb(all.leps ~ spp + year - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

lepmod.spsr1 <- glmmadmb(all.leps ~ spp + site.region - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

lepmod.yrsr1 <- glmmadmb(all.leps ~ year + site.region - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

lepmod.spyrsr1 <- glmmadmb(all.leps ~ spp + year + site.region - 1 + (1|collection), 
                           data = mainlice, family = 'nbinom', zeroInflation = FALSE)

###summaries
summary(lepmodnull1)
summary(lepmod.sp1) 
summary(lepmod.yr1) 
summary(lepmod.sr1)
summary(lepmod.spyr1) 
summary(lepmod.spsr1) 
summary(lepmod.yrsr1) 
summary(lepmod.spyrsr1) 

AICtab(lepmod.spyr1, lepmodnull1)
AICtab(lepmod.spyr1, lepmod.yr1)
AICtab(lepmod.spyr1, lepmod.sp1)
AICtab(lepmod.spyr1, lepmod.spsr1)
AICtab(lepmod.spyr1, lepmod.sr1)
AICtab(lepmod.spyr1, lepmod.yrsr1)
AICtab(lepmod.spyr1, lepmod.spyrsr1)

summary(calmodnull1)
summary(calmod.sp1)
summary(calmod.yr1)
summary(calmod.sr1)
summary(calmod.spyr1)
summary(calmod.spsr1)
summary(calmod.yrsr1)
summary(calmod.spyrsr1)

AICtab(calmod.spyrsr1, calmod.sp1)
AICtab(calmod.spyrsr1, calmod.yr1)
AICtab(calmod.spyrsr1, calmod.sr1)
AICtab(calmod.spyrsr1, calmod.spyr1)
AICtab(calmod.spyrsr1, calmod.spsr1)
AICtab(calmod.spyrsr1, calmod.yrsr1)
AICtab(calmod.spyrsr1, calmodnull1)

#for Cal model, the last model (spyrsr1) is the lowest AIC 
#so reverse transform coeffs, +/- 2 sds 

coeffscal1 <- c(exp(-1.161), exp(-0.581), exp(-0.686), exp(-0.579),
                exp(-0.778), exp(-0.516), exp(0.402))
minussdcal1 <- c(exp(-1.161 - (2*0.173)), exp(-0.581 - (2*0.158)), 
                 exp(-0.686 - (2*0.153)), exp(-0.579 - (2*0.158)),
                 exp(-0.778 - (2*0.285)), exp(-0.516 - (2*0.177)), 
                 exp(0.402 - (2*0.126)))
plussdcal1 <- c(exp(-1.161 + (2*0.173)), exp(-0.581 + (2*0.158)), 
                exp(-0.686 + (2*0.153)), exp(-0.579 + (2*0.158)),
                exp(-0.778 + (2*0.285)), exp(-0.516 + (2*0.177)), 
                exp(0.402 + (2*0.126)))

coeffscal1
minussdcal1
plussdcal1

#for Lep model, the model with just year and species is the lowest AIC 
#so reverse transform coeffs, +/- 2 sds 

coeffslep1 <- c(exp(-2.074), exp(-1.283), exp(-3.572), exp(-0.662),
                exp(-2.108), exp(-2.138))
minussdlep1 <- c(exp(-2.074 - (2*0.362)), exp(-1.283 - (2*0.324)), 
                 exp(-3.572 - (2*0.403)), exp(-0.662 - (2*0.392)),
                 exp(-2.108 - (2*0.900)), exp(-2.138 - (2*0.518)))
plussdlep1 <- c(exp(-2.074 + (2*0.362)), exp(-1.283 + (2*0.324)), 
                exp(-3.572 + (2*0.403)), exp(-0.662 + (2*0.392)),
                exp(-2.108 + (2*0.900)), exp(-2.138 + (2*0.518))) 
coeffslep1
minussdlep1
plussdlep1

############################ 12/12/18 ###################
# So for the next steps, try and plot some of what's going on using the final model 
# and understanding that things are additive here - then try and plot the results of the two final models
# Also run a model set with region*species, region*year and then a full one with 
# species*year + species*region + year*region 

calmod.spyr2 <- glmmadmb(all.cal ~ spp*year - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

calmod.spsr2 <- glmmadmb(all.cal ~ spp*site.region - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

calmod.yrsr2 <- glmmadmb(all.cal ~ year*site.region - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

calmod.spyrsr2 <- glmmadmb(all.cal ~ spp*year + spp*site.region + year*site.region - 1 + (1|collection), 
                           data = mainlice, family = 'nbinom', zeroInflation = FALSE)

####leps

lepmod.spyr2 <- glmmadmb(all.leps ~ spp*year - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

lepmod.spsr2 <- glmmadmb(all.leps ~ spp*site.region - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

lepmod.yrsr2 <- glmmadmb(all.leps ~ year*site.region - 1 + (1|collection), 
                         data = mainlice, family = 'nbinom', zeroInflation = FALSE)

lepmod.spyrsr2 <- glmmadmb(all.leps ~ spp*year + spp*site.region + year*site.region - 1 + (1|collection), 
                           data = mainlice, family = 'nbinom', zeroInflation = FALSE)

###summaries
summary(lepmod.spyr2) 
summary(lepmod.spsr2) 
summary(lepmod.yrsr2) 
summary(lepmod.spyrsr2) 

AICtab(lepmod.spyr1, lepmod.spyr2)
AICtab(lepmod.spyr1, lepmod.spsr2)
AICtab(lepmod.spyr1, lepmod.yrsr2)
AICtab(lepmod.spyr1, lepmod.spyrsr2)


summary(calmod.spyr2)
summary(calmod.spsr2)
summary(calmod.yrsr2)
summary(calmod.spyrsr2)

AICtab(calmod.spyrsr1, calmod.sp1)
AICtab(calmod.spyrsr1, calmod.yr1)
AICtab(calmod.spyrsr1, calmod.sr1)
AICtab(calmod.spyrsr1, calmod.spyr1)
AICtab(calmod.spyrsr1, calmod.spsr1)
AICtab(calmod.spyrsr1, calmod.yrsr1)
AICtab(calmod.spyrsr1, calmodnull1)
AICtab(calmod.spyrsr1, calmod.spyr2)
AICtab(calmod.spyrsr1, calmod.spsr2)
AICtab(calmod.spyrsr1, calmod.yrsr2)
AICtab(calmod.spyrsr1, calmod.spyrsr2)

################################## Plots with multiple effects 

## Okay so what I'm going to try is to use the emmeans package to get the confidence intervals 
## for my different levels

#turn the model into class 'emmGrid'
lepmod.spyr1.l <- as.list(ref_grid(lepmod.spyr1))
lepmod.spyr1.emm <- as.emmGrid(lepmod.spyr1.l)
#check model is in correct class
class(lepmod.spyr1.emm)
#use confint() to take the model of class 'emmGrid', and calculate confidence intervals
#IMPORTANT: ensure 'adjust = 'non'' because if 'adjust = 'tukey'', then the estimates
#will be wrong, use level = 0.95 for 95% confidence intervals
lepsCIs <- confint(lepmod.spyr1.emm, adjust = 'none',level = 0.95)
lepsCIs
#repeat for caligus
calmod.spyrsr1.l <- as.list(ref_grid(calmod.spyrsr1))
calmod.spyrsr1.emm <- as.emmGrid(calmod.spyrsr1.l)
class(calmod.spyrsr1.emm)
calCIs <- confint(calmod.spyrsr1.emm, adjust = 'none', level = 0.95)
calCIs

## now take the data from the above functions and put them into a bunch of objects/vectors so I can plot

calcoeffs2015D <- c(exp(-1.161), exp(-0.581), exp(-0.686))
calcoeffs2016D <- c(exp(-1.161+(-0.579)), exp(-0.581+(-0.579)), exp(-0.686+(-0.579)))
calcoeffs2017D <- c(exp(-1.161+(-0.778)), exp(-0.581+(-0.778)), exp(-0.686+(-0.778)))
calcoeffs2018D <- c(exp(-1.161+(-0.516)), exp(-0.581+(-0.516)), exp(-0.686+(-0.516)))
calcoeffs2015J <- c(exp(-1.161+(0.402)), exp(-0.581+(0.402)), exp(-0.686+(0.402)))
calcoeffs2016J <- c(exp(-1.161+(-0.579)+(0.402)), exp(-0.581+(-0.579)+(0.402)), exp(-0.686+(-0.579)+(0.402)))
calcoeffs2017J <- c(exp(-1.161+(-0.778)+(0.402)), exp(-0.581+(-0.778)+(0.402)), exp(-0.686+(-0.778)+(0.402)))
calcoeffs2018J <- c(exp(-1.161+(-0.516)+(0.402)), exp(-0.581+(-0.516)+(0.402)), exp(-0.686+(-0.516)+(0.402)))

upcalCIs2015D <- c(exp(-0.8222145), exp(-0.2702810), exp(-0.3872934))
upcalCIs2016D <- c(exp(-1.2754004), exp(-0.7150025), exp(-0.8098796))
upcalCIs2017D <- c(exp(-1.2670000), exp(-0.7138563), exp(-0.7907213))
upcalCIs2018D <- c(exp(-1.1759210), exp(-0.6212535), exp(-0.7381814))
upcalCIs2015J <- c(exp(-0.3494470), exp(0.2069287), exp(0.0857904))
upcalCIs2016J <- c(exp(-0.8274216), exp(-0.2658039), exp(-0.3664922))
upcalCIs2017J <- c(exp(-0.8267532), exp(-0.2725516), exp(-0.3543889))
upcalCIs2018J <- c(exp(-0.7294362), exp(-0.1731821), exp(-0.2937770))

localCIs2015D <- c(exp(-1.5001661), exp(-0.8907664), exp(-0.9855528))
localCIs2016D <- c(exp(-2.2040038), exp(-1.6030685), exp(-1.7199901))
localCIs2017D <- c(exp(-2.6109362), exp(-2.0027467), exp(-2.1376804))
localCIs2018D <- c(exp(-2.1781186), exp(-1.5714528), exp(-1.6663236))
localCIs2015J <- c(exp(-1.1698886), exp(-0.5649311), exp(-0.6555915))
localCIs2016J <- c(exp(-1.8489376), exp(-1.2492221), exp(-1.3603325))
localCIs2017J <- c(exp(-2.2481379), exp(-1.6410063), exp(-1.7709678))
localCIs2018J <- c(exp(-1.8215583), exp(-1.2164792), exp(-1.3076830))

lepscoeffs2015 <- c(exp(-2.074), exp(-1.283), exp(-3.572))
lepscoeffs2016 <- c(exp(-2.074+(-0.662)), exp(-1.283+(-0.662)), exp(-3.572+(-0.662)))
lepscoeffs2017 <- c(exp(-2.074+(-2.108)), exp(-1.283+(-2.108)), exp(-3.572+(-2.108)))
lepscoeffs2018 <- c(exp(-2.074+(-2.138)), exp(-1.283+(-2.138)), exp(-3.572+(-2.138)))

uplepsCIs2015 <- c(exp(-1.364259), exp(-0.648040), exp(-2.781702))
uplepsCIs2016 <- c(exp(-1.792393), exp(-1.031745), exp(-3.159972))
uplepsCIs2017 <- c(exp(-2.208059), exp(-1.486955), exp(-3.620113))
uplepsCIs2018 <- c(exp(-2.816773), exp(-2.084027), exp(-4.350522))

lolepsCIs2015 <- c(exp(-2.783116), exp(-1.918645), exp(-4.363001))
lolepsCIs2016 <- c(exp(-3.679171), exp(-2.859129), exp(-5.308920))
lolepsCIs2017 <- c(exp(-6.154670), exp(-5.295085), exp(-7.739946))
lolepsCIs2018 <- c(exp(-5.606849), exp(-4.758905), exp(-7.070429))

sal <- c("CU","PI","SO")
regD <- c('Discovery Islands', 'Discovery Islands', 'Discovery Islands')
regJ <- c('Johnstone Strait', 'Johnstone Strait', 'Johnstone Strait')
yr2015 <- c('2015', '2015','2015')
yr2016 <- c('2016', '2016','2016')
yr2017 <- c('2017', '2017','2017')
yr2018 <- c('2018', '2018','2018')
regNA <- c('NA', 'NA', 'NA')
licespC <- c('cal','cal', 'cal')
licespL <- c('lep', 'lep', 'lep')

cal2015D <- rename(data.frame(calcoeffs2015D,upcalCIs2015D,localCIs2015D,sal,regD,yr2015,licespC), 
                   coeffs = calcoeffs2015D, upCI = upcalCIs2015D, loCI = localCIs2015D,reg = regD, yr = yr2015, lice = licespC)
cal2016D <- rename(data.frame(calcoeffs2016D,upcalCIs2016D,localCIs2016D,sal,regD,yr2016,licespC), 
                   coeffs = calcoeffs2016D, upCI = upcalCIs2016D, loCI = localCIs2016D,reg = regD, yr = yr2016, lice = licespC) 
cal2017D <- rename(data.frame(calcoeffs2017D,upcalCIs2017D,localCIs2017D,sal,regD,yr2017,licespC), 
                   coeffs = calcoeffs2017D, upCI = upcalCIs2017D, loCI = localCIs2017D,reg = regD, yr = yr2017, lice = licespC)
cal2018D <- rename(data.frame(calcoeffs2018D,upcalCIs2018D,localCIs2018D,sal,regD,yr2018,licespC), 
                   coeffs = calcoeffs2018D, upCI = upcalCIs2018D, loCI = localCIs2018D,reg = regD, yr = yr2018, lice = licespC)
cal2015J <- rename(data.frame(calcoeffs2015J,upcalCIs2015J,localCIs2015J,sal,regJ,yr2015,licespC), 
                   coeffs = calcoeffs2015J, upCI = upcalCIs2015J, loCI = localCIs2015J,reg = regJ, yr = yr2015, lice = licespC)
cal2016J <- rename(data.frame(calcoeffs2016J,upcalCIs2016J,localCIs2016J,sal,regJ,yr2016,licespC), 
                   coeffs = calcoeffs2016J, upCI = upcalCIs2016J, loCI = localCIs2016J,reg = regJ, yr = yr2016, lice = licespC)
cal2017J <- rename(data.frame(calcoeffs2017J,upcalCIs2017J,localCIs2017J,sal,regJ,yr2017,licespC), 
                   coeffs = calcoeffs2017J, upCI = upcalCIs2017J, loCI = localCIs2017J,reg = regJ, yr = yr2017, lice = licespC)
cal2018J <- rename(data.frame(calcoeffs2018J,upcalCIs2018J,localCIs2018J,sal,regJ,yr2018,licespC), 
                   coeffs = calcoeffs2018J, upCI = upcalCIs2018J, loCI = localCIs2018J,reg = regJ, yr = yr2018, lice = licespC)

calall <- rbind(cal2015D,cal2016D,cal2017D,cal2018D,cal2015J,cal2016J,cal2017J,cal2018J)

lep2015 <- rename(data.frame(lepscoeffs2015,uplepsCIs2015,lolepsCIs2015,sal,regNA,yr2015,licespL), 
                  coeffs = lepscoeffs2015, upCI = uplepsCIs2015, loCI = lolepsCIs2015, reg = regNA, yr = yr2015,lice = licespL)
lep2016 <- rename(data.frame(lepscoeffs2016,uplepsCIs2016,lolepsCIs2016,sal,regNA,yr2016,licespL), 
                  coeffs = lepscoeffs2016, upCI = uplepsCIs2016, loCI = lolepsCIs2016, reg = regNA, yr = yr2016,lice = licespL) 
lep2017 <- rename(data.frame(lepscoeffs2017,uplepsCIs2017,lolepsCIs2017,sal,regNA,yr2017,licespL), 
                  coeffs = lepscoeffs2017, upCI = uplepsCIs2017, loCI = lolepsCIs2017, reg = regNA, yr = yr2017,lice = licespL)
lep2018 <- rename(data.frame(lepscoeffs2018,uplepsCIs2018,lolepsCIs2018,sal,regNA,yr2018,licespL), 
                  coeffs = lepscoeffs2018, upCI = uplepsCIs2018, loCI = lolepsCIs2018, reg = regNA, yr = yr2018,lice = licespL)

lepsall <- rbind(lep2015,lep2016,lep2017,lep2018)

leg_title <- 'Salmon Species'
library(ggplot2)
lepsfullmodplot <- lepsall %>% 
  group_by(., yr,sal) %>% 
  ggplot(aes(x = sal, y = coeffs, colour = sal)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI,width = 0), colour = 'Black')+
  geom_point(size = 4) +
  facet_wrap(~yr,nrow=1,strip.position = "bottom")+
  theme(strip.background = element_blank(), strip.placement = "outside") + 
  scale_color_manual(leg_title,values=c('seagreen2', 'hotpink1', 'steelblue2'))+
  labs(title = "L. salmonis Effects Plot", x = 'Salmon Species/Year', y = 'Average Number of Motile Lice Per Fish')
lepsfullmodplot

calfullmodplot <- calall %>% 
  group_by(., yr,sal,reg) %>% 
  ggplot(aes(x = sal, y = coeffs, colour = sal, shape = reg)) +
  scale_shape_manual(values = c(15,17)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI,width = 0), position = position_dodge(width = 0.8),colour = 'Black')+
  geom_point(size = 4,position = position_dodge(width = 0.8)) +
  facet_wrap(~yr,nrow=1,strip.position = "bottom")+
  theme(strip.background = element_blank(), strip.placement = "outside") + 
  scale_color_manual(leg_title,values=c('seagreen2', 'hotpink1', 'steelblue2'))+
  labs(title = "C. clemensi Effects Plot", x = 'Salmon Species/Year', y = 'Average Number of Motile Lice Per Fish',shape = 'Region')+
  guides(shape = guide_legend(override.aes = list(shape = c(0,2)), type = 'b'))
calfullmodplot

####### We're going to try putting all of the lice (both cals and leps) on the same plot and see how that works 
##so could either do this by calling two data sources in the differnt ggplot commands or by putting it all in one dataframe
##lets try putting it all in one dataframe first, I think that will be easier

#new dataframe
allliceforplot <- rbind(lep2015,lep2016,lep2017,lep2018,cal2015D,cal2016D,cal2017D,cal2018D,cal2015J,cal2016J,cal2017J,cal2018J)
#new plot
bothfullmodplot <- allliceforplot %>% 
  group_by(., yr,sal,reg,lice) %>% 
  ggplot(aes(x=sal, y = coeffs, colour = sal, shape = reg)) +
  scale_shape_manual(values = c(9,15,17),labels = c('L. salmonis', 'C. clemensi, Discovery Islands', 'C. clemensi, Johnstone Strait'))+
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~yr, nrow=1, strip.position = 'bottom')+
  theme(strip.background = element_blank(),strip.placement = 'outside') +
  scale_color_manual(leg_title,values=c('seagreen2', 'hotpink1', 'steelblue2'))+
  labs(title = "Lice Effects Plot", x = 'Salmon Species/Year', y = 'Average Number of Motile Lice Per Fish',shape = 'Lice and Region')+
  guides(shape = guide_legend(override.aes = list(shape = c(9,0,2)), type = 'b'))
bothfullmodplot

## Poisson Distribution Models w/ year as fixed effect

chumrmod.calp <- glmmadmb(all.cal ~ site.region + year - 1 + (1|week), 
                          data = chum.region, family = 'poisson', zeroInflation = FALSE)
chumrmod.lepsp <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                           data = chum.region, family = 'poisson', zeroInflation = FALSE)
pinkrmod.calp <- glmmadmb(all.cal ~ site.region + year - 1  + (1|week), 
                          data = pink.region, family = 'poisson', zeroInflation = FALSE)
pinkrmod.lepsp <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                           data = pink.region, family = 'poisson', zeroInflation = FALSE)
sockrmod.calp <- glmmadmb(all.cal ~ site.region + year - 1 + (1|week), 
                          data = sock.region, family = 'poisson', zeroInflation = FALSE)
sockrmod.lepsp <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                           data = sock.region, family = 'poisson', zeroInflation = FALSE)

chumrmod.calp1 <- glmmadmb(all.cal ~ site.region + year - 1 + (1|week), 
                           data = chum.region, family = 'poisson', zeroInflation = TRUE)
chumrmod.lepsp1 <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                            data = chum.region, family = 'poisson', zeroInflation = TRUE)
pinkrmod.calp1 <- glmmadmb(all.cal ~ site.region + year - 1  + (1|week), 
                           data = pink.region, family = 'poisson', zeroInflation = TRUE)
pinkrmod.lepsp1 <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                            data = pink.region, family = 'poisson', zeroInflation = TRUE)
sockrmod.calp1 <- glmmadmb(all.cal ~ site.region + year - 1 + (1|week), 
                           data = sock.region, family = 'poisson', zeroInflation = TRUE)
sockrmod.lepsp1 <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                            data = sock.region, family = 'poisson', zeroInflation = TRUE)

summary(chumrmod.calp)
summary(chumrmod.lepsp)
summary(pinkrmod.calp)
summary(pinkrmod.lepsp)
summary(sockrmod.calp)
summary(sockrmod.lepsp)

#poisson
chumrmod.cal.lp <- as.list(ref_grid(chumrmod.calp))
chumrmod.cal.emmp <- as.emmGrid(chumrmod.cal.lp)
chumrcalCIp <- confint(chumrmod.cal.emmp, adjust = 'none', level = 0.95)

chumrmod.leps.lp <- as.list(ref_grid(chumrmod.lepsp))
chumrmod.leps.emmp <- as.emmGrid(chumrmod.leps.lp)
chumrlepsCIp <- confint(chumrmod.leps.emmp, adjust = 'none', level = 0.95)

pinkrmod.cal.lp <- as.list(ref_grid(pinkrmod.calp))
pinkrmod.cal.emmp <- as.emmGrid(pinkrmod.cal.lp)
pinkrcalCIp <- confint(pinkrmod.cal.emmp, adjust = 'none', level = 0.95)

pinkrmod.leps.lp <- as.list(ref_grid(pinkrmod.lepsp))
pinkrmod.leps.emmp <- as.emmGrid(pinkrmod.leps.lp)
pinkrlepsCIp <- confint(pinkrmod.leps.emmp, adjust = 'none', level = 0.95)

sockrmod.cal.lp <- as.list(ref_grid(sockrmod.calp))
sockrmod.cal.emmp <- as.emmGrid(sockrmod.cal.lp)
sockrcalCIp <- confint(sockrmod.cal.emmp, adjust = 'none', level = 0.95)

sockrmod.leps.lp <- as.list(ref_grid(sockrmod.lepsp))
sockrmod.leps.emmp <- as.emmGrid(sockrmod.leps.lp)
sockrlepsCIp <- confint(sockrmod.leps.emmp, adjust = 'none', level = 0.95)

#make a df to store the estimated values in
chumrcalCIp #taking the values from this object and pasting them into the brackets below
chumcalcoeffsp <- c(exp(-0.8537077),exp(-0.6314040),exp(-1.3667853),exp(-1.1444815),
                    exp(-1.1214345),exp(-0.8991308),exp(-1.2416695),exp(-1.0193658))
chumcalupCIp <- c(exp(-0.4918984),exp(-0.2539737), exp(-0.8786619), exp(-0.6471530),
                  exp(-0.6741818), exp(-0.4313183), exp(-0.7842832), exp(-0.5571315))
chumcalloCIp <- c(exp(-1.215517), exp(-1.008834), exp(-1.854909), exp(-1.641810),
                  exp(-1.568687), exp(-1.36694), exp(-1.699056), exp(-1.481600))
sal1 <- c('Chum', 'Chum','Chum', 'Chum','Chum', 'Chum','Chum', 'Chum')
year <- c('2015', '2015', '2016', '2016', '2017', '2017', '2018', '2018')
Region<- c('DI', 'JS','DI', 'JS','DI', 'JS','DI', 'JS')
lice1 <- c('C. clemensi', 'C. clemensi','C. clemensi', 'C. clemensi',
           'C. clemensi', 'C. clemensi','C. clemensi', 'C. clemensi')
chumsregioncalp <- rename(data.frame(chumcalcoeffsp,chumcalupCIp,chumcalloCIp,sal1,year,Region,lice1),
                          coeffs = chumcalcoeffsp, upCI = chumcalupCIp, loCI = chumcalloCIp,
                          lice = lice1, sal = sal1)

chumrlepsCIp
chumlepcoeffsp <- c(exp(-2.234709),exp(-3.235173),exp(-2.210428),exp(-3.210891),
                    exp(-3.516624),exp(-4.517087),exp(-3.381688),exp(-4.382151))
chumlepupCIp <- c(exp(-1.858122),exp(-2.620665), exp(-1.633015), exp(-2.416651),
                  exp(-2.578996), exp(-3.454700), exp(-2.438212), exp(-3.344734))
chumleploCIp <-c(exp(-2.611297), exp(-3.849680), exp(-2.787840), exp(-4.005131),
                 exp(-4.454251), exp(-5.579474), exp(-4.325164), exp(-5.419569))
sal1 <- c('Chum', 'Chum','Chum', 'Chum','Chum', 'Chum','Chum', 'Chum')
year <- c('2015', '2015', '2016', '2016', '2017', '2017', '2018', '2018')
Region<- c('DI', 'JS','DI', 'JS','DI', 'JS','DI', 'JS')
lice2 <- c('L. salmonis', 'L. salmonis', 'L. salmonis', 'L. salmonis',
           'L. salmonis', 'L. salmonis', 'L. salmonis', 'L. salmonis')
chumsregionlepp <- rename(data.frame(chumlepcoeffsp,chumlepupCIp,chumleploCIp,sal1,year,Region,lice2),
                          coeffs = chumlepcoeffsp, upCI = chumlepupCIp, loCI = chumleploCIp,
                          sal = sal1, lice = lice2)

pinkrcalCIp               
pinkcalcoeffsp <- c(exp(-0.7262930), exp(-0.1996632), exp(-1.1277654),exp(-0.6011355),
                    exp(-0.7940383),exp(-0.2674084),exp(-0.8806231),exp(-0.3539932))
pinkcalupCIp <-c(exp(-0.3820842),exp(0.1403318), exp(-0.6618278), exp(-0.1482391),
                 exp(-0.3283793), exp(0.1882883), exp(-0.4682861), exp(0.0470466))
pinkcalloCIp <- c(exp(-1.0705019), exp(-0.5396581), exp(-1.5937029), exp(-1.0540319),
                  exp(-1.2596972), exp(-0.7231052), exp(-1.2929601), exp(-0.7550331))
sal2 <- c('Pink','Pink','Pink','Pink','Pink','Pink','Pink','Pink')
pinkregioncalp <- rename(data.frame(pinkcalcoeffsp,pinkcalupCIp,pinkcalloCIp,sal2,year,Region,lice1),
                         coeffs = pinkcalcoeffsp, upCI = pinkcalupCIp, loCI = pinkcalloCIp,
                         sal = sal2, lice = lice1)

pinkrlepsCIp 
pinklepscoeffsp <- c(exp(-0.734765), exp(-1.400695), exp(-1.545609),exp(-2.211539),
                     exp(-2.419929),exp(-3.085860),exp(-3.058858),exp(-3.724788))
pinklepsupCIp <-c(exp(-0.4353017),exp(-1.0561334), exp(-1.1220501), exp(-1.7089638),
                  exp(-1.5710980), exp(-2.1790181), exp(-2.2664969), exp(-2.9453005))
pinklepsloCIp <- c(exp(-1.034228), exp(-1.745257), exp(-1.969167), exp(-2.714114),
                   exp(-3.268761), exp(-3.992701), exp(-3.851218), exp(-4.504276))
pinkregionlepsp <- rename(data.frame(pinklepscoeffsp,pinklepsupCIp,pinklepsloCIp,sal2,year,Region,lice2),
                          coeffs = pinklepscoeffsp, upCI = pinklepsupCIp, loCI =  pinklepsloCIp,
                          sal = sal2, lice = lice2)
sockrcalCInb
sockrcalCIp
sockcalcoeffsp <- c(exp(-0.5701516), exp(-0.4069235), exp(-1.0811269),exp(-0.9178989),
                    exp(-1.3695692),exp(-1.2063411),exp(-1.6633525),exp(-1.5001244))
sockcalupCIp <-c(exp(-0.3825889),exp(-0.2208739), exp(-0.8580100), exp(-0.6971893),
                 exp(-1.1046470), exp(-0.9418721), exp(-1.2737680), exp(-1.1102294))
sockcalloCIp <- c(exp(-0.7577143), exp(-0.5929731), exp(-1.3042438), exp(-1.1386084),
                  exp(-1.6344913), exp(-1.4708101), exp(-2.0529370), exp(-1.8900195))
sal3 <- c('Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye')
sockregioncalp <- rename(data.frame(sockcalcoeffsp,sockcalupCIp,sockcalloCIp,sal3,year,Region,lice1),
                         coeffs = sockcalcoeffsp, upCI = sockcalupCIp, loCI = sockcalloCIp,
                         sal = sal3, lice = lice1)

sockrlepsCIp
socklepscoeffsp <- c(exp(-4.376264), exp(-3.620631), exp(-4.804352),exp(-4.048719),
                     exp(-4.895525),exp(-4.139893),exp(-4.798660),exp(-4.043027))
socklepsupCIp <-c(exp(-3.699253),exp(-3.018099), exp(-3.961351), exp(-3.267934),
                  exp(-3.855026), exp(-3.122455), exp(-3.439925), exp(-2.696673))
socklepsloCIp <- c(exp(-5.053274), exp(-4.223163), exp(-5.647353), exp(-4.829505),
                   exp(-5.936025), exp(-5.157331), exp(-6.157395), exp(-5.389381))
sal3 <- c('Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye')
sockregionlepsp <- rename(data.frame(socklepscoeffsp,socklepsupCIp,socklepsloCIp,sal3,year,Region,lice2),
                          coeffs = socklepscoeffsp, upCI = socklepsupCIp, loCI = socklepsloCIp,
                          sal = sal3, lice = lice2)

regionestbyspchump <- rbind(chumsregioncalp, chumsregionlepp)
regionestbysppinkp <- rbind(pinkregioncalp,pinkregionlepsp)
regionestbyspsockp <- rbind(sockregioncalp,sockregionlepsp)

#make some themes so it all looks nice beside each other
fte_themeSock <- function(){
  color.background = 'white'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill=color.background,color = color.background)) +
    theme(plot.background = element_rect(fill = color.background, color = color.background)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 15, vjust = 1.25)) +
    theme(axis.text.x = element_text(size = 10, color = color.axis.text)) + 
    theme(axis.text.y = element_text(size = 10, color = color.axis.text)) + 
    theme(axis.title.x = element_text(size = 12, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 12, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(color="black", size = 0.15)) +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(size = 10))+
    theme(legend.position = 'none')
} 
fte_themePink <- function(){
  color.background = 'white'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill=color.background,color = color.background)) +
    theme(plot.background = element_rect(fill = color.background, color = color.background)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 15, vjust = 1.25)) +
    theme(axis.text.x = element_text(size = 10, color = color.axis.text)) + 
    theme(axis.text.y = element_blank()) + 
    theme(axis.title.x = element_text(size = 12, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 12, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(linetype = 'dashed', color = 'grey75')) +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(size = 10))+
    theme(legend.position = 'none')
} 
fte_themeChum <- function(){
  color.background = 'white'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill=color.background,color = color.background)) +
    theme(plot.background = element_rect(fill = color.background, color = color.background)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 15, vjust = 1.25)) +
    theme(axis.text.x = element_text(size = 10, color = color.axis.text)) + 
    theme(axis.text.y = element_blank()) + 
    theme(axis.title.x = element_text(size = 12, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 12, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(linetype = 'dashed', color = 'grey75')) +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(size = 10))+
    theme(legend.position = c(0.8,0.82),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10))
} 

#make the plots
regionestbyspplotchump <- regionestbyspchump %>% 
  group_by(.,lice,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region, shape = lice)) +
  scale_shape_manual(values = c(17,15),labels = c('C. clemensi','L. salmonis'))+
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Chum Salmon",x = '', y = NULL,shape = 'Lice Species')+
  guides(shape = guide_legend(override.aes = list(shape = c(24,22)), type = 'b'))+
  scale_y_continuous(limits = c(0,1.25)) +
  fte_themeChum()
regionestbyspplotchump
regionestbyspplotpinkp <- regionestbysppinkp %>% 
  group_by(.,lice,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region, shape = lice)) +
  scale_shape_manual(values = c(17,15),labels = c('C. clemensi','L. salmonis'))+
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  theme(strip.background = element_blank(),strip.placement = 'outside') +
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Pink Salmon", x = 'Year & Region', y = NULL,shape = 'Lice Species')+
  guides(shape = guide_legend(override.aes = list(shape = c(24,22)), type = 'b'))+
  scale_y_continuous(limits = c(0,1.25)) +
  fte_themePink()
regionestbyspplotpinkp
regionestbyspplotsockp <- regionestbyspsockp %>% 
  group_by(.,lice,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region, shape = lice)) +
  scale_shape_manual(values = c(17,15),labels = c('C. clemensi','L. salmonis'))+
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  theme(strip.background = element_blank(),strip.placement = 'outside') +
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Sockeye Salmon", x = '', y = 'Average Number of Motile C. clemensi Lice Per Fish',shape = 'Lice Species')+
  guides(shape = guide_legend(override.aes = list(shape = c(24,22)), type = 'b'))+
  scale_y_continuous(limits = c(0,1.25))+
  fte_themeSock()
regionestbyspplotsockp

regionbothp <- plot_grid(regionestbyspplotsockp,regionestbyspplotpinkp,regionestbyspplotchump,
                         nrow = 1, labels = NULL, rel_widths = c(1,1,1)) 
regionbothp

## Now repeat but just by species of lice

#subset the dataframes
regionestbyspchumLEPp <- regionestbyspchump %>% 
  filter(lice == 'L. salmonis')
regionestbyspchumCALp <- regionestbyspchump %>% 
  filter(lice == 'C. clemensi')
regionestbysppinkLEPp <- regionestbysppinkp %>% 
  filter(lice == 'L. salmonis')
regionestbysppinkCALp <- regionestbysppinkp %>% 
  filter(lice == 'C. clemensi')
regionestbyspsockLEPp <- regionestbyspsockp %>% 
  filter(lice == 'L. salmonis')
regionestbyspsockCALp <- regionestbyspsockp %>% 
  filter(lice == 'C. clemensi')

#do lep plots
regionestbyspplotchumLEPp <- regionestbyspchumLEPp %>% 
  group_by(.,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Chum Salmon", x = '', y = NULL, Region = c('Discovery Islands', 'Johnstone Strait'))+
  scale_y_continuous(limits = c(0,0.65)) +
  fte_themeChum()
regionestbyspplotchumLEPp
regionestbyspplotpinkLEPp <- regionestbysppinkLEPp %>% 
  group_by(.,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Pink Salmon", x = 'Year and Region', y = NULL)+
  scale_y_continuous(limits = c(0,0.65)) +
  fte_themePink()
regionestbyspplotpinkLEPp
regionestbyspplotsockLEPp <- regionestbyspsockLEPp %>%
  group_by(.,Region,year) %>%
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  theme(strip.background = element_blank(),strip.placement = 'outside') +
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Sockeye Salmon", x = '', y = 'Average Number of Motile L. salmonis Lice Per Fish')+
  scale_y_continuous(limits = c(0,0.65))+
  fte_themeSock()
regionestbyspplotsockLEPp

### option 1
regionLEPp <- plot_grid(regionestbyspplotsockLEPp,regionestbyspplotpinkLEPp,regionestbyspplotchumLEPp, 
                        nrow = 1, labels = NULL) 
regionLEPp #print 800x450 ---- this option has the sockeye column there even though there's no data

#do Cal plots
regionestbyspplotchumCALp <- regionestbyspchumCALp %>% 
  group_by(.,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Chum Salmon", x = '', y = NULL, Region = c('Discovery Islands', 'Johnstone Strait'))+
  scale_y_continuous(limits = c(0,1.25)) +
  fte_themeChum()
regionestbyspplotchumCALp
regionestbyspplotpinkCALp <- regionestbysppinkCALp %>% 
  group_by(.,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Pink Salmon", x = 'Year and Region', y = NULL)+
  scale_y_continuous(limits = c(0,1.25)) +
  fte_themePink()
regionestbyspplotpinkCALp
regionestbyspplotsockCALp <- regionestbyspsockCALp %>% 
  group_by(.,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  theme(strip.background = element_blank(),strip.placement = 'outside') +
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Sockeye Salmon", x = '', y = 'Average Number of Motile C. clemensi Lice Per Fish')+
  scale_y_continuous(limits = c(0,1.25))+
  fte_themeSock()
regionestbyspplotsockCALp

regionCALp <- plot_grid(regionestbyspplotsockCALp,regionestbyspplotpinkCALp,regionestbyspplotchumCALp,
                        nrow = 1, labels = NULL) 
regionCALp #export 800x450













summary(chumrmod.calp1)
summary(chumrmod.calp)


chumrmod.calp1 <- glmmadmb(all.cal ~ site.region + year - 1 + (1|week), 
                           data = chum.region, family = 'poisson', zeroInflation = TRUE)
chumrmod.lepsp1 <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                            data = chum.region, family = 'poisson', zeroInflation = TRUE)
pinkrmod.calp1 <- glmmadmb(all.cal ~ site.region + year - 1  + (1|week), 
                           data = pink.region, family = 'poisson', zeroInflation = TRUE)
pinkrmod.lepsp1 <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                            data = pink.region, family = 'poisson', zeroInflation = TRUE)
sockrmod.calp1 <- glmmadmb(all.cal ~ site.region + year - 1 + (1|week), 
                           data = sock.region, family = 'poisson', zeroInflation = TRUE)
sockrmod.lepsp1 <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                            data = sock.region, family = 'poisson', zeroInflation = TRUE)

chumrmod.cal.lp1 <- as.list(ref_grid(chumrmod.calp1))
chumrmod.cal.emmp1 <- as.emmGrid(chumrmod.cal.lp1)
chumrcalCIp1 <- confint(chumrmod.cal.emmp1, adjust = 'none', level = 0.95)

chumrmod.leps.lp1 <- as.list(ref_grid(chumrmod.lepsp1))
chumrmod.leps.emmp1 <- as.emmGrid(chumrmod.leps.lp1)
chumrlepsCIp1 <- confint(chumrmod.leps.emmp1, adjust = 'none', level = 0.95)

pinkrmod.cal.lp1 <- as.list(ref_grid(pinkrmod.calp1))
pinkrmod.cal.emmp1 <- as.emmGrid(pinkrmod.cal.lp1)
pinkrcalCIp1 <- confint(pinkrmod.cal.emmp1, adjust = 'none', level = 0.95)

pinkrmod.leps.lp1 <- as.list(ref_grid(pinkrmod.lepsp1))
pinkrmod.leps.emmp1 <- as.emmGrid(pinkrmod.leps.lp1)
pinkrlepsCIp1 <- confint(pinkrmod.leps.emmp1, adjust = 'none', level = 0.95)

sockrmod.cal.lp1 <- as.list(ref_grid(sockrmod.calp1))
sockrmod.cal.emmp1 <- as.emmGrid(sockrmod.cal.lp1)
sockrcalCIp1 <- confint(sockrmod.cal.emmp1, adjust = 'none', level = 0.95)

sockrmod.leps.lp1 <- as.list(ref_grid(sockrmod.lepsp1))
sockrmod.leps.emmp1 <- as.emmGrid(sockrmod.leps.lp1)
sockrlepsCIp1 <- confint(sockrmod.leps.emmp1, adjust = 'none', level = 0.95)

#make a df to store the estimated values in
chumrcalCIp1 #taking the values from this object and pasting them into the brackets below
chumcalcoeffsp1 <- c(exp(-0.5518259),exp(-0.3136274),exp(-1.0682002),exp(-0.8300017),
                     exp(-0.8309125),exp(-0.5927140),exp(-0.9386216),exp(-0.7004231))
chumcalupCIp1 <- c(exp(-0.4918984),exp(-0.2539737), exp(-0.8786619), exp(-0.6471530),
                   exp(-0.6741818), exp(-0.4313183), exp(-0.7842832), exp(-0.5571315))
chumcalloCIp1 <- c(exp(-0.9381152), exp(-0.7190851), exp(-1.6491294), exp(-1.5279982),
                   exp(-1.2947895), exp(-1.0733781), exp(-1.3961352), exp(-1.1829157))
sal1 <- c('Chum', 'Chum','Chum', 'Chum','Chum', 'Chum','Chum', 'Chum')
year <- c('2015', '2015', '2016', '2016', '2017', '2017', '2018', '2018')
Region<- c('DI', 'JS','DI', 'JS','DI', 'JS','DI', 'JS')
lice1 <- c('C. clemensi', 'C. clemensi','C. clemensi', 'C. clemensi',
           'C. clemensi', 'C. clemensi','C. clemensi', 'C. clemensi')
chumsregioncalp1 <- rename(data.frame(chumcalcoeffsp1,chumcalupCIp1,chumcalloCIp1,sal1,year,Region,lice1),
                           coeffs = chumcalcoeffsp1, upCI = chumcalupCIp1, loCI = chumcalloCIp1,
                           lice = lice1, sal = sal1)

chumrlepsCIp1
chumlepcoeffsp1 <- c(exp(-0.0914707),exp(-1.0877587),exp(0.7811800),exp(-0.2151080),
                     exp(-1.1229987),exp(-2.1192867),exp(-1.0109312),exp(-2.0072192))
chumlepupCIp1 <- c(exp(0.4077517),exp(-0.2963448), exp(1.8324355), exp(1.0436184),
                   exp(-0.1776994), exp(-1.0551078), exp(0.0207984), exp(-0.7263421))
chumleploCIp1 <-c(exp(-0.590693), exp(-1.879173), exp(-0.270075), exp(-1.473834),
                  exp(-2.068298), exp(-3.183466), exp(-2.042661), exp(-3.288096))
sal1 <- c('Chum', 'Chum','Chum', 'Chum','Chum', 'Chum','Chum', 'Chum')
year <- c('2015', '2015', '2016', '2016', '2017', '2017', '2018', '2018')
Region<- c('DI', 'JS','DI', 'JS','DI', 'JS','DI', 'JS')
lice2 <- c('L. salmonis', 'L. salmonis', 'L. salmonis', 'L. salmonis',
           'L. salmonis', 'L. salmonis', 'L. salmonis', 'L. salmonis')
chumsregionlepp1 <- rename(data.frame(chumlepcoeffsp1,chumlepupCIp1,chumleploCIp1,sal1,year,Region,lice2),
                           coeffs = chumlepcoeffsp1, upCI = chumlepupCIp1, loCI = chumleploCIp1,
                           sal = sal1, lice = lice2)

pinkrcalCIp1              
pinkcalcoeffsp1 <- c(exp(-0.5757091), exp(-0.0473006), exp(-1.0011845),exp(-0.4727759),
                     exp(-0.6797360),exp(-0.1513274),exp(-0.7466859),exp(-0.2182774))
pinkcalupCIp1 <-c(exp(-0.2064323),exp(0.3166059), exp(-0.4399752), exp(0.1676111),
                  exp(-0.1974211), exp(0.3426733), exp(-0.3247723), exp(0.2177233))
pinkcalloCIp1 <- c(exp(-0.9449860), exp(-0.4112071), exp(-1.5623937), exp(-1.1131628),
                   exp(-1.1620509), exp(-0.6453282), exp(-1.1685995), exp(-0.6542780))
sal2 <- c('Pink','Pink','Pink','Pink','Pink','Pink','Pink','Pink')
pinkregioncalp1 <- rename(data.frame(pinkcalcoeffsp1,pinkcalupCIp1,pinkcalloCIp1,sal2,year,Region,lice1),
                          coeffs = pinkcalcoeffsp1, upCI = pinkcalupCIp1, loCI = pinkcalloCIp1,
                          sal = sal2, lice = lice1)

pinkrlepsCIp1 
pinklepscoeffsp1 <- c(exp(0.4034494), exp(-0.1703766), exp(-0.2828551),exp(-0.8566811),
                      exp(-1.1976583),exp(-1.7714844),exp(-1.8122096),exp(-2.3860356))
pinklepsupCIp1 <-c(exp(0.7322138),exp(0.2545632), exp(0.438184), exp(-0.0803998),
                   exp(-0.2193540), exp(-0.8461572), exp(-1.1014681), exp(-1.6487699))
pinklepsloCIp1 <- c(exp(0.0746850), exp(-0.5953164), exp(-1.0038950), exp(-1.6329625),
                    exp(-2.1759626), exp(-2.6968115), exp(-2.5229510), exp(-3.1233013))
pinkregionlepsp1 <- rename(data.frame(pinklepscoeffsp1,pinklepsupCIp1,pinklepsloCIp1,sal2,year,Region,lice2),
                           coeffs = pinklepscoeffsp1, upCI = pinklepsupCIp1, loCI =  pinklepsloCIp1,
                           sal = sal2, lice = lice2)
sockrcalCIp1
sockcalcoeffsp1 <- c(exp(-0.0903600), exp(0.0741443), exp(-0.5704924),exp(-0.4059880),
                     exp(-0.8678823),exp(-0.7033780),exp(-1.1692455),exp(-1.0047412))
sockcalupCIp1 <-c(exp(0.1031708),exp(0.2657798), exp(-0.3064371), exp(-0.0979652),
                  exp(-0.5854371), exp(-0.4223600), exp(-0.7536440), exp(-0.5875000))
sockcalloCIp1 <- c(exp(-0.2838908), exp(-0.1174911), exp(-0.8345476), exp(-0.7140108),
                   exp(-1.1503276), exp(-0.9843961), exp(-1.5848471), exp(-1.4219824))
sal3 <- c('Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye')
sockregioncalp1 <- rename(data.frame(sockcalcoeffsp1,sockcalupCIp1,sockcalloCIp1,sal3,year,Region,lice1),
                          coeffs = sockcalcoeffsp1, upCI = sockcalupCIp1, loCI = sockcalloCIp1,
                          sal = sal3, lice = lice1)

sockrlepsCIp1
socklepscoeffsp1 <- c(exp(-2.382638), exp(-1.701401), exp(-2.724824),exp(-2.043586),
                      exp(-2.857852),exp(-2.176615),exp(-2.794616),exp(-2.113378))
socklepsupCIp1 <-c(exp(-3.699253),exp(-3.018099), exp(-3.961351), exp(-3.267934),
                   exp(-3.855026), exp(-3.122455), exp(-3.439925), exp(-2.696673))
socklepsloCIp1 <- c(exp(-3.297216), exp(-2.513355), exp(-4.070959), exp(-3.336331),
                    exp(-4.150083), exp(-3.375871), exp(-4.378524), exp(-3.629943))
sal3 <- c('Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye','Sockeye')
sockregionlepsp1 <- rename(data.frame(socklepscoeffsp1,socklepsupCIp1,socklepsloCIp1,sal3,year,Region,lice2),
                           coeffs = socklepscoeffsp1, upCI = socklepsupCIp1, loCI = socklepsloCIp1,
                           sal = sal3, lice = lice2)

regionestbyspchump1 <- rbind(chumsregioncalp1, chumsregionlepp1)
regionestbysppinkp1 <- rbind(pinkregioncalp1,pinkregionlepsp1)
regionestbyspsockp1 <- rbind(sockregioncalp1,sockregionlepsp1)

#make some themes so it all looks nice beside each other
fte_themeSock <- function(){
  color.background = 'white'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill=color.background,color = color.background)) +
    theme(plot.background = element_rect(fill = color.background, color = color.background)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 15, vjust = 1.25)) +
    theme(axis.text.x = element_text(size = 10, color = color.axis.text)) + 
    theme(axis.text.y = element_text(size = 10, color = color.axis.text)) + 
    theme(axis.title.x = element_text(size = 12, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 12, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(color="black", size = 0.15)) +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(size = 10))+
    theme(legend.position = 'none')
} 
fte_themePink <- function(){
  color.background = 'white'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill=color.background,color = color.background)) +
    theme(plot.background = element_rect(fill = color.background, color = color.background)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 15, vjust = 1.25)) +
    theme(axis.text.x = element_text(size = 10, color = color.axis.text)) + 
    theme(axis.text.y = element_blank()) + 
    theme(axis.title.x = element_text(size = 12, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 12, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(linetype = 'dashed', color = 'grey75')) +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(size = 10))+
    theme(legend.position = 'none')
} 
fte_themeChum <- function(){
  color.background = 'white'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill=color.background,color = color.background)) +
    theme(plot.background = element_rect(fill = color.background, color = color.background)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 15, vjust = 1.25)) +
    theme(axis.text.x = element_text(size = 10, color = color.axis.text)) + 
    theme(axis.text.y = element_blank()) + 
    theme(axis.title.x = element_text(size = 12, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 12, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(linetype = 'dashed', color = 'grey75')) +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(size = 10))+
    theme(legend.position = c(0.8,0.82),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10))
} 

#make the plots
regionestbyspplotchump1 <- regionestbyspchump1 %>% 
  group_by(.,lice,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region, shape = lice)) +
  scale_shape_manual(values = c(17,15),labels = c('C. clemensi','L. salmonis'))+
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Chum Salmon",x = '', y = NULL,shape = 'Lice Species')+
  guides(shape = guide_legend(override.aes = list(shape = c(24,22)), type = 'b'))+
  scale_y_continuous(limits = c(0,3.5)) +
  fte_themeChum()
regionestbyspplotchump1
regionestbyspplotpinkp1 <- regionestbysppinkp1 %>% 
  group_by(.,lice,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region, shape = lice)) +
  scale_shape_manual(values = c(17,15),labels = c('C. clemensi','L. salmonis'))+
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  theme(strip.background = element_blank(),strip.placement = 'outside') +
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Pink Salmon", x = 'Year & Region', y = NULL,shape = 'Lice Species')+
  guides(shape = guide_legend(override.aes = list(shape = c(24,22)), type = 'b'))+
  scale_y_continuous(limits = c(0,3.5)) +
  fte_themePink()
regionestbyspplotpinkp1
regionestbyspplotsockp1 <- regionestbyspsockp1 %>% 
  group_by(.,lice,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region, shape = lice)) +
  scale_shape_manual(values = c(17,15),labels = c('C. clemensi','L. salmonis'))+
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  theme(strip.background = element_blank(),strip.placement = 'outside') +
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Sockeye Salmon", x = '', y = 'Average Number of Motile C. clemensi Lice Per Fish',shape = 'Lice Species')+
  guides(shape = guide_legend(override.aes = list(shape = c(24,22)), type = 'b'))+
  scale_y_continuous(limits = c(0,3.5))+
  fte_themeSock()
regionestbyspplotsockp1

regionbothp1 <- plot_grid(regionestbyspplotsockp1,regionestbyspplotpinkp1,regionestbyspplotchump1,
                          nrow = 1, labels = NULL, rel_widths = c(1,1,1)) 
regionbothp1

## Now repeat but just by species of lice

#subset the dataframes
regionestbyspchumLEPp1 <- regionestbyspchump1 %>% 
  filter(lice == 'L. salmonis')
regionestbyspchumCALp1 <- regionestbyspchump1 %>% 
  filter(lice == 'C. clemensi')
regionestbysppinkLEPp1 <- regionestbysppinkp1 %>% 
  filter(lice == 'L. salmonis')
regionestbysppinkCALp1 <- regionestbysppinkp1 %>% 
  filter(lice == 'C. clemensi')
regionestbyspsockLEPp1 <- regionestbyspsockp1 %>% 
  filter(lice == 'L. salmonis')
regionestbyspsockCALp1 <- regionestbyspsockp1 %>% 
  filter(lice == 'C. clemensi')

#do lep plots
regionestbyspplotchumLEPp1 <- regionestbyspchumLEPp1 %>% 
  group_by(.,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Chum Salmon", x = '', y = NULL, Region = c('Discovery Islands', 'Johnstone Strait'))+
  scale_y_continuous(limits = c(0,0.65)) +
  fte_themeChum()
regionestbyspplotchumLEPp1
regionestbyspplotpinkLEPp1 <- regionestbysppinkLEPp1 %>% 
  group_by(.,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Pink Salmon", x = 'Year and Region', y = NULL)+
  scale_y_continuous(limits = c(0,0.65)) +
  fte_themePink()
regionestbyspplotpinkLEPp1
regionestbyspplotsockLEPp1 <- regionestbyspsockLEPp1 %>%
  group_by(.,Region,year) %>%
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  theme(strip.background = element_blank(),strip.placement = 'outside') +
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Sockeye Salmon", x = '', y = 'Average Number of Motile L. salmonis Lice Per Fish')+
  scale_y_continuous(limits = c(0,0.65))+
  fte_themeSock()
regionestbyspplotsockLEPp1

### option 1
regionLEPp1 <- plot_grid(regionestbyspplotsockLEPp1,regionestbyspplotpinkLEPp1,regionestbyspplotchumLEPp1, 
                         nrow = 1, labels = NULL) 
regionLEPp1 #print 800x450 ---- this option has the sockeye column there even though there's no data

#do Cal plots
regionestbyspplotchumCALp <- regionestbyspchumCALp %>% 
  group_by(.,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Chum Salmon", x = '', y = NULL, Region = c('Discovery Islands', 'Johnstone Strait'))+
  scale_y_continuous(limits = c(0,1.25)) +
  fte_themeChum()
regionestbyspplotchumCALp
regionestbyspplotpinkCALp <- regionestbysppinkCALp %>% 
  group_by(.,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Pink Salmon", x = 'Year and Region', y = NULL)+
  scale_y_continuous(limits = c(0,1.25)) +
  fte_themePink()
regionestbyspplotpinkCALp
regionestbyspplotsockCALp <- regionestbyspsockCALp %>% 
  group_by(.,Region,year) %>% 
  ggplot(aes(x=Region, y = coeffs, colour = Region)) +
  geom_errorbar(aes(ymin=loCI, ymax = upCI, width = 0), position = position_dodge(width = 0.8),colour = 'Black') +
  geom_point(size = 4, position = position_dodge(width=0.8)) +
  facet_wrap(~year, nrow=1, strip.position = 'bottom')+
  theme(strip.background = element_blank(),strip.placement = 'outside') +
  scale_color_manual(values=c('grey60', 'grey20'))+
  labs(title = "Sockeye Salmon", x = '', y = 'Average Number of Motile C. clemensi Lice Per Fish')+
  scale_y_continuous(limits = c(0,1.25))+
  fte_themeSock()
regionestbyspplotsockCALp

regionCALp <- plot_grid(regionestbyspplotsockCALp,regionestbyspplotpinkCALp,regionestbyspplotchumCALp,
                        nrow = 1, labels = NULL) 
regionCALp 








