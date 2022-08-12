##############################
# Decision tree models
# 12.08.2022
# Steven R. Talbot
##############################

# Please install these packages before using the script!
library(ggplot2)
library(dplyr)
library(readxl)
library(reshape2)
library(rpart)
library(rpart.plot)
library(rattle)

# Load some data ----------------------------------------------------------
figpath    <- "./figs/"

# Meta data from patients
path       <- ".../DATA_2003-2020_2022-06-17.xlsx"     # path to your data
metadat    <- data.frame(read_excel(path, sheet = 1))

# other data (eGFR)
# adjust the path to your "other" data
cache      <- data.frame(read_excel(".../DWH_Datenlieferung_2022-06-18_FULL_RESHAPED_NO-OUTLIERS-CSA.xlsx", sheet = 1))
df         <- cache[, c("mhhid", "Tag", "Ciclosporin1_Wert", "eGFR1_Wert")]



# Preprocessing -----------------------------------------------------------
# please adjust the rules to your liking (<=90 etc)
# also make sure that all parameters that should be included are within the loop
fulldat <- NULL
dat     <- NULL
for(i in 1:length(unique(metadat$mhhid))){
  ddf <- NULL
  ddf <- NULL
  dat     <- metadat[metadat$mhhid %in% unique(metadat$mhhid)[i], ][1,]
  ddf     <- df[df$mhhid %in% unique(metadat$mhhid)[i], ]
  
  ddf$ckd         <- dat$ckd
  ddf$rrt_postsct <- dat$rrt_postsct
  ddf$akinmax     <- dat$akin_max
  ddf$icu1        <- dat$icu1
  ddf$surv_status <- dat$surv_status
  
  ddf$lab         <- ifelse(ddf$Tag <= 90, "pre","post")
  
  ddf$pre_counts  <- sum(ddf$Tag <= 90 &  ddf$eGFR1_Wert <60, na.rm=TRUE) 
  ddf$post_counts <- sum(ddf$Tag   >90 &  ddf$eGFR1_Wert <60, na.rm=TRUE) 
  
  fulldat <- rbind(fulldat, ddf)
}
head(fulldat)



# Working decision tree ---------------------------------------------------

# data aggregation (i.e., median eGFR1_Wert)
dat <- fulldat %>% 
  group_by(mhhid, ckd, rrt_postsct, icu1,surv_status, akinmax, lab) %>% 
  summarise(median=median(eGFR1_Wert, na.rm=TRUE)) %>% 
  as.data.frame()

# dummy code  lab (split in pre and post data)
sdat <- tidyr::spread(dat, key = lab, value = median)

# prepare factors for dc
#sdat$icu1        <- factor(sdat$icu1)
sdat$rrt_postsct <- factor(sdat$rrt_postsct)
sdat$akinmax[is.na(sdat$akinmax )] <- 0
sdat$akinmax     <- factor(sdat$akinmax)



# The decision tree -------------------------------------------------------
# with 10-fold cross-validation (xval=10)
control     <- rpart.control( 
  minbucket = round(5 / 3),
  minsplit  = 3,
  maxdepth  = 6,
  xval      = 10,
  cp        = 0.01) # complexity parameter

fit <- rpart(ckd ~  pre +  akinmax  + post + rrt_postsct,
             data    = sdat[,!(names(sdat) %in% c( "mhhid" ))],
             method  = 'class',
             control = control )

summary(fit)
printcp(fit) # display the results 
plotcp(fit)  # visualize cross-validation results 
 


# plotting the tree as jpeg
jpeg(filename=paste(figpath, "Decision Tree",".jpg",sep=""), pointsize = 60,
     quality = 400, width = 2200, height = 1800)
#rpart.plot(fit , type= 2, extra=106, digits = 4)
fancyRpartPlot(fit, palettes=c("Greys", "Oranges"), type=1, digits = 3)
dev.off()
# Each node box displays the classification, the probability of each class at 
# that node (i.e. the probability of the class conditioned on the node) and the
# percentage of observations used at that node. 

# do some predictions with the full data!
predict_unseen       <- predict(fit, sdat, type = 'class')

table_mat <- table(sdat$ckd, predict_unseen)
table_mat

caret::confusionMatrix(factor(sdat$ckd), predict_unseen )




# Repeat with leave-one-out validation ------------------------------------
control     <- rpart.control( 
  minbucket = round(5 / 3),
  minsplit  = 3,
  maxdepth  = 4,
  xval      = 10,
  cp        = 0.01) # complexity parameter


patients <- unique(sdat$mhhid)
res      <- NULL
for(i in 1:length(patients)){
  
  traindat <- sdat[!(sdat$mhhid %in% patients[i]), ]
  testdat  <- sdat[sdat$mhhid   %in% patients[i]   , ]
  
  fit      <- rpart(ckd ~  pre +  akinmax  + post + rrt_postsct, ,
                    data    = traindat[,!(names(traindat) %in% c( "mhhid" ))],
                    method  = 'class',
                    control = control )
 
  pred     <- predict(fit, testdat, type = 'class')
  
  res <- rbind(res, data.frame(mhhid     = patients[i],
                               true      = testdat$ckd, 
                               predicted = as.numeric(as.character(pred)),
                               sum       = testdat$ckd + as.numeric(as.character(pred))) )

}

d   <- table(res$true, res$predicted)
acc <- sum(diag(d))/sum(d)
acc


# Again with expanded CV grid ---------------------------------------------
set.seed(123)
cp_grid    <- expand.grid(.cp=seq(0.01,0.5,0.01))

clean      <- subset(sdat, select = -c(mhhid, post) )
clean      <- clean[complete.cases(clean), ]
clean$ckd  <- factor(clean$ckd)

fit      <-  caret::train(ckd ~  pre +  akinmax  +  rrt_postsct,  
                  data     = clean[,!(names(clean) %in% c( "mhhid", "post"  ))],
                  method   = 'rpart',
                  tuneGrid = cp_grid,
                  control  = control )
summary(fit$finalModel)
 

  
# test the variants with splitting the data in training and test sets
# however, I think that there are not enough cases for doing that!











 
 
 



