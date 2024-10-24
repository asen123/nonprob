
## Author : Aditi Sen
## Date : Oct 2024

############################################################################################################
# Topic : IPW, bias-corrected and composite estimates for unknown bias case,
#         using modeling on ABS 3
############################################################################################################

library(caret);library(dplyr);library(lme4);library(ggplot2);
library(usmap);library(tidyverse);library(maps);library(fastDummies)
library(PNWColors);library(gganimate);library(ggeasy);library(ggpubr);
library(patchwork);library(boot);library(rpart);library(rpart.plot);
library(partykit);library(gbm);library(pROC); library(gridExtra);
library(randomForest);library(neuralnet); library(nnet); library(gbm);
#library(fastAdaboost)

o123_p123 <- dat %>% dplyr::select(CASEID,SAMPID,WEIGHT,
                                                       # x variables
                                                       AGECAT3BM, RACETHN, EDUC3,
                                                       GENDER, REGION,
                                                       # y variables
                                                       bm_INSURANCE, bm_BLOODPR, bm_HHPARENT, 
                                                       bm_FOODALLERGY, bm_JOBLASTYR, bm_RETACCT, 
                                                       bm_UNEMPCOMP, bm_WRKRSCOMP, bm_FOODSTAMPS, 
                                                       bm_SOCSEC, bm_UNION, bm_CITIZEN) %>% 
  # remove refused values
  filter(across(c(AGECAT3BM, RACETHN, EDUC3, GENDER, REGION,
                  bm_INSURANCE, bm_BLOODPR, bm_HHPARENT, 
                  bm_FOODALLERGY, bm_JOBLASTYR, bm_RETACCT, 
                  bm_UNEMPCOMP, bm_WRKRSCOMP, bm_FOODSTAMPS, 
                  bm_SOCSEC, bm_UNION, bm_CITIZEN), ~ . != 99)) %>%
  # change 2 to 0
  mutate(across(starts_with("bm"),~recode(.,"2"=0)))

# pivot longer for modelling
o123_p123_long <- o123_p123 %>%
  pivot_longer(cols = starts_with("bm_"),names_to = "qs", names_prefix = "bm_", 
    values_to = "y_val", values_drop_na = TRUE)  %>% 
  # convert to factor
  mutate_at(c('qs','AGECAT3BM','RACETHN','EDUC3','GENDER','REGION'),as.factor)

# convert response to yes no (for requirement in ML models)
o123_p123_long$y_new <- as.factor(ifelse(o123_p123_long$y_val==1,"yes","no"))

# Data preparation: Subset P3 and train-test split (80-20)
p3_data <- o123_p123_long %>% subset(SAMPID == 3)
set.seed(7345)
train <- sample(1:nrow(p3_data), 0.8*nrow(p3_data))
p3_train <- p3_data[train,]; p3_test <- p3_data[-train,]
# check response % in train and test
table(p3_train$AGECAT3BM); table(p3_test$AGECAT3BM)

############  Logistic Regression ############ 
logit_fit <- glm(y_new ~  qs + AGECAT3BM + RACETHN + EDUC3 + GENDER + REGION, 
                 family = binomial("logit"), data = p3_train)
summary(logit_fit)
# predict responses in test
pred_logit <- ifelse(predict(logit_fit, newdata = p3_test,type = "response")>=0.5,"yes","no")

############ Decision Tree ##################  
dt_fit <- rpart(y_new ~  qs + AGECAT3BM + RACETHN + EDUC3 + GENDER + REGION, 
                data = p3_train, method = "class")
rpart.plot(dt_fit, extra = 106)
print(dt_fit)
# predict responses in test
pred_dt <- ifelse(predict(dt_fit, newdata = p3_test, type="class")==1,"yes","no")

############ Random Forest ##################   
rf_fit <- train(y_new ~  qs + AGECAT3BM + RACETHN + EDUC3 + GENDER + REGION, 
            data = p3_train, method = 'rf',
            trControl = trainControl(method = 'cv', # Use cross-validation
                                      number = 5)) # Use 5 folds for cross-validation
rf_fit
# predict responses in test
pred_rf <- predict(rf_fit, newdata = p3_test, type="raw")

############ Gradient boosting ##################  
# For Gradient Boosting as implemented by the `gbm` package, we have to take care of 
# a number of tuning parameters. Now the `expand.grid` is helpful as it creates an object 
# with all possible combinations of our try-out values.

# In order to build a set of prediction models it is helpful to follow the `caret` workflow 
# and first decide how to conduct model tuning. Here we use 5-Fold Cross-Validation, 
# mainly to keep computation time to a minimum. `caret` offers many performance metrics, 
# however, they are stored in different functions that need to be combined first.

# For classification and regression using packages gbm and plyr with tuning parameters:
# Number of Boosting Iterations (n.trees, numeric)
# Max Tree Depth (interaction.depth, numeric)
# Shrinkage (shrinkage, numeric)
# Min. Terminal Node Size (n.minobsinnode, numeric)

evalStats <- function(...) c(twoClassSummary(...), defaultSummary(...), mnLogLoss(...))
grid <- expand.grid(interaction.depth = 1:3,n.trees = c(500, 750, 1000), 
                    shrinkage = c(0.05, 0.01),n.minobsinnode = 10)
# List the tuning grid...
grid
#Now we can specify the `trainControl` object.
ctrl  <- trainControl(method = "cv", number = 5, summaryFunction = evalStats, 
                      classProbs = TRUE)
#...and begin the tuning process.
set.seed(744)
gbm <- train(y_new ~  qs + AGECAT3BM + RACETHN + EDUC3 + GENDER + REGION, 
             data = p3_train, method = "gbm",trControl = ctrl,tuneGrid = grid, 
             metric = "logLoss",distribution = "bernoulli",verbose = FALSE)
plot(gbm)
summary(gbm)
# We can also extract single trees of the GBM ensemble.
pretty.gbm.tree(gbm$finalModel, i.tree = 1)
pretty.gbm.tree(gbm$finalModel, i.tree = 2)
# A quick look at feature importance.
plot(varImp(gbm), top = 15)
# predict responses in test
pred_gbm <- predict(gbm, newdata = p3_test, type="raw")

# Given predicted class membership, we can use the function `postResample` in order to 
# get a short summary of each models' performance in the test set.
postResample(pred = pred_logit, obs = p3_test$y_new)
postResample(pred = pred_dt, obs = p3_test$y_new)
postResample(pred = pred_rf, obs = p3_test$y_new)
postResample(pred = pred_gbm, obs = p3_test$y_new) #best 0.844 accuracy

# Confusion matrix: Accuracy, Sensitivity on test data
confusionMatrix(table(pred_logit, p3_test$y_new))
confusionMatrix(table(pred_dt, p3_test$y_new))
confusionMatrix(table(pred_rf, p3_test$y_new))
confusionMatrix(table(pred_gbm, p3_test$y_new))

# Creating `ROC` objects based on predicted probabilities...
p_dt <- predict(dt_fit, newdata = p3_test, type = "prob")
dt_roc <- roc(p3_test$y_new, p_dt[,2])
p_rf <- predict(rf_fit, newdata = p3_test, type = "prob")
rf_roc <- roc(p3_test$y_new, p_rf$yes)
p_gbm <- predict(gbm, newdata = p3_test, type = "prob")
gbm_roc <- roc(p3_test$y_new, p_gbm$yes)

auc(dt_roc)
auc(rf_roc)
auc(gbm_roc)

#Plotting the ROC curves.
ggroc(list(RF = rf_roc, GBM = gbm_roc)) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
               color="darkgrey", linetype="dashed") + 
  theme_bw() + coord_fixed() + labs(color = "Model")

#+ ggtitle("ROC curve") 

# get benchmark values
bchmk_age_long <- bchmk_age %>% dplyr::select(AGECAT,bm_INSURANCE_1,bm_BLOODPR_1,
                                                bm_HHPARENT_1,bm_FOODALLERGY_1,
                                                bm_JOBLASTYR_1,bm_RETACCT_1,
                                                bm_UNEMPCOMP_1,bm_WRKRSCOMP_1,
                                                bm_FOODSTAMPS_1,bm_SOCSEC_1,
                                                bm_UNION_1,bm_CITIZEN_1) %>% 
  #rename(AGECAT3BM = AGECAT)  %>% 
  rename_at(.vars = vars(ends_with("_1")),.funs = funs(sub("[_]1$", "", .))) %>% 
  pivot_longer(cols = starts_with("bm_"),names_to = "qs", names_prefix = "bm_", 
    values_to = "bnch_val", values_drop_na = TRUE)  

# function for obtaining predicted values and calculating MAE
pred_func <- function(sample_id){
  # subset data by sample id
  my_data <- o123_p123_long %>% subset(SAMPID == sample_id)
  # predict responses from different models
  my_data$pred_prob_lr <- predict(logit_fit, newdata = my_data,type = "response")
  my_data$pred_prob_dt <- predict(dt_fit, newdata = my_data, type = "prob")[,2]
  my_data$pred_prob_rf <- predict(rf_fit, newdata = my_data, type = "prob")[,2]
  my_data$pred_prob_gbm <- predict(gbm, newdata = my_data, type = "prob")[,2]
  # get the IPW weight
  data_wt <- filter(data_KML, SAMPID == sample_id) %>% dplyr::select(CASEID)
  if(sample_id==4){
    data_wt$np_weight <- np1_weight$IPW_weight
  }else if(sample_id==5){
    data_wt$np_weight <- np2_weight$IPW_weight
  }else if(sample_id==6){
    data_wt$np_weight <- np3_weight$IPW_weight
  }
  # join IPW weight
  my_data_new <- my_data %>% left_join(data_wt, by=c('CASEID'))
  my_data_new$AGECAT <- recode_factor(my_data_new$AGECAT3BM, 
                                      `1` = "18-29", `2` = "30-64", `3` = "65+")
  # summarize at age and qs level
  my_data_sum <- my_data_new %>%        
    group_by(AGECAT,qs) %>%
    summarise(logistic = weighted.mean(pred_prob_lr, np_weight), 
              decision_tree = weighted.mean(as.numeric(pred_prob_dt), np_weight),
              random_forest = weighted.mean(as.numeric(pred_prob_rf), np_weight),
              gradient_boost = weighted.mean(pred_prob_gbm, np_weight),
              act_mean = weighted.mean(y_val, np_weight))
  
  # join benchmark data
  my_data_sum1 <- bchmk_age_long %>% left_join(my_data_sum, by=c('AGECAT','qs'))
  
  print(paste0("MAE Logistic: ", round(mean(abs(my_data_sum1$logistic - 
                                                  my_data_sum1$bnch_val)),4)))
  print(paste0("MAE Decision Tree: ", round(mean(abs(my_data_sum1$decision_tree - 
                                                       my_data_sum1$bnch_val)),4)))
  print(paste0("MAE Random Forest: ", round(mean(abs(my_data_sum1$random_forest - 
                                                       my_data_sum1$bnch_val)),4)))
  print(paste0("MAE Gradient Boost: ", round(mean(abs(my_data_sum1$gradient_boost - 
                                                        my_data_sum1$bnch_val)),4)))
  print(paste0("MAE CLW: ", round(mean(abs(my_data_sum1$act_mean - 
                                             my_data_sum1$bnch_val)),4)))
  return(my_data_sum1)
}
# run function on 3 nonprob samples
pred_func(4) # O_1
pred_func(5) # O_2
pred_func(6) # O_3

write.csv(ml_est1,"/Users/aditisen/Downloads/KML work/ml_summary.csv")

##################################################################################
# After prediction model is selected use IPW on O3 with predicted responses
##################################################################################
# pivot longer o3 for prediction
o3_data <- data_KML %>% subset(SAMPID == 6) %>% remove_var_label() 
o3_long <- o3_data %>% pivot_longer(cols = starts_with("bm_"), 
                                    names_to = "qs", names_prefix = "bm_", 
                                    values_to = "y_val", values_drop_na = TRUE) %>% remove_var_label() 
# predict on o3  
o3_long$pred_prob_gbm <- predict(gbm, newdata = o3_long, type = "prob")[,2]
# o3_long$y_pred <- ifelse(o3_long$pred_prob_gbm > 0.5,1,0)
# table(o3_long$y_pred); table(o3_long$y_val) # refused??

# pivot wider back
o3_short <- o3_long %>% pivot_wider(names_from = qs, values_from = pred_prob_gbm)
o3_short <- rename(o3_short, union_mem =  UNION)
o3_short1 <- sqldf("select distinct CASEID, GENDER, AGECAT3BM, RACETHN, REGION, EDUC3,
                   sum(INSURANCE) as bm_INSURANCE, sum(BLOODPR) as bm_BLOODPR, 
                   sum(HHPARENT) as bm_HHPARENT, sum(FOODALLERGY) as bm_FOODALLERGY, 
                   sum(JOBLASTYR) as bm_JOBLASTYR, sum(RETACCT) as bm_RETACCT, 
                   sum(UNEMPCOMP) as bm_UNEMPCOMP, sum(WRKRSCOMP) as bm_WRKRSCOMP, 
                   sum(FOODSTAMPS) as bm_FOODSTAMPS, sum(SOCSEC) as bm_SOCSEC, 
                   sum(union_mem) as bm_UNION, sum(CITIZEN) as bm_CITIZEN 
                   from o3_short group by 1,2,3,4,5,6")
# use survey package to get estimates for p3
data_prob3 <- filter(dat_new, SAMPID == 3) 
data_prob3$new_weight <- data_prob3$WEIGHT*(N_B/nrow(data_prob3))
design_p3 <- survey::svydesign(ids = ~ CASEID, data = data_prob3, 
                               weights = ~ new_weight) 
bm_p1 <- survey::svymean(x = ~ 
                           bm_INSURANCE_1 + bm_BLOODPR_1 + bm_HHPARENT_1 + bm_FOODALLERGY_1 +
                           bm_JOBLASTYR_1 + bm_RETACCT_1 + bm_UNEMPCOMP_1 + bm_WRKRSCOMP_1 + 
                           bm_FOODSTAMPS_1 + bm_SOCSEC_1 + bm_UNION_1 + bm_CITIZEN_1,
                         design = design_p3)
p3_output <- as.data.frame(bm_p1)
# IPW function (with predicted propabilities as target)
set.seed(12325672)
IPW_o3 <- nonprob(selection = ~ AGECAT3BM + GENDER + RACETHN + EDUC3 + REGION, 
                     target = ~ 
                       bm_INSURANCE + bm_BLOODPR + bm_HHPARENT + bm_FOODALLERGY + 
                       bm_JOBLASTYR + bm_RETACCT + bm_UNEMPCOMP + bm_WRKRSCOMP + 
                       bm_FOODSTAMPS + bm_SOCSEC + bm_UNION + bm_CITIZEN, 
                     svydesign = design_p3, data = o3_short1, 
                     method_selection = "logit")
IPW_o3_output <- IPW_o3$output
bchmk_long <- as.data.frame(t(bchmk %>% dplyr::select(ends_with("_1"))))
DF <- np36$mean_SE
nonprob_est <- DF %>% filter(row.names(DF) %in% 
                               c('bm_INSURANCE_1','bm_BLOODPR_1', 'bm_HHPARENT_1',
                                 'bm_FOODALLERGY_1','bm_JOBLASTYR_1','bm_RETACCT_1', 
                                 'bm_UNEMPCOMP_1','bm_WRKRSCOMP_1', 'bm_FOODSTAMPS_1',
                                 'bm_SOCSEC_1','bm_UNION_1','bm_CITIZEN_1'))
# collect estimators and bias
chk <- tibble(qs = rownames(bchmk_long),
              bnchmark = bchmk_long$V1,
              m_CLW_bc = IPW_o3_output$mean,
              m_CLW_bc_bias = round(abs(bnchmark - m_CLW_bc),3),
              v2 = (IPW_o3_output$SE^2)*nrow(o3_short1),
              y1 = p3_output$mean,
              y1_bias = round(abs(bnchmark - y1),3),
              v1 = (p3_output$SE^2)*nrow(data_prob3),
              CLW = nonprob_est$mean,
              CLW_bias = round(abs(bnchmark - CLW),3),
              # bias for EV
              m_bias_est = CLW - y1,
              # EV estimator (where bias is CLW - ps)
              EV = ((m_bias_est^2 + v2)*y1 + v1*CLW)/(m_bias_est^2 + v1 + v2),
              EV_bias = round(abs(bnchmark - EV),3),
              EV_MSE = round(v1*(m_bias_est^2 + v2)/(m_bias_est^2 + v1+v2),3),
              # comb estimator
              comb = (v2/(v1 + v2))*y1 + (v1/(v1 + v2))*m_CLW_bc,
              comb_bias = round(abs(bnchmark - comb),3),
              comb_MSE = round((v1*v2)/(v1+v2),3)
              )

print(paste0("ABS bias ",mean(chk$y1_bias)))  
print(paste0("Original CLW bias ",mean(chk$CLW_bias)))
print(paste0("Model predicted CLW bias ",mean(chk$m_CLW_bc_bias)))
print(paste0("EV bias ",mean(chk$EV_bias)))    
print(paste0("Comb bias ",mean(chk$comb_bias)))    


unknown_final <- chk %>% dplyr::select(qs,bnchmark,EV,EV_MSE,comb,comb_MSE)

# output these tables for latex
table_title1=paste0("Composite Estimates")
print(xtable(unknown_final, caption = table_title1,digits=3),
      caption.placement = 'top',include.colnames =T)


write.csv(o3_long,"/Users/aditisen/Downloads/KML work/o3_long.csv")    
    
    