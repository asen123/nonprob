
## Author : Aditi Sen
## Date : Oct 2024

############################################################################################################
# Topic : IPW, bias-corrected and composite estimates for the known bias case
############################################################################################################

library(survey); library(nonprobsvy); library(dplyr); library(fastDummies)
library(haven); library(labelled); library(ggplot2); library(tidyr); library(sqldf)
library(caret);library(dplyr);library(lme4);library(tidyverse);library(fastDummies)
library(PNWColors);library(gganimate);library(ggeasy);library(ggpubr); library(xtable)
library(patchwork);library(boot);library(rpart);library(rpart.plot);
library(partykit);library(gbm);library(pROC); library(gridExtra);
library(randomForest);library(neuralnet); library(nnet); library(gbm);

dat <- read_sav("/Users/aditisen/Downloads/KML work/2021_benchmarking_study_public_release.sav")
dat2 <- to_factor(dat) # convert all of the labelled variables to factors

bchmk <- read.csv("/Users/aditisen/Downloads/KML work/benchmark1.csv")
bchmk_age <- read.csv("/Users/aditisen/Downloads/KML work/benchmark_age1.csv")
bchmk_race <- read.csv("/Users/aditisen/Downloads/KML work/benchmark_race1.csv")
bchmk_edu <- read.csv("/Users/aditisen/Downloads/KML work/benchmark_edu1.csv")
qs_lookup <- read.csv("/Users/aditisen/Downloads/KML work/qs_lookup.csv")

# Keep relevant variables only (including 12 benchmarks) and 
# remove refused values of auxiliary variables
data_KML <- dat %>%  
  dplyr::select(SAMPID, CASEID, GENDER, AGECAT3BM, RACETHN, REGION, EDUC3, WEIGHT, 
                # 12 benchmark variables
         bm_INSURANCE, bm_BLOODPR, bm_HHPARENT, bm_FOODALLERGY, bm_JOBLASTYR, bm_RETACCT,     
         bm_UNEMPCOMP, bm_WRKRSCOMP, bm_FOODSTAMPS, bm_SOCSEC, bm_UNION, bm_CITIZEN)  %>%
  # remove refused values
  filter(across(c(AGECAT3BM, RACETHN, EDUC3, GENDER, REGION), ~ . != 99)) %>% 
  mutate_at(c('AGECAT3BM','RACETHN','EDUC3','GENDER','REGION'),as.factor) %>% 
  # change 2 to 0
  mutate(across(starts_with("bm"),~recode(.,"2"=0)))

# check distribution of auxiliary variables
dist_func <- function(sample_id){
  data_surv <- filter(data_KML, SAMPID == sample_id)
  age_dist <- data_surv %>% group_by(AGECAT3BM) %>% 
    summarise(ct=n(),pct=n()*100/nrow(data_surv))
  race_dist <- data_surv %>% group_by(RACETHN) %>% 
    summarise(ct=n(),pct=n()*100/nrow(data_surv))
  print(paste("Distribution of age for sample",sample_id, sep = " "))
  print(age_dist)
  print(paste("Distribution of race for sample", sample_id, sep = " "))
  print(race_dist)
}
dist_func(1); dist_func(2); dist_func(3)
dist_func(4); dist_func(5); dist_func(6)

# mention population size
N_B <- 251588972
# select benchmark variables
bm_col_sub <- colnames(dplyr::select(data_KML,starts_with("bm")))
# create dummy variables and drop columns with 0 and 99 
dat_bm <- dummy_cols(data_KML, select_columns = as.character(bm_col_sub))
new_col<- setdiff(colnames(dat_bm), colnames(dplyr::select(dat_bm,ends_with(c("99")))))
dat_bm <- dat_bm %>% dplyr::select(all_of(new_col))
# dataframe with dummy variables (dropping the original benchmark variables and 
# keeping dummy variables)
dat_new <- dat_bm[,-c(9:ncol(data_KML))]

####################################################################################
# function to create IPW weights and get CLW estimates (using R package nonprobsvy)
####################################################################################
IPW_func <- function(prob_surv, non_prob_surv)
{
 #prob_surv=2; non_prob_surv=5
  # probability survey
  data_prob <- filter(dat_new, SAMPID == prob_surv) 
  # calibrate weight to population total
  data_prob$new_weight <- data_prob$WEIGHT*(N_B/nrow(data_prob))
  # define design and survey estimate from the probability survey
  design_p <- survey::svydesign(ids = ~ CASEID, data = data_prob, weights = ~ new_weight) 
  # non probability survey  
  data_nonprob <- filter(dat_new, SAMPID == non_prob_surv)
  #Variable selection for IPW using SCAD penalty
  # set seed to get same result
  set.seed(12325672)
  IPW_logit <- nonprob(selection = ~ AGECAT3BM + GENDER + RACETHN + EDUC3 + REGION, 
    target = ~ 
      bm_INSURANCE_0 + bm_INSURANCE_1 + # insurance
      bm_BLOODPR_0 + bm_BLOODPR_1 + # blood pressure
      bm_HHPARENT_0 + bm_HHPARENT_1 +# parent
      bm_FOODALLERGY_0 + bm_FOODALLERGY_1 +# food allergy
      bm_JOBLASTYR_0 + bm_JOBLASTYR_1 + # job last year
      bm_RETACCT_0 + bm_RETACCT_1 + # retirement account
      bm_UNEMPCOMP_0 + bm_UNEMPCOMP_1 +# unemployment compensation
      bm_WRKRSCOMP_0 + bm_WRKRSCOMP_1 + # workers compensation
      bm_FOODSTAMPS_0 + bm_FOODSTAMPS_1 +# food stamps
      bm_SOCSEC_0 + bm_SOCSEC_1 + # social security
      bm_UNION_0 + bm_UNION_1 + # union
      bm_CITIZEN_0 + bm_CITIZEN_1, 
    svydesign = design_p, data = data_nonprob, 
    method_selection = "logit")
  # return the inverse propensity weights
  data_nonprob$IPW_weight <- IPW_logit$weights
  # produce estimates and SE by category
  func_cat <- function(input_category, cat_level){
    #input_category = "AGECAT3BM"; cat_level = 1
    data_prob_cat <- filter(data_prob, get(input_category) == cat_level)
    data_nonprob_cat <- filter(data_nonprob, get(input_category) == cat_level)
    # unweighted
    data_cat_uw <- data_nonprob_cat %>% dplyr::select(starts_with("bm"))
    # define design
    data_prob_cat$new_new_weight <- data_prob_cat$new_weight*(N_B/nrow(data_prob_cat))
    design_p1c <- survey::svydesign(ids = ~ CASEID, data = data_prob_cat, 
                                    weights = ~ new_new_weight) 
    # mean and SE for 12 binary benchmark variables
    set.seed(12325672)
    if(input_category == "AGECAT3BM"){
      IPW_logit_c <- nonprob(selection = ~ RACETHN + EDUC3 + REGION, 
        target = ~ 
          bm_INSURANCE_0 + bm_INSURANCE_1 + bm_BLOODPR_0 + bm_BLOODPR_1 + 
          bm_HHPARENT_0 + bm_HHPARENT_1 + bm_FOODALLERGY_0 + bm_FOODALLERGY_1 +
          bm_JOBLASTYR_0 + bm_JOBLASTYR_1 + bm_RETACCT_0 + bm_RETACCT_1 +
          bm_UNEMPCOMP_0 + bm_UNEMPCOMP_1 + bm_WRKRSCOMP_0 + bm_WRKRSCOMP_1 + 
          bm_FOODSTAMPS_0 + bm_FOODSTAMPS_1 + bm_SOCSEC_0 + bm_SOCSEC_1 +
          bm_UNION_0 + bm_UNION_1 + bm_CITIZEN_0 + bm_CITIZEN_1, 
        svydesign = design_p1c, data = data_nonprob_cat, method_selection = "logit")
    }else if(input_category == "RACETHN"){
      IPW_logit_c <- nonprob(selection = ~ AGECAT3BM + GENDER + EDUC3 + REGION, 
        target = ~ 
          bm_INSURANCE_0 + bm_INSURANCE_1 + bm_BLOODPR_0 + bm_BLOODPR_1 + 
          bm_HHPARENT_0 + bm_HHPARENT_1 + bm_FOODALLERGY_0 + bm_FOODALLERGY_1 +
          bm_JOBLASTYR_0 + bm_JOBLASTYR_1 + bm_RETACCT_0 + bm_RETACCT_1 +
          bm_UNEMPCOMP_0 + bm_UNEMPCOMP_1 + bm_WRKRSCOMP_0 + bm_WRKRSCOMP_1 + 
          bm_FOODSTAMPS_0 + bm_FOODSTAMPS_1 + bm_SOCSEC_0 + bm_SOCSEC_1 +
          bm_UNION_0 + bm_UNION_1 + bm_CITIZEN_0 + bm_CITIZEN_1, 
        svydesign = design_p1c, data = data_nonprob_cat, method_selection = "logit")
    }else if(input_category == "EDUC3"){
      IPW_logit_c <- nonprob(selection = ~ AGECAT3BM + RACETHN + GENDER + REGION, 
        target = ~ 
          bm_INSURANCE_0 + bm_INSURANCE_1 + bm_BLOODPR_0 + bm_BLOODPR_1 + 
          bm_HHPARENT_0 + bm_HHPARENT_1 + bm_FOODALLERGY_0 + bm_FOODALLERGY_1 +
          bm_JOBLASTYR_0 + bm_JOBLASTYR_1 + bm_RETACCT_0 + bm_RETACCT_1 +
          bm_UNEMPCOMP_0 + bm_UNEMPCOMP_1 + bm_WRKRSCOMP_0 + bm_WRKRSCOMP_1 + 
          bm_FOODSTAMPS_0 + bm_FOODSTAMPS_1 + bm_SOCSEC_0 + bm_SOCSEC_1 +
          bm_UNION_0 + bm_UNION_1 + bm_CITIZEN_0 + bm_CITIZEN_1, 
        svydesign = design_p1c, data = data_nonprob_cat, method_selection = "logit")
    }
    otpt <- list()
    otpt$w <- as.data.frame(IPW_logit_c$output)
    otpt$uw <- as.data.frame(t(colSums(data_cat_uw,na.rm = TRUE)/nrow(data_cat_uw)))
    return(otpt)
  }
  # age
  age1 <- func_cat("AGECAT3BM",1); age2 <- func_cat("AGECAT3BM",2); age3 <- func_cat("AGECAT3BM",3)
  bm_age1 <- age1[[1]]; bm_age1_un <- age1[[2]]
  bm_age2 <- age2[[1]]; bm_age2_un <- age2[[2]]
  bm_age3 <- age3[[1]]; bm_age3_un <- age3[[2]]
  # race
  race1 <- func_cat("RACETHN",1); race2 <- func_cat("RACETHN",2); race3 <- func_cat("RACETHN",3)
  race4 <- func_cat("RACETHN",4); race5 <- func_cat("RACETHN",5)
  bm_race1 <- race1[[1]]; bm_race1_un <- race1[[2]]
  bm_race2 <- race2[[1]]; bm_race2_un <- race2[[2]]
  bm_race3 <- race3[[1]]; bm_race3_un <- race3[[2]]
  bm_race4 <- race4[[1]]; bm_race4_un <- race4[[2]]
  bm_race5 <- race5[[1]]; bm_race5_un <- race5[[2]]
  # educ
  edu1 <- func_cat("EDUC3",1); edu2 <- func_cat("EDUC3",2); edu3 <- func_cat("EDUC3",3)
  bm_edu1 <- edu1[[1]]; bm_edu1_un <- edu1[[2]]
  bm_edu2 <- edu2[[1]]; bm_edu2_un <- edu2[[2]]
  bm_edu3 <- edu3[[1]]; bm_edu3_un <- edu3[[2]]
  # save estimate of mean and standard error
  age_cat <- as.data.frame(rbind(t(bm_age1$mean),t(bm_age2$mean),t(bm_age3$mean)))
    age_se_cat <- as.data.frame(rbind(t(bm_age1$SE),t(bm_age2$SE),t(bm_age3$SE)))
    data_age_uw <- rbind(bm_age1_un, bm_age2_un, bm_age3_un)
  race_cat <- as.data.frame(rbind(t(bm_race1$mean),t(bm_race2$mean),t(bm_race3$mean)))
    race_se_cat <- as.data.frame(rbind(t(bm_race1$SE),t(bm_race2$SE),t(bm_race3$SE)))
    data_race_uw <- rbind(bm_race1_un, bm_race2_un, bm_race3_un)
  edu_cat <- as.data.frame(rbind(t(bm_edu1$mean),t(bm_edu2$mean),t(bm_edu3$mean)))
    edu_se_cat <- as.data.frame(rbind(t(bm_edu1$SE),t(bm_edu2$SE),t(bm_edu3$SE)))
    data_edu_uw <- rbind(bm_edu1_un, bm_edu2_un, bm_edu3_un)
  # set column names of the output tables
  colnames(age_cat) <- colnames(age_se_cat) <- 
    colnames(race_cat) <- colnames(race_se_cat) <- 
      colnames(edu_cat) <- colnames(edu_se_cat) <-
      c('bm_INSURANCE_0','bm_INSURANCE_1','bm_BLOODPR_0','bm_BLOODPR_1','bm_HHPARENT_0',
        'bm_HHPARENT_1','bm_FOODALLERGY_0','bm_FOODALLERGY_1','bm_JOBLASTYR_0','bm_JOBLASTYR_1',
        'bm_RETACCT_0','bm_RETACCT_1','bm_UNEMPCOMP_0','bm_UNEMPCOMP_1',
        'bm_WRKRSCOMP_0','bm_WRKRSCOMP_1','bm_FOODSTAMPS_0','bm_FOODSTAMPS_1','bm_SOCSEC_0',
        'bm_SOCSEC_1','bm_UNION_0','bm_UNION_1','bm_CITIZEN_0','bm_CITIZEN_1')
  # output data
  otpt_IPW <- list()
  otpt_IPW$mean_SE <- IPW_logit$output
  otpt_IPW$otpt_data <- data_nonprob[,c('CASEID','IPW_weight')]
  otpt_IPW$age_cat <- age_cat; otpt_IPW$age_se_cat <- age_se_cat; otpt_IPW$data_age_uw <- data_age_uw
  otpt_IPW$race_cat <- race_cat; otpt_IPW$race_se_cat <- race_se_cat; otpt_IPW$data_race_uw <- data_race_uw
  otpt_IPW$edu_cat <- edu_cat; otpt_IPW$edu_se_cat <- edu_se_cat; otpt_IPW$data_edu_uw <- data_edu_uw
  return(otpt_IPW)
}
# Run function
# NP 1 with P 1 as reference
np14 <- IPW_func(1, 4); np1_weight <- np14[[2]][,2]; 
  hist(np1_weight$IPW_weight); sum(np1_weight$IPW_weight) 
# NP 2 with P 2 as reference
np25 <- IPW_func(2, 5); np2_weight <- np25[[2]][,2]; 
  hist(np2_weight$IPW_weight); sum(np2_weight$IPW_weight) 
# NP 3 with P 3 as reference
np36 <- IPW_func(3, 6); np3_weight <- np36[[2]][,2]; 
  hist(np3_weight$IPW_weight); sum(np3_weight$IPW_weight) 

###################################################################################
# function to calculate bias, takes sample id, weight, benchmark as input 
###################################################################################
IPW_func2 <- function(sample_id, np_weight, b_all, b_age, b_race, b_edu)
{
  #sample_id = 4; np_weight <- np1_weight$IPW_weight
  # filter out the sample needed by the sample ID
  data_surv <- filter(dat_new, SAMPID == sample_id)
  # since weight was adjusted to sample size we multiply weight with calibration 
  # factor to align to the population size
  adj_fac <- N_B/nrow(data_surv); data_surv$KML_weight <- (data_surv$WEIGHT)*(adj_fac)
  if(sample_id < 4) {# for probability sample use given weight
    data_surv$new_weight <- data_surv$KML_weight
  }else {# for nonprobability sample use IPW weight and given weight using function
    data_surv$new_weight <- np_weight
  }
  #####################################################################################
  # define function to get weighted mean for 12 binary benchmark variables 
  # overall and by category
  all_sum <- function(input_weight, input_category){
    #input_weight = "new_weight"; input_category = "AGECAT3BM"
    # for nonprobability samples considering IPW weights get CLW estimate and variance 
    # from previous function IPW_func() outputs 
    if(sample_id == 4 & input_weight == "new_weight"){ bm_p1 <- np14$mean_SE
    } else if(sample_id == 5 & input_weight == "new_weight"){ bm_p1 <- np25$mean_SE
    } else if(sample_id == 6 & input_weight == "new_weight"){ bm_p1 <- np36$mean_SE
    } else{
    # for nonprobability samples considering Pew weights and for probability samples 
    # use survey package
      design_p1 <- survey::svydesign(ids = ~ CASEID, data = data_surv, weights =~ get(input_weight))
      bm_p1 <- survey::svymean(x = ~ 
                                 bm_INSURANCE_0 + bm_INSURANCE_1 + bm_BLOODPR_0 + bm_BLOODPR_1 + 
                                 bm_HHPARENT_0 + bm_HHPARENT_1 + bm_FOODALLERGY_0 + bm_FOODALLERGY_1 +
                                 bm_JOBLASTYR_0 + bm_JOBLASTYR_1 + bm_RETACCT_0 + bm_RETACCT_1 +
                                 bm_UNEMPCOMP_0 + bm_UNEMPCOMP_1 + bm_WRKRSCOMP_0 + bm_WRKRSCOMP_1 + 
                                 bm_FOODSTAMPS_0 + bm_FOODSTAMPS_1 + bm_SOCSEC_0 + bm_SOCSEC_1 +
                                 bm_UNION_0 + bm_UNION_1 + bm_CITIZEN_0 + bm_CITIZEN_1,
                               design = design_p1)
    }
    bm_p11 <- as.data.frame(bm_p1)  
    # save estimate of mean and standard error
    col_sum <- as.data.frame(t(bm_p11$mean)); col_se <- as.data.frame(t(bm_p11$SE))
    colnames(col_sum) <- colnames(col_se) <- 
      c('bm_INSURANCE_0','bm_INSURANCE_1','bm_BLOODPR_0','bm_BLOODPR_1','bm_HHPARENT_0',
        'bm_HHPARENT_1','bm_FOODALLERGY_0','bm_FOODALLERGY_1','bm_JOBLASTYR_0','bm_JOBLASTYR_1',
        'bm_RETACCT_0','bm_RETACCT_1','bm_UNEMPCOMP_0','bm_UNEMPCOMP_1',
        'bm_WRKRSCOMP_0','bm_WRKRSCOMP_1','bm_FOODSTAMPS_0','bm_FOODSTAMPS_1','bm_SOCSEC_0',
        'bm_SOCSEC_1','bm_UNION_0','bm_UNION_1','bm_CITIZEN_0','bm_CITIZEN_1')
    ############################################################################################
    # function for weighted and unweighted sum by category
    func_cat <- function(cat_level){
      #cat_level <- 1 
        # weighted sum
        data_cat <- filter(data_surv, get(input_category) == cat_level)
        # define design
        design_p1c <- survey::svydesign(ids = ~ CASEID, data = data_cat, weights =~get(input_weight)) 
        # mean and SE for 12 binary benchmark variables
        bm_p1c <- survey::svymean(x = ~ 
                                    bm_INSURANCE_0 + bm_INSURANCE_1 + bm_BLOODPR_0 + bm_BLOODPR_1 + 
                                    bm_HHPARENT_0 + bm_HHPARENT_1 + bm_FOODALLERGY_0 + bm_FOODALLERGY_1 +
                                    bm_JOBLASTYR_0 + bm_JOBLASTYR_1 + bm_RETACCT_0 + bm_RETACCT_1 +
                                    bm_UNEMPCOMP_0 + bm_UNEMPCOMP_1 + bm_WRKRSCOMP_0 + bm_WRKRSCOMP_1 + 
                                    bm_FOODSTAMPS_0 + bm_FOODSTAMPS_1 + bm_SOCSEC_0 + bm_SOCSEC_1 +
                                    bm_UNION_0 + bm_UNION_1 + bm_CITIZEN_0 + bm_CITIZEN_1,
                                 design = design_p1c)
      bm_p11c <- as.data.frame(bm_p1c)  
      # save estimate of mean and standard error
      col_sum_cat <- as.data.frame(t(bm_p11c$mean)); col_se_cat <- as.data.frame(t(bm_p11c$SE))
      colnames(col_sum_cat) <- colnames(col_se_cat) <- 
        c('bm_INSURANCE_0','bm_INSURANCE_1','bm_BLOODPR_0','bm_BLOODPR_1','bm_HHPARENT_0',
          'bm_HHPARENT_1','bm_FOODALLERGY_0','bm_FOODALLERGY_1','bm_JOBLASTYR_0','bm_JOBLASTYR_1',
          'bm_RETACCT_0','bm_RETACCT_1','bm_UNEMPCOMP_0','bm_UNEMPCOMP_1',
          'bm_WRKRSCOMP_0','bm_WRKRSCOMP_1','bm_FOODSTAMPS_0','bm_FOODSTAMPS_1','bm_SOCSEC_0',
          'bm_SOCSEC_1','bm_UNION_0','bm_UNION_1','bm_CITIZEN_0','bm_CITIZEN_1')
      # unweighted sum
      data_cat_uw <- data_cat %>% dplyr::select(starts_with("bm"))
      col_sum_cat_uw <- as.data.frame(t(colSums(data_cat_uw,na.rm = TRUE)/nrow(data_cat_uw)))
      # function output
      output_cat <- list()
      output_cat$col_sum_cat <- col_sum_cat; output_cat$col_se_cat <- col_se_cat
      output_cat$col_sum_cat_uw <- col_sum_cat_uw; #output_cat$col_se_cat_uw <- col_se_cat_uw
      return(output_cat)
    }
    # for nonprobability samples considering IPW weights get CLW estimate and variance 
    # from previous function IPW_func() outputs
    if(sample_id == 4 & input_weight == "new_weight" & input_category == "AGECAT3BM"){ 
      data_cat <- np14$age_cat; data_cat_se <- np14$age_se_cat; data_cat_uw <- np14$data_age_uw
    } else if(sample_id == 5 & input_weight == "new_weight" & input_category == "AGECAT3BM"){ 
      data_cat <- np25$age_cat; data_cat_se <- np25$age_se_cat; data_cat_uw <- np25$data_age_uw
    } else if(sample_id == 6 & input_weight == "new_weight" & input_category == "AGECAT3BM"){ 
      data_cat <- np36$age_cat; data_cat_se <- np36$age_se_cat; data_cat_uw <- np36$data_age_uw
    } else if(sample_id == 4 & input_weight == "new_weight" & input_category == "RACETHN"){ 
      data_cat <- np14$race_cat; data_cat_se <- np14$race_se_cat; data_cat_uw <- np14$data_race_uw
    } else if(sample_id == 5 & input_weight == "new_weight" & input_category == "RACETHN"){ 
      data_cat <- np25$race_cat; data_cat_se <- np25$race_se_cat; data_cat_uw <- np25$data_race_uw
    } else if(sample_id == 6 & input_weight == "new_weight" & input_category == "RACETHN"){ 
      data_cat <- np36$race_cat; data_cat_se <- np36$race_se_cat; data_cat_uw <- np36$data_race_uw
    } else if(sample_id == 4 & input_weight == "new_weight" & input_category == "EDUC3"){ 
      data_cat <- np14$edu_cat; data_cat_se <- np14$edu_se_cat; data_cat_uw <- np14$data_edu_uw
    } else if(sample_id == 5 & input_weight == "new_weight" & input_category == "EDUC3"){ 
      data_cat <- np25$edu_cat; data_cat_se <- np25$edu_se_cat; data_cat_uw <- np25$data_edu_uw
    } else if(sample_id == 6 & input_weight == "new_weight" & input_category == "EDUC3"){ 
      data_cat <- np36$edu_cat; data_cat_se <- np36$edu_se_cat; data_cat_uw <- np36$data_edu_uw
    }
    else{
    # for probability samples considering Pew weights use function func_cat()
      cat1 <- func_cat(1); cat2 <- func_cat(2); cat3 <- func_cat(3) 
      # weighted sum combining all categories
      data_cat <- as.data.frame(rbind(cat1$col_sum_cat, cat2$col_sum_cat, cat3$col_sum_cat))
      data_cat_se <- as.data.frame(rbind(cat1$col_se_cat,cat2$col_se_cat,cat3$col_se_cat))
      # unweighted sum combining all categories
      data_cat_uw <- as.data.frame(rbind(cat1$col_sum_cat_uw, cat2$col_sum_cat_uw,cat3$col_sum_cat_uw))
    }
    # function output
    output_all_sum <- list()
    output_all_sum$col_sum <- col_sum; output_all_sum$col_se <- col_se
    output_all_sum$data_cat <- data_cat; output_all_sum$data_cat_se <- data_cat_se
    output_all_sum$data_cat_uw <- data_cat_uw
    return(output_all_sum)
  }
  # KML weighted estimates
  f_age_KML <- all_sum("KML_weight","AGECAT3BM") # age
  f_race_KML <- all_sum("KML_weight","RACETHN") # race
  f_edu_KML <- all_sum("KML_weight","EDUC3") # edu
  # IPW weights estimates
  f_age_IPW <- all_sum("new_weight","AGECAT3BM") # age
  f_race_IPW <- all_sum("new_weight","RACETHN") # race
  f_edu_IPW <- all_sum("new_weight","EDUC3") # edu
  # condition for probability and nonprobability samples
  if(sample_id < 4)
  {# for probability sample use given weight as KML, there is no IPW
    col_sum_KML <- col_sum_IPW <- f_age_KML$col_sum # overall mean
    col_se_KML <- col_se_IPW <- f_age_KML$col_se # overall se
    # age
    data_age_uw <- f_age_KML$data_cat_uw 
    data_age_KML <- data_age_IPW <- f_age_KML$data_cat # mean by age (KML/IPW weight is same)
    data_age_KML_se <- data_age_IPW_se <- f_age_KML$data_cat_se # se by age
    # race
    data_race_uw <- f_race_KML$data_cat_uw 
    data_race_KML <- data_race_IPW <- f_race_KML$data_cat  # mean by race 
    data_race_KML_se <- data_race_IPW_se <- f_race_KML$data_cat_se  # se by race 
    # education
    data_edu_uw <- f_edu_KML$data_cat_uw 
    data_edu_IPW <- data_edu_KML <- f_edu_KML$data_cat  # mean by edu
    data_edu_IPW_se <- data_edu_KML_se <- f_edu_KML$data_cat_se # mean by edu 
  }else
  {
    # for nonprobability sample use IPW weight and given weight using function
    col_sum_IPW <- f_age_IPW$col_sum; col_se_IPW <- f_age_IPW$col_se 
    col_sum_KML <- f_age_KML$col_sum; col_se_KML <- f_age_KML$col_se
    # age
    data_age_uw <- f_age_KML$data_cat_uw 
    data_age_IPW <- f_age_IPW$data_cat; data_age_IPW_se <- f_age_IPW$data_cat_se
    data_age_KML <- f_age_KML$data_cat; data_age_KML_se <- f_age_KML$data_cat_se
    # race
    data_race_uw <- f_race_KML$data_cat_uw 
    data_race_IPW <- f_race_IPW$data_cat; data_race_IPW_se <- f_race_IPW$data_cat_se
    data_race_KML <- f_race_KML$data_cat; data_race_KML_se <- f_race_KML$data_cat_se
    # education
    data_edu_uw <- f_edu_KML$data_cat_uw 
    data_edu_IPW <- f_edu_IPW$data_cat; data_edu_IPW_se <- f_edu_IPW$data_cat_se
    data_edu_KML <- f_edu_KML$data_cat; data_edu_KML_se <- f_edu_KML$data_cat_se
  }
  # unweighted
  col_sum_uw <- as.data.frame(t(colSums(data.frame(data_surv[,-c(1:8,ncol(data_surv)-1,
                                                                 ncol(data_surv))]),
                                        na.rm = TRUE)/nrow(data_surv)))
  # drop refused columns
  new_col <- setdiff(colnames(col_sum_KML), colnames(dplyr::select(col_sum_KML,ends_with(c("NA","99")))))
  # select other than refused columns 
  col_sum_1_uw <- col_sum_uw %>% dplyr::select(all_of(new_col));
  col_sum_1_KML <- col_sum_KML %>% dplyr::select(all_of(new_col)); 
    col_se_1_KML <- col_se_KML %>% dplyr::select(all_of(new_col))
  col_sum_1_IPW <- col_sum_IPW %>% dplyr::select(all_of(new_col)); 
    col_se_1_IPW <- col_se_IPW %>% dplyr::select(all_of(new_col))
  col_sum_age_uw <- data_age_uw %>% dplyr::select(all_of(new_col))
  col_sum_age_KML <- data_age_KML %>% dplyr::select(all_of(new_col)); 
    col_se_age_KML <- data_age_KML_se %>% dplyr::select(all_of(new_col))
  col_sum_age_IPW <- data_age_IPW %>% dplyr::select(all_of(new_col)); 
    col_se_age_IPW <- data_age_IPW_se %>% dplyr::select(all_of(new_col))
  col_sum_race_uw <- data_race_uw %>% dplyr::select(all_of(new_col))
  col_sum_race_KML <- data_race_KML %>% dplyr::select(all_of(new_col)); 
    col_se_race_KML <- data_race_KML_se %>% dplyr::select(all_of(new_col))
  col_sum_race_IPW <- data_race_IPW %>% dplyr::select(all_of(new_col)); 
    col_se_race_IPW <- data_race_IPW_se %>% dplyr::select(all_of(new_col))
  col_sum_edu_uw <- data_edu_uw %>% dplyr::select(all_of(new_col))
  col_sum_edu_KML <- data_edu_KML %>% dplyr::select(all_of(new_col)); 
    col_se_edu_KML <- data_edu_KML_se %>% dplyr::select(all_of(new_col))
  col_sum_edu_IPW <- data_edu_IPW %>% dplyr::select(all_of(new_col)); 
    col_se_edu_IPW <- data_edu_IPW_se %>% dplyr::select(all_of(new_col))
  # mean from all benchmark variables
  # create a blank dataframe
  yc_bar_uw <- yc_bar_KML <- yc_bar_KML_se <- yc_bar_IPW <- yc_bar_IPW_se <-
    yc_bar_age_uw <- yc_bar_age_KML <-yc_bar_age_KML_se <- 
    yc_bar_age_IPW <- yc_bar_age_IPW_se <-
    yc_bar_race_uw <- yc_bar_race_KML <- yc_bar_race_KML_se <- 
    yc_bar_race_IPW <- yc_bar_race_IPW_se <- 
    yc_bar_edu_uw <- yc_bar_edu_KML <- yc_bar_edu_KML_se <- 
    yc_bar_edu_IPW <- yc_bar_edu_IPW_se <- 
    data.frame(matrix(nrow = 0, ncol = length(bm_col_sub))) 
  # set the column names
  colnames(yc_bar_uw) <- 
    colnames(yc_bar_KML) <- colnames(yc_bar_KML_se) <- 
    colnames(yc_bar_IPW) <- colnames(yc_bar_IPW_se) <- 
    colnames(yc_bar_age_uw) <- 
    colnames(yc_bar_age_KML) <- colnames(yc_bar_age_KML_se) <-  
    colnames(yc_bar_age_IPW) <- colnames(yc_bar_age_IPW_se) <-
    colnames(yc_bar_race_uw) <- 
    colnames(yc_bar_race_KML) <- colnames(yc_bar_race_KML_se) <- 
    colnames(yc_bar_race_IPW) <- colnames(yc_bar_race_IPW_se) <- 
    colnames(yc_bar_edu_uw) <- 
    colnames(yc_bar_edu_KML) <- colnames(yc_bar_edu_KML_se) <- 
    colnames(yc_bar_edu_IPW) <- colnames(yc_bar_edu_IPW_se) <- 
    bm_col_sub
  # create blank dataframe
  x_i <- y_i_uw <- y_i_KML <- y_i_KML_se <- y_i_IPW <- y_i_IPW_se <-
    bias_i_uw <- bias_i_KML  <- bias_i_IPW <- 
    data.frame(matrix(nrow = 1, ncol = 0)) 
  x_i_age <- y_i_age_uw <- y_i_age_KML <- y_i_age_KML_se <- y_i_age_IPW <- y_i_age_IPW_se <-
    bias_i_age_uw <- bias_i_age_KML <- bias_i_age_IPW <-
    x_i_race <- y_i_race_uw <- y_i_race_KML <- y_i_race_KML_se <- y_i_race_IPW <- y_i_race_IPW_se <-
    bias_i_race_uw <- bias_i_race_KML <- bias_i_race_IPW <- 
    x_i_edu <- y_i_edu_uw <- y_i_edu_KML <- y_i_edu_KML_se <- y_i_edu_IPW <- y_i_edu_IPW_se <-
    bias_i_edu_uw <- bias_i_edu_KML <- bias_i_edu_IPW <-
    data.frame(matrix(nrow = 3, ncol = 0)) 
  
  for(i in bm_col_sub){
    #i = "bm_INSURANCE"
    # identify number of levels of each benchmark variable
    k <- length(dplyr::select(col_sum_1_KML,starts_with(i)))
    # calculate difference with benchmark variables
    diff_uw <- dplyr::select(col_sum_1_uw,starts_with(i)) - dplyr::select(b_all,starts_with(i))
    diff_KML <- dplyr::select(col_sum_1_KML,starts_with(i)) - dplyr::select(b_all,starts_with(i))
    diff_IPW <- dplyr::select(col_sum_1_IPW,starts_with(i)) - dplyr::select(b_all,starts_with(i))
    # calculate absolute difference and summarize at variable level
    yc_bar_uw[1,which(bm_col_sub==i)] <- sum(abs(diff_uw))/k
    yc_bar_KML[1,which(bm_col_sub==i)] <- sum(abs(diff_KML))/k
    yc_bar_IPW[1,which(bm_col_sub==i)] <- sum(abs(diff_IPW))/k
    # summarize standard error at variable level
    yc_bar_KML_se[1,which(bm_col_sub==i)] <- 
      rowMeans(dplyr::select(col_se_1_KML,starts_with(i)))
    yc_bar_IPW_se[1,which(bm_col_sub==i)] <- 
      rowMeans(as.matrix(dplyr::select(col_se_1_IPW,starts_with(i))))
    # calculate difference by age and summarize at variable level
    diff_age_uw <- dplyr::select(col_sum_age_uw,starts_with(i)) - 
      dplyr::select(b_age[,-1],starts_with(i))
    diff_age_KML <- dplyr::select(col_sum_age_KML,starts_with(i)) - 
      dplyr::select(b_age[,-1],starts_with(i))
    diff_age_IPW <- dplyr::select(col_sum_age_IPW,starts_with(i)) - 
      dplyr::select(b_age[,-1],starts_with(i))
    yc_bar_age_uw[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_age_uw)))/k
    yc_bar_age_KML[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_age_KML)))/k
    yc_bar_age_IPW[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_age_IPW)))/k
    # summarize standard error at variable level
    yc_bar_age_KML_se[c(1:3),which(bm_col_sub==i)] <- 
      rowMeans(dplyr::select(col_se_age_KML,starts_with(i)))
    yc_bar_age_IPW_se[c(1:3),which(bm_col_sub==i)] <- 
      rowMeans(dplyr::select(col_se_age_IPW,starts_with(i)))
    # calculate difference by race and summarize at variable level  
    diff_race_uw <- dplyr::select(col_sum_race_uw,starts_with(i)) - 
      dplyr::select(bchmk_race[,-1],starts_with(i))
    diff_race_KML <- dplyr::select(col_sum_race_KML,starts_with(i)) - 
      dplyr::select(bchmk_race[,-1],starts_with(i))
    diff_race_IPW <- dplyr::select(col_sum_race_IPW,starts_with(i)) - 
      dplyr::select(bchmk_race[,-1],starts_with(i))
    yc_bar_race_uw[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_race_uw)))/k
    yc_bar_race_KML[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_race_KML)))/k
    yc_bar_race_IPW[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_race_IPW)))/k
    # summarize standard error at variable level
    yc_bar_race_KML_se[c(1:3),which(bm_col_sub==i)] <- 
      rowMeans(dplyr::select(col_se_race_KML,starts_with(i)))
    yc_bar_race_IPW_se[c(1:3),which(bm_col_sub==i)] <- 
      rowMeans(dplyr::select(col_se_race_IPW,starts_with(i)))
    # calculate difference by education and summarize at variable level 
    diff_edu_uw <- dplyr::select(col_sum_edu_uw,starts_with(i)) - 
      dplyr::select(bchmk_edu[,-1],starts_with(i))
    diff_edu_KML <- dplyr::select(col_sum_edu_KML,starts_with(i)) - 
      dplyr::select(bchmk_edu[,-1],starts_with(i))
    diff_edu_IPW <- dplyr::select(col_sum_edu_IPW,starts_with(i)) - 
      dplyr::select(bchmk_edu[,-1],starts_with(i))
    yc_bar_edu_uw[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_edu_uw)))/k
    yc_bar_edu_KML[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_edu_KML)))/k
    yc_bar_edu_IPW[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_edu_IPW)))/k
    # summarize standard error at variable level 
    yc_bar_edu_KML_se[c(1:3),which(bm_col_sub==i)] <- 
      rowMeans(dplyr::select(col_se_edu_KML,starts_with(i)))
    yc_bar_edu_IPW_se[c(1:3),which(bm_col_sub==i)] <- 
      rowMeans(dplyr::select(col_se_edu_IPW,starts_with(i)))
    # store values of estimates and bias 
    if(k == 2){
      y_i_uw <- cbind(y_i_uw, dplyr::select(col_sum_1_uw, starts_with(i)))
      y_i_KML <- cbind(y_i_KML, dplyr::select(col_sum_1_KML, starts_with(i)))
        y_i_KML_se <- cbind(y_i_KML_se, dplyr::select(col_se_1_KML, starts_with(i)))
      y_i_IPW <- cbind(y_i_IPW, dplyr::select(col_sum_1_IPW, starts_with(i)))
        y_i_IPW_se <- cbind(y_i_IPW_se, dplyr::select(col_se_1_IPW, starts_with(i)))
      x_i <- cbind(x_i, dplyr::select(bchmk, starts_with(i)))
      bias_i_uw <- cbind(bias_i_uw, diff_uw)
      bias_i_KML <- cbind(bias_i_KML, diff_KML)
      bias_i_IPW <- cbind(bias_i_IPW, diff_IPW)
      # age
      y_i_age_uw <- cbind(y_i_age_uw, dplyr::select(col_sum_age_uw, starts_with(i)))
      y_i_age_KML <- cbind(y_i_age_KML, dplyr::select(col_sum_age_KML, starts_with(i)))
        y_i_age_KML_se <- cbind(y_i_age_KML_se, dplyr::select(col_se_age_KML, starts_with(i)))
      y_i_age_IPW <- cbind(y_i_age_IPW, dplyr::select(col_sum_age_IPW, starts_with(i)))
        y_i_age_IPW_se <- cbind(y_i_age_IPW_se, dplyr::select(col_se_age_IPW, starts_with(i)))
      x_i_age <- cbind(x_i_age, dplyr::select(b_age[,-1], starts_with(i)))
      bias_i_age_uw <- cbind(bias_i_age_uw, diff_age_uw)
      bias_i_age_KML <- cbind(bias_i_age_KML, diff_age_KML)
      bias_i_age_IPW <- cbind(bias_i_age_IPW, diff_age_IPW)
      # race
      y_i_race_uw <- cbind(y_i_race_uw, dplyr::select(col_sum_race_uw, starts_with(i)))
      y_i_race_KML <- cbind(y_i_race_KML, dplyr::select(col_sum_race_KML, starts_with(i)))
        y_i_race_KML_se <- cbind(y_i_race_KML_se, dplyr::select(col_se_race_KML, starts_with(i)))
      y_i_race_IPW <- cbind(y_i_race_IPW, dplyr::select(col_sum_race_IPW, starts_with(i)))
        y_i_race_IPW_se <- cbind(y_i_race_IPW_se, dplyr::select(col_se_race_IPW, starts_with(i)))
      x_i_race <- cbind(x_i_race, dplyr::select(bchmk_race[,-1],starts_with(i)))
      bias_i_race_uw <- cbind(bias_i_race_uw, diff_race_uw)
      bias_i_race_KML <- cbind(bias_i_race_KML, diff_race_KML)
      bias_i_race_IPW <- cbind(bias_i_race_IPW, diff_race_IPW)
      # education
      y_i_edu_uw <- cbind(y_i_edu_uw, dplyr::select(col_sum_edu_uw, starts_with(i)))
      y_i_edu_KML <- cbind(y_i_edu_KML, dplyr::select(col_sum_edu_KML, starts_with(i)))
        y_i_edu_KML_se <- cbind(y_i_edu_KML_se, dplyr::select(col_se_edu_KML, starts_with(i)))
      y_i_edu_IPW <- cbind(y_i_edu_IPW, dplyr::select(col_sum_edu_IPW, starts_with(i)))
        y_i_edu_IPW_se <- cbind(y_i_edu_IPW_se, dplyr::select(col_se_edu_IPW, starts_with(i)))
      x_i_edu <- cbind(x_i_edu, dplyr::select(bchmk_edu[,-1],starts_with(i)))
      bias_i_edu_uw <- cbind(bias_i_edu_uw, diff_edu_uw)
      bias_i_edu_KML <- cbind(bias_i_edu_KML, diff_edu_KML)
      bias_i_edu_IPW <- cbind(bias_i_edu_IPW, diff_edu_IPW)
    }
  }
  # output dataset
  out_put_val <- list()
  # calculate mean of difference for 12 benchmark variables  
  out_put_val$mean_uw <- rowSums(yc_bar_uw)/length(bm_col_sub)
  out_put_val$mean_KML <- rowSums(yc_bar_KML)/length(bm_col_sub)
  out_put_val$mean_IPW <- rowSums(yc_bar_IPW)/length(bm_col_sub)
  # group by age
  out_put_val$mean_age_uw <- rowSums(yc_bar_age_uw)/length(bm_col_sub)
  out_put_val$mean_age_KML <- rowSums(yc_bar_age_KML)/length(bm_col_sub)
  out_put_val$mean_age_IPW <- rowSums(yc_bar_age_IPW)/length(bm_col_sub)
  # group by race
  out_put_val$mean_race_uw <- rowSums(yc_bar_race_uw)/length(bm_col_sub)
  out_put_val$mean_race_KML <- rowSums(yc_bar_race_KML)/length(bm_col_sub)
  out_put_val$mean_race_IPW <- rowSums(yc_bar_race_IPW)/length(bm_col_sub)
  # group by education
  out_put_val$mean_edu_uw <- rowSums(yc_bar_edu_uw)/length(bm_col_sub)
  out_put_val$mean_edu_KML <- rowSums(yc_bar_edu_KML)/length(bm_col_sub)
  out_put_val$mean_edu_IPW <- rowSums(yc_bar_edu_IPW)/length(bm_col_sub)
  # original values, benchmark and bias  
  out_put_val$y_i_uw <- y_i_uw %>% dplyr::select(ends_with("1"))
  out_put_val$y_i_KML <- y_i_KML %>% dplyr::select(ends_with("1"))
    out_put_val$y_i_KML_se <- y_i_KML_se %>% dplyr::select(ends_with("1"))
  out_put_val$y_i_IPW <- y_i_IPW %>% dplyr::select(ends_with("1"))
    out_put_val$y_i_IPW_se <- y_i_IPW_se %>% dplyr::select(ends_with("1"))
  out_put_val$x_i <- x_i %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_uw <- bias_i_uw %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_KML <- bias_i_KML %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_IPW <- bias_i_IPW %>% dplyr::select(ends_with("1"))
  # age
  out_put_val$y_i_age_uw <- y_i_age_uw %>% dplyr::select(ends_with("1"))
  out_put_val$y_i_age_KML <- y_i_age_KML %>% dplyr::select(ends_with("1"))
    out_put_val$y_i_age_KML_se <- y_i_age_KML_se %>% dplyr::select(ends_with("1"))
  out_put_val$y_i_age_IPW <- y_i_age_IPW %>% dplyr::select(ends_with("1"))
  out_put_val$y_i_age_IPW_se <- y_i_age_IPW_se %>% dplyr::select(ends_with("1"))
  out_put_val$x_i_age <- x_i_age %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_age_uw <- bias_i_age_uw %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_age_KML <- bias_i_age_KML %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_age_IPW <- bias_i_age_IPW %>% dplyr::select(ends_with("1"))
  # race
  out_put_val$y_i_race_uw <- y_i_race_uw %>% dplyr::select(ends_with("1"))
  out_put_val$y_i_race_KML <- y_i_race_KML %>% dplyr::select(ends_with("1"))
    out_put_val$y_i_race_KML_se <- y_i_race_KML_se %>% dplyr::select(ends_with("1"))
  out_put_val$y_i_race_IPW <- y_i_race_IPW %>% dplyr::select(ends_with("1"))
    out_put_val$y_i_race_IPW_se <- y_i_race_IPW_se %>% dplyr::select(ends_with("1"))
  out_put_val$x_i_race <- x_i_race %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_race_uw <- bias_i_race_uw %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_race_KML <- bias_i_race_KML %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_race_IPW <- bias_i_race_IPW %>% dplyr::select(ends_with("1"))
  # education
  out_put_val$y_i_edu_uw <- y_i_edu_uw %>% dplyr::select(ends_with("1"))
  out_put_val$y_i_edu_KML <- y_i_edu_KML %>% dplyr::select(ends_with("1"))
    out_put_val$y_i_edu_KML_se <- y_i_edu_KML_se %>% dplyr::select(ends_with("1"))
  out_put_val$y_i_edu_IPW <- y_i_edu_IPW %>% dplyr::select(ends_with("1"))
    out_put_val$y_i_edu_IPW_se <- y_i_edu_IPW_se %>% dplyr::select(ends_with("1"))
  out_put_val$x_i_edu <- x_i_edu %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_edu_uw <- bias_i_edu_uw %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_edu_KML <- bias_i_edu_KML %>% dplyr::select(ends_with("1"))
  out_put_val$bias_i_edu_IPW <- bias_i_edu_IPW %>% dplyr::select(ends_with("1"))
  return(out_put_val)
} 

diff_opt_1 <- IPW_func2(4, np1_weight, b_all = bchmk, b_age = bchmk_age, 
                        b_race = bchmk_race, b_edu = bchmk_edu) # NP 1
diff_opt_2 <- IPW_func2(5, np2_weight, b_all = bchmk, b_age = bchmk_age, 
                        b_race = bchmk_race, b_edu = bchmk_edu) # NP 2
diff_opt_3 <- IPW_func2(6, np3_weight, b_all = bchmk, b_age = bchmk_age, 
                        b_race = bchmk_race, b_edu = bchmk_edu) # NP 3
diff_prob_1 <- IPW_func2(1, np1_weight, b_all = bchmk, b_age = bchmk_age, 
                         b_race = bchmk_race, b_edu = bchmk_edu) # P 1
diff_prob_2 <- IPW_func2(2, np2_weight, b_all = bchmk, b_age = bchmk_age, 
                         b_race = bchmk_race, b_edu = bchmk_edu) # P 2
diff_prob_3 <- IPW_func2(3, np3_weight, b_all = bchmk, b_age = bchmk_age, 
                         b_race = bchmk_race, b_edu = bchmk_edu) # P 3

####################################################################################
# Create table to output MAE values at survey level 
####################################################################################
sampid <- c("P1","P2","P3","O1","O2","O3")
diff_genpop <- matrix(NA, ncol = length(sampid), nrow = 3)
colnames(diff_genpop) <- sampid
row.names(diff_genpop) <- c('Mean (unweighted)', 'Mean (KML weight) ', 'Mean (IPW)')
for(i in 1:3){
  # mean
  diff_mean <- round(c(diff_prob_1[[i]], diff_prob_2[[i]], diff_prob_3[[i]], 
                 diff_opt_1[[i]], diff_opt_2[[i]], diff_opt_3[[i]])*100,1)
  diff_genpop[i,] <- diff_mean
}
rm(i, diff_mean)
####################################################################################
# Create table to output MAE values group by age, race and education 
####################################################################################
age_cat_vec <- c("18-29","30-64","65+")
race_cat_vec <- c("White non-Hispanic", "Black non-Hispanic", "Hispanic")
edu_cat_vec <- c("HS or less", "Some college", "College grad")
diff_age <- matrix(NA, ncol = length(sampid)+1, nrow = length(age_cat_vec)*3)
diff_race <- matrix(NA, ncol = length(sampid)+1, nrow = length(race_cat_vec)*3)
diff_edu <- matrix(NA, ncol = length(sampid)+1, nrow = length(edu_cat_vec)*3)
colnames(diff_age) <- c("age", sampid)
colnames(diff_race) <- c("race", sampid)
colnames(diff_edu) <- c("edu", sampid)
row.names(diff_age) <- 
  c('Mean (unweighted) 18-29 yrs', 'Mean (unweighted) 30-64 yrs', 'Mean (unweighted) 65+ yrs',
  'Mean (KML weight)  18-29 yrs', 'Mean (KML weight) 30-64 yrs', 'Mean (KML weight) 65+ yrs', 
  'Mean (IPW)  18-29 yrs', 'Mean (IPW) 30-64 yrs', 'Mean (IPW) 65+ yrs')
row.names(diff_race) <- 
  c('Mean (unweighted) White non-Hispanic', 'Mean (unweighted) Black non-Hispanic', 'Mean (unweighted) Hispanic',
    'Mean (KML weight)  White non-Hispanic', 'Mean (KML weight) Black non-Hispanic', 'Mean (KML weight) Hispanic', 
    'Mean (IPW)  White non-Hispanic', 'Mean (IPW) Black non-Hispanic', 'Mean (IPW) Hispanic')
row.names(diff_edu) <- 
  c('Mean (unweighted) HS or less', 'Mean (unweighted) Some college', 'Mean (unweighted) College grad',
  'Mean (KML weight)  HS or less', 'Mean (KML weight) Some college', 'Mean (KML weight) College grad', 
  'Mean (IPW)  HS or less', 'Mean (IPW) Some college', 'Mean (IPW) College grad')
counter <- 0
for (i in 1:3) {
  # age
  diff_age_mean <- data.frame(cbind(age_cat_vec,
                                    round(diff_prob_1[[3+i]]*100,1), round(diff_prob_2[[3+i]]*100,1), 
                                    round(diff_prob_3[[3+i]]*100,1), round(diff_opt_1[[3+i]]*100,1), 
                                    round(diff_opt_2[[3+i]]*100,1), round(diff_opt_3[[3+i]]*100,1)
                                    ))
  diff_age_mean <- diff_age_mean %>% mutate(across(where(is.numeric), round, 3))
  diff_age[ c((counter+1):(counter+length(age_cat_vec))) ,] <- as.matrix(diff_age_mean)
  # race
  diff_race_mean <- data.frame(cbind(race_cat_vec,
                                     round(diff_prob_1[[6+i]]*100,1), round(diff_prob_2[[6+i]]*100,1), 
                                     round(diff_prob_3[[6+i]]*100,1), round(diff_opt_1[[6+i]]*100,1), 
                                     round(diff_opt_2[[6+i]]*100,1), round(diff_opt_3[[6+i]]*100,1)
                                                                                     ))
  diff_race[ c((counter+1):(counter+length(race_cat_vec))),] <- as.matrix(diff_race_mean)
  # education
  diff_edu_mean <- data.frame(cbind(edu_cat_vec,
                                    round(diff_prob_1[[9+i]]*100,1), round(diff_prob_2[[9+i]]*100,1), 
                                    round(diff_prob_3[[9+i]]*100,1), round(diff_opt_1[[9+i]]*100,1), 
                                    round(diff_opt_2[[9+i]]*100,1), round(diff_opt_3[[9+i]]*100,1)
                                    ))
  diff_edu[ c((counter+1):(counter+length(edu_cat_vec))) ,] <- as.matrix(diff_edu_mean)
  counter <- counter + length(race_cat_vec)
  rm(diff_age_mean, diff_race_mean, diff_edu_mean)
}
rm(i, counter)

# output these tables for latex
table_title1=paste0("MAE for general population")
print(xtable(diff_genpop, caption = table_title1,digits=2),
      caption.placement = 'top',include.colnames =T)

table_title2=paste0("MAE for age groups")
print(xtable(diff_age, caption = table_title2,digits=2),
      caption.placement = 'top',include.colnames =T)

table_title3=paste0("MAE for race groups")
print(xtable(diff_race, caption = table_title3,digits=2),
      caption.placement = 'top',include.colnames =T)

table_title4=paste0("MAE for education groups")
print(xtable(diff_edu, caption = table_title4,digits=2),
      caption.placement = 'top',include.colnames =T)

################################################################################################################
############ function for calculating composite estimators ############ 
################################################################################################################

model_func <- function(data_surv, samplid){
  ncol_val <- ncol(data_surv$bias_i_age_IPW)
  m_data <- tibble(
    samplid = samplid,
    age = c(rep("18-29",ncol_val),rep("30-64",ncol_val),rep("65+",ncol_val)),
    qs = c(rep(colnames(data_surv$bias_i_age_IPW),3)),
    # age
    x_i_age = c(t(data_surv$x_i_age[1,]),t(data_surv$x_i_age[2,]),t(data_surv$x_i_age[3,])),
    y_i_age_IPW = c(t(data_surv$y_i_age_IPW[1,]),t(data_surv$y_i_age_IPW[2,]),
                    t(data_surv$y_i_age_IPW[3,])), 
    y_i_age_IPW_se = c(t(data_surv$y_i_age_IPW_se[1,]),t(data_surv$y_i_age_IPW_se[2,]),
                       t(data_surv$y_i_age_IPW_se[3,])),
    bias_i_age_IPW = c(t(data_surv$bias_i_age_IPW[1,]),t(data_surv$bias_i_age_IPW[2,]),
                       t(data_surv$bias_i_age_IPW[3,])),
    # race
    race = c(rep("White non-Hispanic",ncol_val),rep("Black non-Hispanic",ncol_val),rep("Hispanic",ncol_val)),
    x_i_race = c(t(data_surv$x_i_race[1,]),t(data_surv$x_i_race[2,]),t(data_surv$x_i_race[3,])),
    y_i_race_IPW = c(t(data_surv$y_i_race_IPW[1,]),t(data_surv$y_i_race_IPW[2,]),
                     t(data_surv$y_i_race_IPW[3,])),
    y_i_race_IPW_se = c(t(data_surv$y_i_race_IPW_se[1,]),t(data_surv$y_i_race_IPW_se[2,]),
                        t(data_surv$y_i_race_IPW_se[3,])),
    bias_i_race_IPW = c(t(data_surv$bias_i_race_IPW[1,]),t(data_surv$bias_i_race_IPW[2,]),
                        t(data_surv$bias_i_race_IPW[3,])),
    # education
    edu = c(rep("HS or less",ncol_val),rep("Some college",ncol_val),rep("College grad",ncol_val)),
    x_i_edu = c(t(data_surv$x_i_edu[1,]),t(data_surv$x_i_edu[2,]),t(data_surv$x_i_edu[3,])),
    y_i_edu_IPW = c(t(data_surv$y_i_edu_IPW[1,]),t(data_surv$y_i_edu_IPW[2,]),
                    t(data_surv$y_i_edu_IPW[3,])),
    y_i_edu_IPW_se = c(t(data_surv$y_i_edu_IPW_se[1,]),t(data_surv$y_i_edu_IPW_se[2,]),
                       t(data_surv$y_i_edu_IPW_se[3,])),
    bias_i_edu_IPW = c(t(data_surv$bias_i_edu_IPW[1,]),t(data_surv$bias_i_edu_IPW[2,]),
                       t(data_surv$bias_i_edu_IPW[3,]))
                  )
  
  m_all <- tibble(mu_IPW = t(data_surv$y_i_IPW[1,]), 
                  mu_IPW_se = t(data_surv$y_i_IPW_se[1,]),
                  mu_KML = t(data_surv$y_i_KML[1,]), 
                  mu_KML_se = t(data_surv$y_i_KML_se[1,]),
                  bench = t(data_surv$x_i))
  colnames(m_all) <- c("mu_IPW", "mu_IPW_se", "mu_KML", "mu_KML_se","bench")
  output = list()
  output$m_data <- m_data
  output$m_all <- m_all
  return(output)
}
# nonprobability samples
m1 <- model_func(diff_opt_1,"Opt1"); m1_data <- m1$m_data; m1_est <- m1$m_all
m2 <- model_func(diff_opt_2,"Opt2"); m2_data <- m2$m_data; m2_est <- m2$m_all
m3 <- model_func(diff_opt_3,"Opt3"); m3_data <- m3$m_data; m3_est <- m3$m_all
# probability samples
mp1 <- model_func(diff_prob_1,"Prob1"); mp1_data <- mp1$m_data; mp1_est <- mp1$m_all
mp2 <- model_func(diff_prob_2,"Prob2"); mp2_data <- mp2$m_data; mp2_est <- mp2$m_all
mp3 <- model_func(diff_prob_3,"Prob3"); mp3_data <- mp3$m_data; mp3_est <- mp3$m_all

# bias in the known case is the average of the differences between probability 
# and nonprobability surveys
bias_est <- ((m1_data$y_i_age_IPW - mp1_data$y_i_age_IPW) + # O1 - P1
  (m1_data$y_i_age_IPW - mp2_data$y_i_age_IPW) + # O1 - P2
    (m2_data$y_i_age_IPW - mp1_data$y_i_age_IPW) + # O2 - P1
    (m2_data$y_i_age_IPW - mp2_data$y_i_age_IPW) # O2 - P2
  )/4
# bias corrected CLW
m3_data$CLW_bc <- m3_data$y_i_age_IPW - bias_est

##############################################################################################
## Table of estimators in the known bias case
##############################################################################################

m_data <- tibble(age = mp3_data$age,
                 qs = mp3_data$qs,
                 bnchmark = mp3_data$x_i_age,
                 # ABS
                 ABS = mp3_data$y_i_age_IPW,
                 ABS_bias = round(abs(mp3_data$bias_i_age_IPW),3),
                 # CLW
                 CLW = m3_data$y_i_age_IPW,
                 CLW_bias = round(abs(m3_data$bias_i_age_IPW),3),
                 # bias corrected CLW
                 CLW_bc = m3_data$CLW_bc,
                 CLW_bc_bias = round(abs(bnchmark - CLW_bc),3), 
                 v1 = ((mp3_data$y_i_age_IPW_se)^2)*nrow(filter(dat_new, SAMPID == 3)),
                 v2 = ((m3_data$y_i_age_IPW_se)^2)*nrow(np3_weight),
                 # EV estimator
                 EV = ((bias_est^2 + v2)*(mp3_data$y_i_age_IPW) + 
                         (v1)*m3_data$y_i_age_IPW)/(bias_est^2 + v1 + v2),
                 EV_bias = round(abs(bnchmark - EV),3),
                 # comb estimator
                 comb = (v2/(v1 + v2))*(mp3_data$y_i_age_IPW) + 
                   (v1/(v1 + v2))*m3_data$CLW_bc,
                 comb_bias = round(abs(bnchmark - comb),3)
                 )
# MAE EV and Comb
m_data %>% summarise(ABS_bias = mean(ABS_bias),
                                       CLW_bias = mean(CLW_bias),
                                       CLW_bc_bias = mean(CLW_bc_bias),
                                       EV_bias = mean(EV_bias), 
                                       comb_bias = mean(comb_bias))
# group by age
m_data %>% group_by(age) %>% summarise(ABS_bias = mean(ABS_bias),
                                       CLW_bias = mean(CLW_bias),
                                       CLW_bc_bias = mean(CLW_bc_bias),
                                       EV_bias = mean(EV_bias), 
                                       comb_bias = mean(comb_bias))

print(paste0("ABS bias ",mean(m_data$ABS_bias)))  
print(paste0("Original CLW bias ",mean(m_data$CLW_bias)))
print(paste0("Bias corrected CLW bias ",mean(m_data$CLW_bc_bias)))
print(paste0("EV bias ",mean(m_data$EV_bias)))    
print(paste0("Comb bias ",mean(m_data$comb_bias)))    


mean(m_data$EV_bias); mean(m_data$comb_bias)

write.csv(m_data,"/Users/aditisen/Downloads/KML work/m_data.csv")
write.csv(m3_data,"/Users/aditisen/Downloads/KML work/m3_data.csv")


