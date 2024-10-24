
## Author : Aditi Sen
## Date : Oct 2024

############################################################################################################
# Topic: Comparison of composite estimator and estimator form probability sample
#         by sampling from ABS 3
############################################################################################################

library(pps); library(gtsummary)
 
# consider ps P3
P3_frame <- filter(data_KML, SAMPID == 3)
hist(P3_frame$WEIGHT)
# create stratum ID for each stratum based important variables
P3_frame$stratum_ID <- 
  ifelse(P3_frame$WEIGHT < 0.5,"0-<0.5",
         ifelse((P3_frame$WEIGHT >= 0.5) & (P3_frame$WEIGHT < 1), "0.5-<1",
                ifelse((P3_frame$WEIGHT >= 1) & (P3_frame$WEIGHT < 1.5), "1-<1.5",
                       ifelse((P3_frame$WEIGHT >= 1.5) & (P3_frame$WEIGHT < 2), 
                              "1.5-<2",">=2"))))
# check distribution (modify here by combining groups if sample size is too small)
P3_frame %>% group_by(stratum_ID) %>% summarise(ct=n(),pct=n()*100/nrow(P3_frame))
# population size 
N <- nrow(P3_frame) 
# create a numeric household ID as the function doesn't take non_numeric
P3_frame$ID <- seq(1:nrow(P3_frame))
# check how many strata are getting created
length(unique(P3_frame$stratum_ID))
# population frame needs to be ordered wrt to the stratum id before sampling
P3_frame_ord <- arrange(P3_frame, stratum_ID)
# vector of stratum population sizes
N_h <- P3_frame_ord %>% group_by(stratum_ID) %>% summarise(n = n())
N_h; sum(N_h$n)
# vector of weights
W_h <- (N_h$n)/N; sum(W_h)

### Universe Estimates for Raking Variables
P3_frame$AGECAT3BM <- as.factor(P3_frame$AGECAT3BM)
P3_frame$RACETHN <- as.factor(P3_frame$RACETHN)
P3_frame$EDUC3 <- as.factor(P3_frame$EDUC3)

# weighted count of age
age.dist  <- dplyr::count(x = P3_frame, AGECAT3BM, wt = WEIGHT)
colnames(age.dist) <- c("AGECAT3BM","Freq")
# weighted count of race
race.dist <- dplyr::count(x = P3_frame, RACETHN, wt = WEIGHT)
colnames(race.dist) <- c("RACETHN","Freq")
# weighted count of edu
edu.dist  <- dplyr::count(x = P3_frame, EDUC3, wt = WEIGHT)
colnames(edu.dist) <- c("EDUC3","Freq")

# function of sample size 
samp_gen <- function(n_opt){
  # vector of stratum sample sizes (round off to integer)
  n_h <- round(n_opt*W_h, 0)
  # function to output unit id to be dplyr::selected in the sample
  set.seed(138942087)
  sample_id <- ppssstrat(P3_frame_ord$ID, P3_frame_ord$stratum_ID, n_h)
  # final P3 sample of size n_opt
  P3_sample_data <- filter(P3_frame, ID %in% sample_id)
  ### Perform Raking
  Data.svy.unwgtd <- svydesign(ids=~1, weights = ~WEIGHT, data = P3_sample_data)
  data.svy.rake <- rake(Data.svy.unwgtd,list(~AGECAT3BM,~RACETHN,~EDUC3),
                         list(age.dist, race.dist, edu.dist))
  ### Weight Trimming
  #data.svy.rake.trim <- trimWeights(data.svy.rake, lower=.4, upper=5, strict=TRUE)
  P3_sample_data$wgt  <- weights(data.svy.rake)
  return(P3_sample_data)
}

samp100 <- samp_gen(100); samp200 <- samp_gen(200)
samp500 <- samp_gen(500); samp1000 <- samp_gen(1000)

nrow(samp100); sum(samp100$wgt); hist(samp100$wgt)
nrow(samp200); sum(samp200$wgt); hist(samp200$wgt)
nrow(samp500); sum(samp500$wgt); hist(samp500$wgt)
nrow(samp1000); sum(samp1000$wgt); hist(samp1000$wgt)

# There's no observation for worker's compensation for samples
samp1000 %>% group_by(AGECAT3BM,bm_WRKRSCOMP) %>% summarise(ct=n())
P3_frame %>% group_by(AGECAT3BM,bm_WRKRSCOMP) %>% summarise(ct=n())

results100 <- 
  survey::svydesign(~ 1, data = samp100, weights = ~ wgt) %>%
  tbl_svysummary(include = c(AGECAT3BM, RACETHN, EDUC3),
    statistic = list(all_categorical() ~ "{n_unweighted} ({p_unweighted}%) ({p}%)"))
results100

results500 <-
  survey::svydesign(~ 1, data = samp500, weights = ~ wgt) %>%
  tbl_svysummary(include = c(AGECAT3BM, RACETHN, EDUC3),
    statistic = list(all_categorical() ~ "{n_unweighted} ({p_unweighted}%) ({p}%)"))
results500

results1000 <-
  survey::svydesign(~ 1, data = samp1000, weights = ~ wgt) %>%
  tbl_svysummary(include = c(AGECAT3BM, RACETHN, EDUC3),
    statistic = list(all_categorical() ~ "{n_unweighted} ({p_unweighted}%) ({p}%)"))
results1000

results2 <-
  survey::svydesign(~ 1, data = P3_frame, weights = ~ WEIGHT) %>%
  tbl_svysummary(include = c(AGECAT3BM, RACETHN, EDUC3),
    statistic = list(all_categorical() ~ "{n_unweighted} ({p_unweighted})%  ({p}%)"))
results2

#=======================================================================================================
# check which columns are missing
nc100 <- setdiff(colnames(dummy_cols(data_KML, select_columns = as.character(bm_col_sub)) %>% dplyr::select(starts_with("bm"),-ends_with("99"))), 
        colnames(dummy_cols(samp100, select_columns = as.character(bm_col_sub)))) 
nc200 <-  setdiff(colnames(dummy_cols(data_KML, select_columns = as.character(bm_col_sub)) %>% dplyr::select(starts_with("bm"),-ends_with("99"))), 
                  colnames(dummy_cols(samp200, select_columns = as.character(bm_col_sub)))) 
nc500 <- setdiff(colnames(dummy_cols(data_KML, select_columns = as.character(bm_col_sub)) %>% dplyr::select(starts_with("bm"),-ends_with("99"))), 
                 colnames(dummy_cols(samp500, select_columns = as.character(bm_col_sub)))) 
nc1000 <- setdiff(colnames(dummy_cols(data_KML, select_columns = as.character(bm_col_sub)) %>% dplyr::select(starts_with("bm"),-ends_with("99"))), 
                  colnames(dummy_cols(samp1000, select_columns = as.character(bm_col_sub)))) 

# correct the benchmark dataframe
bchmk100 <- bchmk %>% dplyr::select(-all_of(nc100))
bchmk_age100 <- bchmk_age %>% dplyr::select(-all_of(nc100)) 
bchmk_race100 <- bchmk_race %>% dplyr::select(-all_of(nc100)) 
bchmk_edu100 <- bchmk_edu %>% dplyr::select(-all_of(nc100))

bchmk500 <- bchmk %>% dplyr::select(-all_of(nc500))
bchmk_age500 <- bchmk_age %>% dplyr::select(-all_of(nc500)) 
bchmk_race500 <- bchmk_race %>% dplyr::select(-all_of(nc500)) 
bchmk_edu500 <- bchmk_edu %>% dplyr::select(-all_of(nc500))

bchmk1000 <- bchmk %>% dplyr::select(-all_of(nc1000))
bchmk_age1000 <- bchmk_age %>% dplyr::select(-all_of(nc1000)) 
bchmk_race1000 <- bchmk_race %>% dplyr::select(-all_of(nc1000)) 
bchmk_edu1000 <- bchmk_edu %>% dplyr::select(-all_of(nc1000))

# Function :
IPW_func3 <- function(data_input, b_all, b_age, b_race, b_edu)
{
  bm_col_sub <- colnames(dplyr::select(data_input,starts_with("bm")))
  dat_bm <- dummy_cols(data_input, select_columns = as.character(bm_col_sub))
  new_col<- setdiff(colnames(dat_bm), colnames(dplyr::select(dat_bm,ends_with(c("99")))))
  dat_bm <- dat_bm %>% dplyr::select(all_of(new_col))
  data_surv <- dat_bm[,-c(9:ncol(data_input))]

  # calibration factor to multiply weight with
  adj_fac <- N_B/nrow(data_surv); data_surv$KML_weight <- (data_input$wgt)*(adj_fac)
  data_surv$new_weight <- (data_input$wgt)*(adj_fac)
  ############################################################################################
  # define function to multiply weight with dummy variables and sum overall and by category
  all_sum <- function(input_weight, input_category){
    #input_weight = "KML_weight"; input_category = "AGECAT3BM"
    # define design
    design_p1 <- survey::svydesign(ids = ~ CASEID, data = data_surv, weights =~ get(input_weight)) 
    bm_p1 <- survey::svymean( data_surv[9:(ncol(data_surv)-2)], design = design_p1)
    bm_p11 <- as.data.frame(bm_p1)  
    col_sum <- as.data.frame(t(bm_p11$mean)); col_se <- as.data.frame(t(bm_p11$SE))
    colnames(col_sum) <- colnames(col_se) <- 
      c('bm_INSURANCE_0','bm_INSURANCE_1','bm_BLOODPR_0','bm_BLOODPR_1','bm_HHPARENT_0',
               'bm_HHPARENT_1','bm_FOODALLERGY_0','bm_FOODALLERGY_1','bm_JOBLASTYR_0','bm_JOBLASTYR_1',
               'bm_RETACCT_0','bm_RETACCT_1','bm_UNEMPCOMP_0','bm_UNEMPCOMP_1',
               'bm_WRKRSCOMP_0','bm_WRKRSCOMP_1','bm_FOODSTAMPS_0','bm_FOODSTAMPS_1','bm_SOCSEC_0',
               'bm_SOCSEC_1','bm_UNION_0','bm_UNION_1','bm_CITIZEN_0','bm_CITIZEN_1')
    
    ############################################################################################
    # weighted and unweighted sum by cat
    func_cat <- function(cat_level){
      #cat_level <- 1 
      # weighted sum
      data_cat <- filter(data_surv, get(input_category) == cat_level)
      # define design
      design_p1c <- survey::svydesign(ids = ~ CASEID, data = data_cat, weights =~get(input_weight)) 
      # mean and SE for all variables
      bm_p1c <- survey::svymean( data_cat[9:(ncol(data_cat)-2)], design = design_p1c)
      bm_p11c <- as.data.frame(bm_p1c)  
      
      col_sum_cat <- as.data.frame(t(bm_p11c$mean)); col_se_cat <- as.data.frame(t(bm_p11c$SE))
      colnames(col_sum_cat) <- colnames(col_se_cat) <- 
        c('bm_INSURANCE_0','bm_INSURANCE_1','bm_BLOODPR_0','bm_BLOODPR_1','bm_HHPARENT_0',
                 'bm_HHPARENT_1','bm_FOODALLERGY_0','bm_FOODALLERGY_1','bm_JOBLASTYR_0','bm_JOBLASTYR_1',
                 'bm_RETACCT_0','bm_RETACCT_1','bm_UNEMPCOMP_0','bm_UNEMPCOMP_1',
                 'bm_WRKRSCOMP_0','bm_WRKRSCOMP_1','bm_FOODSTAMPS_0','bm_FOODSTAMPS_1','bm_SOCSEC_0',
                 'bm_SOCSEC_1','bm_UNION_0','bm_UNION_1','bm_CITIZEN_0','bm_CITIZEN_1')
      
      # unweighted sum
      data_cat_uw <- filter(data_surv, get(input_category) == cat_level) %>% dplyr::select(starts_with("bm"))
      col_sum_cat_uw <- as.data.frame(t(colSums(data_cat_uw,na.rm = TRUE)/nrow(data_cat_uw)))
      
      output_cat <- list()
      output_cat$col_sum_cat <- col_sum_cat; output_cat$col_se_cat <- col_se_cat
      output_cat$col_sum_cat_uw <- col_sum_cat_uw; #output_cat$col_se_cat_uw <- col_se_cat_uw
      return(output_cat)
    }
    ############################################################################################
    cat1 <- func_cat(1); cat2 <- func_cat(2); cat3 <- func_cat(3) 
    # weighted sum by category
    col_sum_cat1 <- cat1$col_sum_cat; col_sum_cat2 <- cat2$col_sum_cat; col_sum_cat3 <- cat3$col_sum_cat
    col_se_cat1 <- cat1$col_se_cat; col_se_cat2 <- cat2$col_se_cat; col_se_cat3 <- cat3$col_se_cat
    # unweighted sum by category
    col_sum_cat1_uw <- cat1$col_sum_cat_uw; col_sum_cat2_uw <- cat2$col_sum_cat_uw; col_sum_cat3_uw <- cat3$col_sum_cat_uw
    
    # combine all categories
    data_cat <- as.data.frame(rbind(col_sum_cat1,col_sum_cat2,col_sum_cat3))
    data_cat_se <- as.data.frame(rbind(col_se_cat1,col_se_cat2,col_se_cat3))
    data_cat_uw <- as.data.frame(rbind(col_sum_cat1_uw,col_sum_cat2_uw,col_sum_cat3_uw))
    
    output_all_sum <- list()
    output_all_sum$col_sum <- col_sum; output_all_sum$col_se <- col_se
    output_all_sum$data_cat <- data_cat; output_all_sum$data_cat_se <- data_cat_se
    output_all_sum$data_cat_uw <- data_cat_uw
    
    return(output_all_sum)
  }
  ############################################################################################
  f_age_KML <- all_sum("KML_weight","AGECAT3BM") # age
  f_race_KML <- all_sum("KML_weight","RACETHN") # race
  f_edu_KML <- all_sum("KML_weight","EDUC3") # edu
  
  f_age_IPW <- all_sum("new_weight","AGECAT3BM") # age
  f_race_IPW <- all_sum("new_weight","RACETHN") # race
  f_edu_IPW <- all_sum("new_weight","EDUC3") # edu
  ############################################################################################
  # for probability sample use given weight as KML, there is no IPW
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
  
  ############################################################################################
  # unweighted
  col_sum_uw <- as.data.frame(t(colSums(data.frame(data_surv[,-c(1:8,ncol(data_surv)-1,ncol(data_surv))]),
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
    yc_bar_age_uw <- yc_bar_age_KML <-yc_bar_age_KML_se <- yc_bar_age_IPW <- yc_bar_age_IPW_se <-
    yc_bar_race_uw <- yc_bar_race_KML <- yc_bar_race_KML_se <- yc_bar_race_IPW <- yc_bar_race_IPW_se <- 
    yc_bar_edu_uw <- yc_bar_edu_KML <- yc_bar_edu_KML_se <- yc_bar_edu_IPW <- yc_bar_edu_IPW_se <- 
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
  
  # data required for modelling
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
    yc_bar_KML_se[1,which(bm_col_sub==i)] <- rowMeans(dplyr::select(col_se_1_KML,starts_with(i)))
    yc_bar_IPW_se[1,which(bm_col_sub==i)] <- rowMeans(dplyr::select(col_se_1_IPW,starts_with(i)))
    ############################################################################################
    # calculate difference by age and summarize at variable level
    diff_age_uw <- dplyr::select(col_sum_age_uw,starts_with(i)) - dplyr::select(b_age[,-1],starts_with(i))
    diff_age_KML <- dplyr::select(col_sum_age_KML,starts_with(i)) - dplyr::select(b_age[,-1],starts_with(i))
    diff_age_IPW <- dplyr::select(col_sum_age_IPW,starts_with(i)) - dplyr::select(b_age[,-1],starts_with(i))
    
    yc_bar_age_uw[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_age_uw)))/k
    yc_bar_age_KML[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_age_KML)))/k
    yc_bar_age_IPW[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_age_IPW)))/k
    
    # summarize standard error at variable level 
    yc_bar_age_KML_se[c(1:3),which(bm_col_sub==i)] <- rowMeans(dplyr::select(col_se_age_KML,starts_with(i)))
    yc_bar_age_IPW_se[c(1:3),which(bm_col_sub==i)] <- rowMeans(dplyr::select(col_se_age_IPW,starts_with(i)))
    ############################################################################################
    # calculate difference by race and summarize at variable level  
    diff_race_uw <- dplyr::select(col_sum_race_uw,starts_with(i)) - dplyr::select(b_race[,-1],starts_with(i))
    diff_race_KML <- dplyr::select(col_sum_race_KML,starts_with(i)) - dplyr::select(b_race[,-1],starts_with(i))
    diff_race_IPW <- dplyr::select(col_sum_race_IPW,starts_with(i)) - dplyr::select(b_race[,-1],starts_with(i))
    
    yc_bar_race_uw[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_race_uw)))/k
    yc_bar_race_KML[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_race_KML)))/k
    yc_bar_race_IPW[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_race_IPW)))/k
    
    # summarize standard error at variable level 
    yc_bar_race_KML_se[c(1:3),which(bm_col_sub==i)] <- rowMeans(dplyr::select(col_se_race_KML,starts_with(i)))
    yc_bar_race_IPW_se[c(1:3),which(bm_col_sub==i)] <- rowMeans(dplyr::select(col_se_race_IPW,starts_with(i)))
    ############################################################################################
    # calculate difference by education and summarize at variable level 
    diff_edu_uw <- dplyr::select(col_sum_edu_uw,starts_with(i)) - dplyr::select(b_edu[,-1],starts_with(i))
    diff_edu_KML <- dplyr::select(col_sum_edu_KML,starts_with(i)) - dplyr::select(b_edu[,-1],starts_with(i))
    diff_edu_IPW <- dplyr::select(col_sum_edu_IPW,starts_with(i)) - dplyr::select(b_edu[,-1],starts_with(i))
    
    yc_bar_edu_uw[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_edu_uw)))/k
    yc_bar_edu_KML[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_edu_KML)))/k
    yc_bar_edu_IPW[c(1:3),which(bm_col_sub==i)] <- rowSums(as.matrix(abs(diff_edu_IPW)))/k
    
    # summarize standard error at variable level 
    yc_bar_edu_KML_se[c(1:3),which(bm_col_sub==i)] <- rowMeans(dplyr::select(col_se_edu_KML,starts_with(i)))
    yc_bar_edu_IPW_se[c(1:3),which(bm_col_sub==i)] <- rowMeans(dplyr::select(col_se_edu_IPW,starts_with(i)))
    ############################################################################################
    # store bias for modelling
    if(k == 2){
      y_i_uw <- cbind(y_i_uw, dplyr::select(col_sum_1_uw, starts_with(i)))
      y_i_KML <- cbind(y_i_KML, dplyr::select(col_sum_1_KML, starts_with(i)))
      y_i_KML_se <- cbind(y_i_KML_se, dplyr::select(col_se_1_KML, starts_with(i)))
      y_i_IPW <- cbind(y_i_IPW, dplyr::select(col_sum_1_IPW, starts_with(i)))
      y_i_IPW_se <- cbind(y_i_IPW_se, dplyr::select(col_se_1_IPW, starts_with(i)))
      x_i <- cbind(x_i,select(b_all, starts_with(i)))
      
      bias_i_uw <- cbind(bias_i_uw, diff_uw)
      bias_i_KML <- cbind(bias_i_KML, diff_KML)
      bias_i_IPW <- cbind(bias_i_IPW, diff_IPW)
      
      y_i_age_uw <- cbind(y_i_age_uw, dplyr::select(col_sum_age_uw, starts_with(i)))
      y_i_age_KML <- cbind(y_i_age_KML, dplyr::select(col_sum_age_KML, starts_with(i)))
      y_i_age_KML_se <- cbind(y_i_age_KML_se, dplyr::select(col_se_age_KML, starts_with(i)))
      y_i_age_IPW <- cbind(y_i_age_IPW, dplyr::select(col_sum_age_IPW, starts_with(i)))
      y_i_age_IPW_se <- cbind(y_i_age_IPW_se, dplyr::select(col_se_age_IPW, starts_with(i)))
      
      x_i_age <- cbind(x_i_age, dplyr::select(b_age[,-1], starts_with(i)))
      
      bias_i_age_uw <- cbind(bias_i_age_uw, diff_age_uw)
      bias_i_age_KML <- cbind(bias_i_age_KML, diff_age_KML)
      bias_i_age_IPW <- cbind(bias_i_age_IPW, diff_age_IPW)
      
      y_i_race_uw <- cbind(y_i_race_uw, dplyr::select(col_sum_race_uw, starts_with(i)))
      y_i_race_KML <- cbind(y_i_race_KML, dplyr::select(col_sum_race_KML, starts_with(i)))
      y_i_race_KML_se <- cbind(y_i_race_KML_se, dplyr::select(col_se_race_KML, starts_with(i)))
      y_i_race_IPW <- cbind(y_i_race_IPW, dplyr::select(col_sum_race_IPW, starts_with(i)))
      y_i_race_IPW_se <- cbind(y_i_race_IPW_se, dplyr::select(col_se_race_IPW, starts_with(i)))
      
      x_i_race <- cbind(x_i_race,select(b_race[,-1],starts_with(i)))
      
      bias_i_race_uw <- cbind(bias_i_race_uw, diff_race_uw)
      bias_i_race_KML <- cbind(bias_i_race_KML, diff_race_KML)
      bias_i_race_IPW <- cbind(bias_i_race_IPW, diff_race_IPW)
      
      y_i_edu_uw <- cbind(y_i_edu_uw, dplyr::select(col_sum_edu_uw, starts_with(i)))
      y_i_edu_KML <- cbind(y_i_edu_KML, dplyr::select(col_sum_edu_KML, starts_with(i)))
      y_i_edu_KML_se <- cbind(y_i_edu_KML_se, dplyr::select(col_se_edu_KML, starts_with(i)))
      y_i_edu_IPW <- cbind(y_i_edu_IPW, dplyr::select(col_sum_edu_IPW, starts_with(i)))
      y_i_edu_IPW_se <- cbind(y_i_edu_IPW_se, dplyr::select(col_se_edu_IPW, starts_with(i)))
      
      x_i_edu <- cbind(x_i_edu,select(b_edu[,-1],starts_with(i)))
      
      bias_i_edu_uw <- cbind(bias_i_edu_uw, diff_edu_uw)
      bias_i_edu_KML <- cbind(bias_i_edu_KML, diff_edu_KML)
      bias_i_edu_IPW <- cbind(bias_i_edu_IPW, diff_edu_IPW)
      
    }
  }
  
  out_put_val <- list()
  
  # calculate mean of difference for 12 benchmark variables  
  out_put_val$mean_uw <- rowSums(yc_bar_uw)/12
  out_put_val$mean_KML <- rowSums(yc_bar_KML)/12
  out_put_val$mean_IPW <- rowSums(yc_bar_IPW)/12
  
  # how to summarize SE at dataset level????
  out_put_val$se_uw <-  sd(yc_bar_uw)/5 # need to correct this
  out_put_val$se_KML <-  rowMeans(yc_bar_KML_se) 
  out_put_val$se_IPW <-  rowMeans(yc_bar_IPW_se) 
  
  # group by age
  out_put_val$mean_age_uw <- rowSums(yc_bar_age_uw)/12
  out_put_val$mean_age_KML <- rowSums(yc_bar_age_KML)/12
  out_put_val$mean_age_IPW <- rowSums(yc_bar_age_IPW)/12
  
  out_put_val$se_age_uw <- c(sd(yc_bar_age_uw[1,])/5,sd(yc_bar_age_uw[2,])/5,sd(yc_bar_age_uw[3,])/5) # need to correct this
  out_put_val$se_age_KML <- rowMeans(yc_bar_age_KML_se)
  out_put_val$se_age_IPW <- rowMeans(yc_bar_age_IPW_se) 
  
  # group by race
  out_put_val$mean_race_uw <- rowSums(yc_bar_race_uw)/12
  out_put_val$mean_race_KML <- rowSums(yc_bar_race_KML)/12
  out_put_val$mean_race_IPW <- rowSums(yc_bar_race_IPW)/12
  
  out_put_val$se_race_uw <-  c(sd(yc_bar_race_uw[1,])/5,sd(yc_bar_race_uw[2,])/5,sd(yc_bar_race_uw[3,])/5) # need to correct this
  out_put_val$se_race_KML <-  rowMeans(yc_bar_race_KML_se)
  out_put_val$se_race_IPW <-  rowMeans(yc_bar_race_IPW_se)  
  
  # group by education
  out_put_val$mean_edu_uw <- rowSums(yc_bar_edu_uw)/12
  out_put_val$mean_edu_KML <- rowSums(yc_bar_edu_KML)/12
  out_put_val$mean_edu_IPW <- rowSums(yc_bar_edu_IPW)/12
  
  out_put_val$se_edu_uw <-  c(sd(yc_bar_edu_uw[1,])/5,sd(yc_bar_edu_uw[2,])/5,sd(yc_bar_edu_uw[3,])/5) # need to correct this
  out_put_val$se_edu_KML <-  rowMeans(yc_bar_edu_KML_se) 
  out_put_val$se_edu_IPW <-  rowMeans(yc_bar_edu_IPW_se)  
  
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
ABS_3_sample100 <- IPW_func3(data_input = samp100, b_all = bchmk100, 
                             b_age = bchmk_age100, b_race = bchmk_race100, 
                             b_edu = bchmk_edu100)
ABS_3_sample500 <- IPW_func3(data_input = samp500, b_all = bchmk500, 
                             b_age = bchmk_age500, b_race = bchmk_race500, 
                             b_edu = bchmk_edu500)
ABS_3_sample1000 <- IPW_func3(data_input = samp1000, b_all = bchmk1000, 
                              b_age = bchmk_age1000, b_race = bchmk_race1000, 
                              b_edu = bchmk_edu1000)

ABS_3_data100 <- model_func(ABS_3_sample100,"ABS 3 sample")
ABS_3_data500 <- model_func(ABS_3_sample500,"ABS 3 sample")
ABS_3_data1000 <- model_func(ABS_3_sample1000,"ABS 3 sample")

ABS_3_data_1_100 <- ABS_3_data100$m_data
ABS_3_data_1_500 <- ABS_3_data500$m_data
ABS_3_data_1_1000 <- ABS_3_data1000$m_data

P3_sample_data <- tibble(age = mp3_data$age,
                 qs = mp3_data$qs,
                 bnchmark = mp3_data$x_i_age,
                 # ABS
                 ABS_100 = ABS_3_data_1_100$y_i_age_IPW,
                 ABS_100_bias = round(abs(ABS_3_data_1_100$bias_i_age_IPW),3),
                 ABS_500 = ABS_3_data_1_500$y_i_age_IPW,
                 ABS_500_bias = round(abs(ABS_3_data_1_500$bias_i_age_IPW),3),
                 ABS_1000 = ABS_3_data_1_1000$y_i_age_IPW,
                 ABS_1000_bias = round(abs(ABS_3_data_1_1000$bias_i_age_IPW),3),
                 # CLW
                 CLW = m3_data$y_i_age_IPW,
                 CLW_bias = round(abs(m3_data$bias_i_age_IPW),3),
                 # bias corrected CLW
                 CLW_bc = m3_data$CLW_bc,
                 CLW_bc_bias = round(abs(bnchmark - CLW_bc),3), 
                 v1_100 = ((ABS_3_data_1_100$y_i_age_IPW_se)^2)*100,
                 v1_500 = ((ABS_3_data_1_500$y_i_age_IPW_se)^2)*500,
                 v1_1000 = ((ABS_3_data_1_1000$y_i_age_IPW_se)^2)*1000,
                 v2 = ((m3_data$y_i_age_IPW_se)^2)*nrow(np3_weight),
                 # comb estimator
                 comb_100 = (v2/(v1_100 + v2))*ABS_100 + (v1_100/(v1_100 + v2))*CLW_bc,
                 comb_500 = (v2/(v1_500 + v2))*ABS_500 + (v1_500/(v1_500 + v2))*CLW_bc,
                 comb_1000 = (v2/(v1_1000 + v2))*ABS_1000 + (v1_1000/(v1_1000 + v2))*CLW_bc,
                 comb_100_bias = round(abs(bnchmark - comb_100),3),
                 comb_500_bias = round(abs(bnchmark - comb_500),3),
                 comb_1000_bias = round(abs(bnchmark - comb_1000),3)
)
# MAE EV and Comb
P3_sample_data %>% summarise(ABS_100_bias = mean(ABS_100_bias),
                     ABS_500_bias = mean(ABS_500_bias),
                     ABS_1000_bias = mean(ABS_1000_bias),
                     comb_100_bias = mean(comb_100_bias),
                     comb_500_bias = mean(comb_500_bias),
                     comb_1000_bias = mean(comb_1000_bias)
                     )

write.csv(ABS_3_data_1_100,"/Users/aditisen/Downloads/KML work/ABS_3_data_1_100.csv")
write.csv(ABS_3_data_1_500,"/Users/aditisen/Downloads/KML work/ABS_3_data_1_500.csv")
write.csv(ABS_3_data_1_1000,"/Users/aditisen/Downloads/KML work/ABS_3_data_1_1000.csv")




