library(plyr)
library(dplyr)
library(mice)
library(miceadds)

# To normalize numeric data
z_transform = function(impList_item){
  
  df_num = impList_item[unlist(lapply(impList_item, is.numeric))]
  df_fac = impList_item[unlist(lapply(impList_item, is.factor))]
  
  means = lapply(df_num, mean)
  sds = lapply(df_num, sd)
  
  df_num_transf = data.frame(lapply(df_num, scale, center = TRUE, scale = TRUE))
  df_new = cbind(df_num_transf, df_fac)
  
  return(list(df_new, means, sds))
}

# Imputing data
impute = function(imp_df){
  numeric_cols = c('ventricles','hippocampus', 'icv', 'educ',
                   'age','abeta', 'tau', 'bmi', 'fdg', 'mmse') 
  factor_cols = c('dx_bl', 'sex', 'alc_cat', 'smok_cat','hyp_ind', 'apoe4', 'cardio')
  
  imp_df_sel = cbind(imp_df[, numeric_cols],
                     colwise(factor, ordered=TRUE)(imp_df[factor_cols[!factor_cols%in%c('sex', 'apoe4')]]),
                     sex=factor(imp_df$sex), apoe4 = factor(imp_df$apoe4) #Sex should not be an ordered factor.
  )
  
  impData = mice(imp_df_sel, m=1, seed=42)
  imputed_df = complete(impData, action = 1) #Get imputed dataset
  
  # Z-transform imputed data
  imputed_df_t = z_transform(imputed_df)
  imp_t_df = imputed_df_t[[1]] #get transformed values
  
  #Save means and sds of all imputed variables. 
  #The function returns the means/sds before the z-transformation
  imp_t_means = imputed_df_t[[2]]
  imp_t_sds = imputed_df_t[[3]]
  
  #Revalue dx_bl for every impList item
  imp_t_df$dx_bl = revalue(imp_t_df$dx_bl, c("1"="0", "2"="1", "3"="1", "4"="1"))
  
  return(list(imp_t_df, imp_t_means, imp_t_sds))
}

# Extract SEM parameters
extract_semparams = function(parTable_df){
  lhs_betas = parTable_df$lhs[parTable_df$op=='~']
  
  model_strs = c()
  variances = c()
  ths = c()
  ths_df = data.frame()
  
  for(endov in unique(lhs_betas)){
    df = parTable_df[parTable_df$lhs==endov,]
    model_str = ''
    
    for(row in 1:nrow(df)){
      rhs = df[row, 'rhs']
      lhs = df[row, 'lhs']
      op = df[row, 'op']
      param = df[row, 'est']  #Take the mean of all fits  --->>13/4/21 changed to est instead of 'Mean' because I am currently not working with mean imputations
      
      if(row==1){
        if(op=='~'){
          model_str = paste(model_str, param, '*', rhs)}}
      if(row>1 & op=='~'){
        model_str = paste(model_str, '+ (', param, ') *', rhs )}
      if(op=='~1'){
        model_str = paste(model_str, '+ (', param, ')' )}
      
      #save the variance into the variance vector
      if(op=='~~'){
        if(rhs==lhs){
          variances = c(variances, param)}
      }
      #variances = c(variances, ifelse(rhs==lhs, param, 1))
      #variances = c(variances, param)}
      
      #Attach threshold
      ths = ifelse(op=='|', param, NA)
      ths_df = rbind(ths_df, data.frame(endov=endov, ths=ths))
    }
    model_strs = c(model_strs, model_str)
  }
  
  genData_df = data.frame(endov=unique(lhs_betas), variances, model_strs)
  
  return(list(genData_df, na.omit(ths_df)))
}

# Generate endogenous variables
generate_endogenous = function(genData_df, ths_df = ths_df, a, ap, s, n_samples, numeric_cols=numeric_cols){
  for(ii in genData_df$endov){
    if(ii%in%numeric_cols){
      variance = genData_df$variances[genData_df$endov==ii]
      if(variance<0){variance=variance+1}
      std_dev = sqrt(variance)
      
      rand_noise = rnorm(n_new_samples, mean = 0, sd=std_dev)
      assign(ii, eval(parse(text=genData_df$model_strs[genData_df$endov==ii])) + rand_noise)
      
      
    }else{ #If categorical variable
      variance = genData_df$variances[genData_df$endov==ii]
      std_dev = sqrt(variance)
      
      rand_noise = rnorm(n_new_samples, mean = 0, sd=std_dev)
      sim0 = eval(parse(text=genData_df$model_strs[genData_df$endov==ii])) + rand_noise
      ths = ths_df$ths[which(ths_df$endov == ii)]
      assign(ii, mapply(FUN=function(sim0, ths){ifelse(sim0>ths,1,0)}, sim0, ths))
    }
  }
  
  all_sim = data.frame(a=a, ap=ap, s = s, dx= factor(dx), me=me, hi=hi, educ = educ, ab=ab, ve = ve,
                       bmi = bmi, hyp = factor(hyp), alc=factor(alc), smo = factor(smo), tau = tau,
                       icv = icv, fdg = fdg, cardio = factor(cardio)) 
  
  return(all_sim)
}


# Initialize report
init_report = function(){
  
  report <- data.frame(sim = NA_real_,
                       ICI0_logistic_all = NA_real_,
                       ICI0_logistic_par = NA_real_,
                       ICI0_logistic_child = NA_real_,
                       
                       ICI1_logistic_all = NA_real_,
                       ICI1_logistic_par = NA_real_,
                       ICI1_logistic_child = NA_real_,
                       
                       ICI2_logistic_all = NA_real_,
                       ICI2_logistic_par = NA_real_,
                       ICI2_logistic_child = NA_real_,
                       
                       ICI3_logistic_all = NA_real_,
                       ICI3_logistic_par = NA_real_,
                       ICI3_logistic_child = NA_real_,
                       
                       ICI4_logistic_all = NA_real_,
                       ICI4_logistic_par = NA_real_,
                       ICI4_logistic_child = NA_real_,
                       
                       ICI0_elastnet_all = NA_real_,
                       ICI0_elastnet_par = NA_real_,
                       ICI0_elastnet_child = NA_real_,
                       
                       ICI1_elastnet_all = NA_real_,
                       ICI1_elastnet_par = NA_real_,
                       ICI1_elastnet_child = NA_real_,
                       
                       ICI2_elastnet_all = NA_real_,
                       ICI2_elastnet_par = NA_real_,
                       ICI2_elastnet_child = NA_real_,
                       
                       ICI3_elastnet_all = NA_real_,
                       ICI3_elastnet_par = NA_real_,
                       ICI3_elastnet_child = NA_real_,
                       
                       ICI4_elastnet_all = NA_real_,
                       ICI4_elastnet_par = NA_real_,
                       ICI4_elastnet_child = NA_real_,
                       
                       ICI0_rf_all = NA_real_,
                       ICI0_rf_par = NA_real_,
                       ICI0_rf_child = NA_real_,
                       
                       ICI1_rf_all = NA_real_,
                       ICI1_rf_par = NA_real_,
                       ICI1_rf_child = NA_real_,
                       ICI2_rf_all = NA_real_,
                       ICI2_rf_par = NA_real_,
                       ICI2_rf_child = NA_real_,
                       ICI3_rf_all = NA_real_,
                       ICI3_rf_par = NA_real_,
                       ICI3_rf_child = NA_real_,
                       ICI4_rf_all = NA_real_,
                       ICI4_rf_par = NA_real_,
                       ICI4_rf_child = NA_real_,
                       
                       ICI0_gbm_all = NA_real_,
                       ICI0_gbm_par = NA_real_,
                       ICI0_gbm_child = NA_real_,
                       ICI1_gbm_all = NA_real_,
                       ICI1_gbm_par = NA_real_,
                       ICI1_gbm_child = NA_real_,
                       ICI2_gbm_all = NA_real_,
                       ICI2_gbm_par = NA_real_,
                       ICI2_gbm_child = NA_real_,
                       ICI3_gbm_all = NA_real_,
                       ICI3_gbm_par = NA_real_,
                       ICI3_gbm_child = NA_real_,
                       ICI4_gbm_all = NA_real_,
                       ICI4_gbm_par = NA_real_,
                       ICI4_gbm_child = NA_real_,
                       
                       Brier0_logistic_all = NA_real_,
                       Brier0_logistic_par = NA_real_,
                       Brier0_logistic_child = NA_real_,
                       Brier1_logistic_all = NA_real_,
                       Brier1_logistic_par = NA_real_,
                       Brier1_logistic_child = NA_real_,
                       Brier2_logistic_all = NA_real_,
                       Brier2_logistic_par = NA_real_,
                       Brier2_logistic_child = NA_real_,
                       Brier3_logistic_all = NA_real_,
                       Brier3_logistic_par = NA_real_,
                       Brier3_logistic_child = NA_real_,
                       Brier4_logistic_all = NA_real_,
                       Brier4_logistic_par = NA_real_,
                       Brier4_logistic_child = NA_real_,
                       
                       Brier0_elastnet_all = NA_real_,
                       Brier0_elastnet_par = NA_real_,
                       Brier0_elastnet_child = NA_real_,
                       
                       Brier1_elastnet_all = NA_real_,
                       Brier1_elastnet_par = NA_real_,
                       Brier1_elastnet_child = NA_real_,
                       
                       Brier2_elastnet_all = NA_real_,
                       Brier2_elastnet_par = NA_real_,
                       Brier2_elastnet_child = NA_real_,
                       
                       Brier3_elastnet_all = NA_real_,
                       Brier3_elastnet_par = NA_real_,
                       Brier3_elastnet_child = NA_real_,
                       Brier4_elastnet_all = NA_real_,
                       Brier4_elastnet_par = NA_real_,
                       Brier4_elastnet_child = NA_real_,
                       
                       Brier0_rf_all = NA_real_,
                       Brier0_rf_par = NA_real_,
                       Brier0_rf_child = NA_real_,
                       Brier1_rf_all = NA_real_,
                       Brier1_rf_par = NA_real_,
                       Brier1_rf_child = NA_real_,
                       Brier2_rf_all = NA_real_,
                       Brier2_rf_par = NA_real_,
                       Brier2_rf_child = NA_real_,
                       Brier3_rf_all = NA_real_,
                       Brier3_rf_par = NA_real_,
                       Brier3_rf_child = NA_real_,
                       Brier4_rf_all = NA_real_,
                       Brier4_rf_par = NA_real_,
                       Brier4_rf_child = NA_real_,
                       
                       Brier0_gbm_all = NA_real_,
                       Brier0_gbm_par = NA_real_,
                       Brier0_gbm_child = NA_real_,
                      
                       Brier1_gbm_all = NA_real_,
                       Brier1_gbm_par = NA_real_,
                       Brier1_gbm_child = NA_real_,
                       
                       Brier2_gbm_all = NA_real_,
                       Brier2_gbm_par = NA_real_,
                       Brier2_gbm_child = NA_real_,
                       
                       Brier3_gbm_all = NA_real_,
                       Brier3_gbm_par = NA_real_,
                       Brier3_gbm_child = NA_real_,
                       Brier4_gbm_all = NA_real_,
                       Brier4_gbm_par = NA_real_,
                       Brier4_gbm_child = NA_real_,
                        
                       no_lassopred_all = NA_real_,
                       no_lassopred_par = NA_real_,
                       no_lassopred_child = NA_real_,
                       lambda_all = NA_real_,
                       lambda_par = NA_real_,
                       lambda_child = NA_real_,
                       
                       sum_dev_num = NA_real_,
                       sum_dev_cat = NA_real_,
                       sum_test_num = NA_real_,
                       sum_test_cat = NA_real_,
                       sum_intap_num = NA_real_,
                       sum_intap_cat = NA_real_,
                       sum_ints_num = NA_real_,
                       sum_ints_cat = NA_real_,
                       sum_inttau_num = NA_real_,
                       sum_inttau_cat = NA_real_,
                       
                       elast_all_a = NA_real_,
                       elast_all_ap = NA_real_,
                       elast_all_s = NA_real_,
                       elast_all_educ = NA_real_,
                       elast_all_bmi = NA_real_,
                       elast_all_hyp = NA_real_,
                       elast_all_alc = NA_real_,
                       elast_all_smo = NA_real_,
                       elast_all_cardio = NA_real_,
                       elast_all_tau = NA_real_,
                       elast_all_ab = NA_real_,
                       elast_all_hi = NA_real_,
                       elast_all_ve = NA_real_,
                       elast_all_icv = NA_real_,
                       elast_all_fdg = NA_real_,
                       elast_all_me = NA_real_
                       )
  
  return(report)
}


grid_search_elastnet = function(dev_df, outcome='dx'){
  #Get best lambda for Lasso-Regression (alpha=1)
  library(dplyr)
  library(glmnet)
  
  X = data.matrix(dplyr::select(dev_df,-outcome))
  y = data.matrix(dev_df[, outcome])
    
  cv <- cv.glmnet(X, y, family = "binomial", nfold = 10, 
                    type.measure = "deviance", alpha = 1)
    
  lambda.1se = cv$lambda[cv$lambda == cv$lambda.1se]
  lambda.min = cv$lambda[cv$lambda == cv$lambda.min]
  
  return(lambda.min)
}


# Calculate calibration metrics 
ICI <-function(outcome, prediction){
  P.calibrate <-predict(loess(outcome~prediction))
  (ICI <-mean(abs(P.calibrate-prediction)))}

brierscore = function(out, pred){
  bins = unique(c(0, as.vector(quantile(pred, probs=seq(0, 1, 0.1))), 1))
  brier = data.frame(BrierDecomp(p=pred, y=out, bins = bins, bias.corrected=TRUE))
  reliability = brier$REL[1]
  return(reliability)
}


# Predict functions
predict_logistic = function(logistic_mod, test_sim_t=test_sim_t, 
                            synth_df_inta=synth_df_inta, 
                            synth_df_intap=synth_df_intap, synth_df_ints=synth_df_ints, synth_df_inttau=inttau_num){
  #Note: It is not necessary to pass the prediction 'set' here, as in rf or gbm predict-functions. The model that is passed, includes this info already
  
  pred_0 = predict(logistic_mod, type = "response", newdata = test_sim_t)
  pred_1 = predict(logistic_mod, type = "response", newdata = synth_df_inta)
  pred_2 = predict(logistic_mod, type = "response", newdata = synth_df_intap)
  pred_3 = predict(logistic_mod, type = "response", newdata = synth_df_ints)
  pred_4 = predict(logistic_mod, type = "response", newdata = synth_df_inttau)
  
  preds = data.frame(pred_0, pred_1, pred_2, pred_3, pred_4)
  return(preds)
}

predict_elastic = function(elastic_mod, set, test_sim_t=test_sim_t, 
                            synth_df_inta=synth_df_inta, 
                            synth_df_intap=synth_df_intap, synth_df_ints=synth_df_ints, 
                            synth_df_inttau=inttau_num){
  s = elastic_mod$lambda
  
  pred_0 = c(predict(elastic_mod, type = "response", s=s, newx = data.matrix(test_sim_t[, set])))
  pred_1 = c(predict(elastic_mod, type = "response", s=s, newx = data.matrix(synth_df_inta[, set])))
  pred_2 = c(predict(elastic_mod, type = "response", s=s, newx = data.matrix(synth_df_intap[, set])))
  pred_3 = c(predict(elastic_mod, type = "response", s=s, newx = data.matrix(synth_df_ints[, set])))
  pred_4 = c(predict(elastic_mod, type = "response", s=s, newx = data.matrix(synth_df_inttau[, set])))
  
  preds = data.frame(pred_0, pred_1, pred_2, pred_3, pred_4)
  return(preds)
}

predict_rf = function(rf_mod, set, test_data_predictors, 
                      inta_predictors, 
                      intap_predictors, ints_predictors, 
                      inttau_predictors=inttau_num){
  
  #pred_0 = rf_mod$test$votes[, 2]
  pred_0 = predict(rf_mod, test_data_predictors[, set], type='prob', predict.all=FALSE, norm.votes=TRUE)[, 2]
  pred_1 = predict(rf_mod, inta_predictors[, set], type='prob', predict.all=FALSE, norm.votes=TRUE)[, 2]
  pred_2 = predict(rf_mod, intap_predictors[, set], type='prob', predict.all=FALSE, norm.votes=TRUE)[, 2]
  pred_3 = predict(rf_mod, ints_predictors[, set], type='prob', predict.all=FALSE, norm.votes=TRUE)[, 2]
  pred_4 = predict(rf_mod, inttau_predictors[, set], type='prob', predict.all=FALSE, norm.votes=TRUE)[, 2]
  
  preds = data.frame(pred_0, pred_1, pred_2, pred_3, pred_4)
  return(preds)
}

predict_gbm = function(gbm_mod, set, 
                       train_data_predictors=dev_sim_t,
                       test_data_predictors=test_sim_t, 
                       inta_predictors=synth_df_inta, 
                       intap_predictors=synth_df_intap, ints_predictors=synth_df_ints, 
                       inttau_predictors=synth_df_inttau){
  
  pred_0 = predict.gbm(gbm_mod, test_data_predictors[, set], type='response') 
  pred_1 = predict.gbm(gbm_mod, inta_predictors[, set], type='response')
  pred_2 = predict.gbm(gbm_mod, intap_predictors[, set], type='response')
  pred_3 = predict.gbm(gbm_mod, ints_predictors[, set], type='response')
  pred_4 = predict.gbm(gbm_mod, inttau_predictors[, set], type='response')
  
  preds = data.frame(pred_0, pred_1, pred_2, pred_3, pred_4)
  return(preds)
}


#Trial to calibrate Random Forest: (Not published)
# Code by Zoe
dev_summary <- function(data, lev = NULL, model = NULL) {
  low_prob <- 0.000001
  high_prob <- 0.999999
  is_class1 <- ifelse(data$obs == lev[1], 1, 0)
  prob_class1 <- data[, lev[1]]
  prob_class1[prob_class1==0] <- low_prob
  prob_class1[prob_class1==1] <- high_prob
  c(deviance = -2*sum(is_class1*log(prob_class1) + ((1-is_class1)*log(1-prob_class1))), twoClassSummary(data, lev = lev))}


