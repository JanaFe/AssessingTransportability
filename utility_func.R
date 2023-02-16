#Copy from version in 'ToPublish' / Github

library(plyr)
library(dplyr)
library(mice)
library(miceadds)
library(pROC)

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
impute = function(imp_df, seed=42){
  numeric_cols = c('ventricles','hippocampus', 'icv', 'educ',
                   'age','abeta', 'tau', 'bmi', 'fdg', 'mmse') 
  factor_cols = c('dx_bl', 'sex', 'alc_cat', 'smok_cat','hyp_ind', 'apoe4', 'cardio')
  
  imp_df_sel = cbind(imp_df[, numeric_cols],
                     colwise(factor, ordered=TRUE)(imp_df[factor_cols[!factor_cols%in%c('sex', 'apoe4')]]),
                     sex=factor(imp_df$sex), apoe4 = factor(imp_df$apoe4) #Sex should not be an ordered factor.
  )
  
  impData = mice(imp_df_sel, m=1, seed=seed)
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
      param = df[row, 'est'] 
      
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
  cat_cols = c('hyp', 'alc', 'smo', 'cardio', 'dx')
  
  for(ii in genData_df$endov){
    #if(ii%in%numeric_cols){
      variance = genData_df$variances[genData_df$endov==ii]
      if(variance<0){variance=variance+1}
      std_dev = sqrt(variance)
      
      rand_noise = rnorm(n_new_samples, mean = 0, sd=std_dev)
      assign(ii, eval(parse(text=genData_df$model_strs[genData_df$endov==ii])) + rand_noise)
      
      
    #}else{ #If categorical variable
      #variance = genData_df$variances[genData_df$endov==ii]
      #std_dev = sqrt(variance)
      
      #rand_noise = rnorm(n_new_samples, mean = 0, sd=std_dev)
     # sim0 = eval(parse(text=genData_df$model_strs[genData_df$endov==ii])) + rand_noise
      #ths = ths_df$ths[which(ths_df$endov == ii)]
      #assign(ii, mapply(FUN=function(sim0, ths){ifelse(sim0>ths,1,0)}, sim0, ths))
      #assign(ii, eval(parse(text=genData_df$model_strs[genData_df$endov==ii])) + rand_noise)
    #}
  }
  
  ### APPLY THRESHOLD AFTER LOOP AT THE END:
  for(cc in cat_cols){
    sim0 = get(cc)
    ths = ths_df$ths[which(ths_df$endov == cc)]
    assign(cc, mapply(FUN=function(sim0, ths){ifelse(sim0>ths,1,0)}, sim0, ths))
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
                       ICI0_logistic_exo = NA_real_,
                       
                       ICI1_logistic_all = NA_real_,
                       ICI1_logistic_par = NA_real_,
                       ICI1_logistic_child = NA_real_,
                       ICI1_logistic_exo = NA_real_,
                       
                       ICI2_logistic_all = NA_real_,
                       ICI2_logistic_par = NA_real_,
                       ICI2_logistic_child = NA_real_,
                       ICI2_logistic_exo = NA_real_,
                       
                       ICI3_logistic_all = NA_real_,
                       ICI3_logistic_par = NA_real_,
                       ICI3_logistic_child = NA_real_,
                       ICI3_logistic_exo = NA_real_,
                       
                       ICI4_logistic_all = NA_real_,
                       ICI4_logistic_par = NA_real_,
                       ICI4_logistic_child = NA_real_,
                       ICI4_logistic_exo = NA_real_,
                       
                       ICI5_logistic_all = NA_real_,
                       ICI5_logistic_par = NA_real_,
                       ICI5_logistic_child = NA_real_,
                       ICI5_logistic_exo = NA_real_,
                       
                       ICI0_elastnet_all = NA_real_,
                       ICI0_elastnet_par = NA_real_,
                       ICI0_elastnet_child = NA_real_,
                       ICI0_elastnet_exo = NA_real_,
                       
                       ICI1_elastnet_all = NA_real_,
                       ICI1_elastnet_par = NA_real_,
                       ICI1_elastnet_child = NA_real_,
                       ICI1_elastnet_exo = NA_real_,
                       
                       ICI2_elastnet_all = NA_real_,
                       ICI2_elastnet_par = NA_real_,
                       ICI2_elastnet_child = NA_real_,
                       ICI2_elastnet_exo = NA_real_,
                       
                       ICI3_elastnet_all = NA_real_,
                       ICI3_elastnet_par = NA_real_,
                       ICI3_elastnet_child = NA_real_,
                       ICI3_elastnet_exo = NA_real_,
                       
                       ICI4_elastnet_all = NA_real_,
                       ICI4_elastnet_par = NA_real_,
                       ICI4_elastnet_child = NA_real_,
                       ICI4_elastnet_exo = NA_real_,
                       
                       ICI5_elastnet_all = NA_real_,
                       ICI5_elastnet_par = NA_real_,
                       ICI5_elastnet_child = NA_real_,
                       ICI5_elastnet_exo = NA_real_,
                       
                       ICI0_rf_all = NA_real_,
                       ICI0_rf_par = NA_real_,
                       ICI0_rf_child = NA_real_,
                       ICI0_rf_exo = NA_real_,
                       
                       ICI1_rf_all = NA_real_,
                       ICI1_rf_par = NA_real_,
                       ICI1_rf_child = NA_real_,
                       ICI1_rf_exo = NA_real_,
                       
                       ICI2_rf_all = NA_real_,
                       ICI2_rf_par = NA_real_,
                       ICI2_rf_child = NA_real_,
                       ICI2_rf_exo = NA_real_,
                       
                       ICI3_rf_all = NA_real_,
                       ICI3_rf_par = NA_real_,
                       ICI3_rf_child = NA_real_,
                       ICI3_rf_exo = NA_real_,
                       
                       ICI4_rf_all = NA_real_,
                       ICI4_rf_par = NA_real_,
                       ICI4_rf_child = NA_real_,
                       ICI4_rf_exo = NA_real_,
                       
                       ICI5_rf_all = NA_real_,
                       ICI5_rf_par = NA_real_,
                       ICI5_rf_child = NA_real_,
                       ICI5_rf_exo = NA_real_,
                       
                       ICI0_gbm_all = NA_real_,
                       ICI0_gbm_par = NA_real_,
                       ICI0_gbm_child = NA_real_,
                       ICI0_gbm_exo = NA_real_,
                       
                       ICI1_gbm_all = NA_real_,
                       ICI1_gbm_par = NA_real_,
                       ICI1_gbm_child = NA_real_,
                       ICI1_gbm_exo = NA_real_,
                       
                       ICI2_gbm_all = NA_real_,
                       ICI2_gbm_par = NA_real_,
                       ICI2_gbm_child = NA_real_,
                       ICI2_gbm_exo = NA_real_,
                       
                       ICI3_gbm_all = NA_real_,
                       ICI3_gbm_par = NA_real_,
                       ICI3_gbm_child = NA_real_,
                       ICI3_gbm_exo = NA_real_,
                       
                       ICI4_gbm_all = NA_real_,
                       ICI4_gbm_par = NA_real_,
                       ICI4_gbm_child = NA_real_,
                       ICI4_gbm_exo = NA_real_,
                       
                       ICI5_gbm_all = NA_real_,
                       ICI5_gbm_par = NA_real_,
                       ICI5_gbm_child = NA_real_,
                       ICI5_gbm_exo = NA_real_,
                       
                       Brier0_logistic_all = NA_real_,
                       Brier0_logistic_par = NA_real_,
                       Brier0_logistic_child = NA_real_,
                       Brier0_logistic_exo = NA_real_,
                       
                       Brier1_logistic_all = NA_real_,
                       Brier1_logistic_par = NA_real_,
                       Brier1_logistic_child = NA_real_,
                       Brier1_logistic_exo = NA_real_,
                       
                       Brier2_logistic_all = NA_real_,
                       Brier2_logistic_par = NA_real_,
                       Brier2_logistic_child = NA_real_,
                       Brier2_logistic_exo = NA_real_,
                       
                       Brier3_logistic_all = NA_real_,
                       Brier3_logistic_par = NA_real_,
                       Brier3_logistic_child = NA_real_,
                       Brier3_logistic_exo = NA_real_,
                       
                       Brier4_logistic_all = NA_real_,
                       Brier4_logistic_par = NA_real_,
                       Brier4_logistic_child = NA_real_,
                       Brier4_logistic_exo = NA_real_,
                       
                       Brier5_logistic_all = NA_real_,
                       Brier5_logistic_par = NA_real_,
                       Brier5_logistic_child = NA_real_,
                       Brier5_logistic_exo = NA_real_,
                       
                       Brier0_elastnet_all = NA_real_,
                       Brier0_elastnet_par = NA_real_,
                       Brier0_elastnet_child = NA_real_,
                       Brier0_elastnet_exo = NA_real_,
                       
                       Brier1_elastnet_all = NA_real_,
                       Brier1_elastnet_par = NA_real_,
                       Brier1_elastnet_child = NA_real_,
                       Brier1_elastnet_exo = NA_real_,
                       
                       Brier2_elastnet_all = NA_real_,
                       Brier2_elastnet_par = NA_real_,
                       Brier2_elastnet_child = NA_real_,
                       Brier2_elastnet_exo = NA_real_,
                       
                       Brier3_elastnet_all = NA_real_,
                       Brier3_elastnet_par = NA_real_,
                       Brier3_elastnet_child = NA_real_,
                       Brier3_elastnet_exo = NA_real_,
                       
                       Brier4_elastnet_all = NA_real_,
                       Brier4_elastnet_par = NA_real_,
                       Brier4_elastnet_child = NA_real_,
                       Brier4_elastnet_exo = NA_real_,
                       
                       Brier5_elastnet_all = NA_real_,
                       Brier5_elastnet_par = NA_real_,
                       Brier5_elastnet_child = NA_real_,
                       Brier5_elastnet_exo = NA_real_,
                       
                       Brier0_rf_all = NA_real_,
                       Brier0_rf_par = NA_real_,
                       Brier0_rf_child = NA_real_,
                       Brier0_rf_exo = NA_real_,
                       
                       Brier1_rf_all = NA_real_,
                       Brier1_rf_par = NA_real_,
                       Brier1_rf_child = NA_real_,
                       Brier1_rf_exo = NA_real_,
                       
                       Brier2_rf_all = NA_real_,
                       Brier2_rf_par = NA_real_,
                       Brier2_rf_child = NA_real_,
                       Brier2_rf_exo = NA_real_,
                       
                       Brier3_rf_all = NA_real_,
                       Brier3_rf_par = NA_real_,
                       Brier3_rf_child = NA_real_,
                       Brier3_rf_exo = NA_real_,
                       
                       Brier4_rf_all = NA_real_,
                       Brier4_rf_par = NA_real_,
                       Brier4_rf_child = NA_real_,
                       Brier4_rf_exo = NA_real_,
                       
                       Brier5_rf_all = NA_real_,
                       Brier5_rf_par = NA_real_,
                       Brier5_rf_child = NA_real_,
                       Brier5_rf_exo = NA_real_,
                       
                       Brier0_gbm_all = NA_real_,
                       Brier0_gbm_par = NA_real_,
                       Brier0_gbm_child = NA_real_,
                       Brier0_gbm_exo = NA_real_,
                      
                       Brier1_gbm_all = NA_real_,
                       Brier1_gbm_par = NA_real_,
                       Brier1_gbm_child = NA_real_,
                       Brier1_gbm_exo = NA_real_,
                       
                       Brier2_gbm_all = NA_real_,
                       Brier2_gbm_par = NA_real_,
                       Brier2_gbm_child = NA_real_,
                       Brier2_gbm_exo = NA_real_,
                       
                       Brier3_gbm_all = NA_real_,
                       Brier3_gbm_par = NA_real_,
                       Brier3_gbm_child = NA_real_,
                       Brier3_gbm_exo = NA_real_,
                       
                       Brier4_gbm_all = NA_real_,
                       Brier4_gbm_par = NA_real_,
                       Brier4_gbm_child = NA_real_,
                       Brier4_gbm_exo = NA_real_,
                       
                       Brier5_gbm_all = NA_real_,
                       Brier5_gbm_par = NA_real_,
                       Brier5_gbm_child = NA_real_,
                       Brier5_gbm_exo = NA_real_,
                       
                       tn0_logistic_all = NA_real_,
                       tn0_logistic_par = NA_real_,
                       tn0_logistic_child = NA_real_,
                       tn0_logistic_exo = NA_real_,
                       tp0_logistic_all = NA_real_,
                       tp0_logistic_par = NA_real_,
                       tp0_logistic_child = NA_real_,
                       tp0_logistic_exo = NA_real_,
                       tn1_logistic_all = NA_real_,
                       tn1_logistic_par = NA_real_,
                       tn1_logistic_child = NA_real_,
                       tn1_logistic_exo = NA_real_,
                       tp1_logistic_all = NA_real_,
                       tp1_logistic_par = NA_real_,
                       tp1_logistic_child = NA_real_,
                       tp1_logistic_exo = NA_real_,
                       tn2_logistic_all = NA_real_,
                       tn2_logistic_par = NA_real_,
                       tn2_logistic_child = NA_real_,
                       tn2_logistic_exo = NA_real_,
                       tp2_logistic_all = NA_real_,
                       tp2_logistic_par = NA_real_,
                       tp2_logistic_child = NA_real_,
                       tp2_logistic_exo = NA_real_,
                       tn3_logistic_all = NA_real_,
                       tn3_logistic_par = NA_real_,
                       tn3_logistic_child = NA_real_,
                       tn3_logistic_exo = NA_real_,
                       tp3_logistic_all = NA_real_,
                       tp3_logistic_par = NA_real_,
                       tp3_logistic_child = NA_real_,
                       tp3_logistic_exo = NA_real_,
                       tn4_logistic_all = NA_real_,
                       tn4_logistic_par = NA_real_,
                       tn4_logistic_child = NA_real_,
                       tn4_logistic_exo = NA_real_,
                       tp4_logistic_all = NA_real_,
                       tp4_logistic_par = NA_real_,
                       tp4_logistic_child = NA_real_,
                       tp4_logistic_exo = NA_real_,
                       tn5_logistic_all = NA_real_,
                       tn5_logistic_par = NA_real_,
                       tn5_logistic_child = NA_real_,
                       tn5_logistic_exo = NA_real_,
                       tp5_logistic_all = NA_real_,
                       tp5_logistic_par = NA_real_,
                       tp5_logistic_child = NA_real_,
                       tp5_logistic_exo = NA_real_,
                    
                       balacc0_logistic_all = NA_real_,
                       balacc0_logistic_par = NA_real_,
                       balacc0_logistic_child = NA_real_,
                       balacc0_logistic_exo = NA_real_,
                       
                       balacc1_logistic_all = NA_real_,
                       balacc1_logistic_par = NA_real_,
                       balacc1_logistic_child = NA_real_,
                       balacc1_logistic_exo = NA_real_,
                       
                       balacc2_logistic_all = NA_real_,
                       balacc2_logistic_par = NA_real_,
                       balacc2_logistic_child = NA_real_,
                       balacc2_logistic_exo = NA_real_,
                       
                       balacc3_logistic_all = NA_real_,
                       balacc3_logistic_par = NA_real_,
                       balacc3_logistic_child = NA_real_,
                       balacc3_logistic_exo = NA_real_,
                       
                       balacc4_logistic_all = NA_real_,
                       balacc4_logistic_par = NA_real_,
                       balacc4_logistic_child = NA_real_,
                       balacc4_logistic_exo = NA_real_,
                       
                       balacc5_logistic_all = NA_real_,
                       balacc5_logistic_par = NA_real_,
                       balacc5_logistic_child = NA_real_,
                       balacc5_logistic_exo = NA_real_,
                       
                       auc0_logistic_all = NA_real_,
                       auc0_logistic_par = NA_real_,
                       auc0_logistic_child = NA_real_,
                       auc0_logistic_exo = NA_real_,
                       
                       auc1_logistic_all = NA_real_,
                       auc1_logistic_par = NA_real_,
                       auc1_logistic_child = NA_real_,
                       auc1_logistic_exo = NA_real_,
                       
                       auc2_logistic_all = NA_real_,
                       auc2_logistic_par = NA_real_,
                       auc2_logistic_child = NA_real_,
                       auc2_logistic_exo = NA_real_,
                       
                       auc3_logistic_all = NA_real_,
                       auc3_logistic_par = NA_real_,
                       auc3_logistic_child = NA_real_,
                       auc3_logistic_exo = NA_real_,
                       
                       auc4_logistic_all = NA_real_,
                       auc4_logistic_par = NA_real_,
                       auc4_logistic_child = NA_real_,
                       auc4_logistic_exo = NA_real_,
                       
                       auc5_logistic_all = NA_real_,
                       auc5_logistic_par = NA_real_,
                       auc5_logistic_child = NA_real_,
                       auc5_logistic_exo = NA_real_,
                       
                       Fone0_logistic_all = NA_real_,
                       Fone0_logistic_par = NA_real_,
                       Fone0_logistic_child = NA_real_,
                       Fone0_logistic_exo = NA_real_,
                       
                       Fone1_logistic_all = NA_real_,
                       Fone1_logistic_par = NA_real_,
                       Fone1_logistic_child = NA_real_,
                       Fone1_logistic_exo = NA_real_,
                       
                       Fone2_logistic_all = NA_real_,
                       Fone2_logistic_par = NA_real_,
                       Fone2_logistic_child = NA_real_,
                       Fone2_logistic_exo = NA_real_,
                       
                       Fone3_logistic_all = NA_real_,
                       Fone3_logistic_par = NA_real_,
                       Fone3_logistic_child = NA_real_,
                       Fone3_logistic_exo = NA_real_,
                       
                       Fone4_logistic_all = NA_real_,
                       Fone4_logistic_par = NA_real_,
                       Fone4_logistic_child = NA_real_,
                       Fone4_logistic_exo = NA_real_,
                       
                       Fone5_logistic_all = NA_real_,
                       Fone5_logistic_par = NA_real_,
                       Fone5_logistic_child = NA_real_,
                       Fone5_logistic_exo = NA_real_,
                       
                       tn0_elastnet_all = NA_real_,
                       tn0_elastnet_par = NA_real_,
                       tn0_elastnet_child = NA_real_,
                       tn0_elastnet_exo = NA_real_,
                       tp0_elastnet_all = NA_real_,
                       tp0_elastnet_par = NA_real_,
                       tp0_elastnet_child = NA_real_,
                       tp0_elastnet_exo = NA_real_,
                       tn1_elastnet_all = NA_real_,
                       tn1_elastnet_par = NA_real_,
                       tn1_elastnet_child = NA_real_,
                       tn1_elastnet_exo = NA_real_,
                       tp1_elastnet_all = NA_real_,
                       tp1_elastnet_par = NA_real_,
                       tp1_elastnet_child = NA_real_,
                       tp1_elastnet_exo = NA_real_,
                       tn2_elastnet_all = NA_real_,
                       tn2_elastnet_par = NA_real_,
                       tn2_elastnet_child = NA_real_,
                       tn2_elastnet_exo = NA_real_,
                       tp2_elastnet_all = NA_real_,
                       tp2_elastnet_par = NA_real_,
                       tp2_elastnet_child = NA_real_,
                       tp2_elastnet_exo = NA_real_,
                       tn3_elastnet_all = NA_real_,
                       tn3_elastnet_par = NA_real_,
                       tn3_elastnet_child = NA_real_,
                       tn3_elastnet_exo = NA_real_,
                       tp3_elastnet_all = NA_real_,
                       tp3_elastnet_par = NA_real_,
                       tp3_elastnet_child = NA_real_,
                       tp3_elastnet_exo = NA_real_,
                       tn4_elastnet_all = NA_real_,
                       tn4_elastnet_par = NA_real_,
                       tn4_elastnet_child = NA_real_,
                       tn4_elastnet_exo = NA_real_,
                       tp4_elastnet_all = NA_real_,
                       tp4_elastnet_par = NA_real_,
                       tp4_elastnet_child = NA_real_,
                       tp4_elastnet_exo = NA_real_,
                       tn5_elastnet_all = NA_real_,
                       tn5_elastnet_par = NA_real_,
                       tn5_elastnet_child = NA_real_,
                       tn5_elastnet_exo = NA_real_,
                       tp5_elastnet_all = NA_real_,
                       tp5_elastnet_par = NA_real_,
                       tp5_elastnet_child = NA_real_,
                       tp5_elastnet_exo = NA_real_,
                       
                       
                       balacc0_elastnet_all = NA_real_,
                       balacc0_elastnet_par = NA_real_,
                       balacc0_elastnet_child = NA_real_,
                       balacc0_elastnet_exo = NA_real_,
                       
                       balacc1_elastnet_all = NA_real_,
                       balacc1_elastnet_par = NA_real_,
                       balacc1_elastnet_child = NA_real_,
                       balacc1_elastnet_exo = NA_real_,
                       
                       balacc2_elastnet_all = NA_real_,
                       balacc2_elastnet_par = NA_real_,
                       balacc2_elastnet_child = NA_real_,
                       balacc2_elastnet_exo = NA_real_,
                       
                       balacc3_elastnet_all = NA_real_,
                       balacc3_elastnet_par = NA_real_,
                       balacc3_elastnet_child = NA_real_,
                       balacc3_elastnet_exo = NA_real_,
                       
                       balacc4_elastnet_all = NA_real_,
                       balacc4_elastnet_par = NA_real_,
                       balacc4_elastnet_child = NA_real_,
                       balacc4_elastnet_exo = NA_real_,
                       
                       balacc5_elastnet_all = NA_real_,
                       balacc5_elastnet_par = NA_real_,
                       balacc5_elastnet_child = NA_real_,
                       balacc5_elastnet_exo = NA_real_,
                       
                       auc0_elastnet_all = NA_real_,
                       auc0_elastnet_par = NA_real_,
                       auc0_elastnet_child = NA_real_,
                       auc0_elastnet_exo = NA_real_,
                       
                       auc1_elastnet_all = NA_real_,
                       auc1_elastnet_par = NA_real_,
                       auc1_elastnet_child = NA_real_,
                       auc1_elastnet_exo = NA_real_,
                       
                       auc2_elastnet_all = NA_real_,
                       auc2_elastnet_par = NA_real_,
                       auc2_elastnet_child = NA_real_,
                       auc2_elastnet_exo = NA_real_,
                       
                       auc3_elastnet_all = NA_real_,
                       auc3_elastnet_par = NA_real_,
                       auc3_elastnet_child = NA_real_,
                       auc3_elastnet_exo = NA_real_,
                       
                       auc4_elastnet_all = NA_real_,
                       auc4_elastnet_par = NA_real_,
                       auc4_elastnet_child = NA_real_,
                       auc4_elastnet_exo = NA_real_,
                       
                       auc5_elastnet_all = NA_real_,
                       auc5_elastnet_par = NA_real_,
                       auc5_elastnet_child = NA_real_,
                       auc5_elastnet_exo = NA_real_,
                       
                       Fone0_elastnet_all = NA_real_,
                       Fone0_elastnet_par = NA_real_,
                       Fone0_elastnet_child = NA_real_,
                       Fone0_elastnet_exo = NA_real_,
                       
                       Fone1_elastnet_all = NA_real_,
                       Fone1_elastnet_par = NA_real_,
                       Fone1_elastnet_child = NA_real_,
                       Fone1_elastnet_exo = NA_real_,
                       
                       Fone2_elastnet_all = NA_real_,
                       Fone2_elastnet_par = NA_real_,
                       Fone2_elastnet_child = NA_real_,
                       Fone2_elastnet_exo = NA_real_,
                       
                       Fone3_elastnet_all = NA_real_,
                       Fone3_elastnet_par = NA_real_,
                       Fone3_elastnet_child = NA_real_,
                       Fone3_elastnet_exo = NA_real_,
                       
                       Fone4_elastnet_all = NA_real_,
                       Fone4_elastnet_par = NA_real_,
                       Fone4_elastnet_child = NA_real_,
                       Fone4_elastnet_exo = NA_real_,
                       
                       Fone5_elastnet_all = NA_real_,
                       Fone5_elastnet_par = NA_real_,
                       Fone5_elastnet_child = NA_real_,
                       Fone5_elastnet_exo = NA_real_,
                      
                       tn0_rf_all = NA_real_,
                       tn0_rf_par = NA_real_,
                       tn0_rf_child = NA_real_,
                       tn0_rf_exo = NA_real_,
                       tp0_rf_all = NA_real_,
                       tp0_rf_par = NA_real_,
                       tp0_rf_child = NA_real_,
                       tp0_rf_exo = NA_real_,
                       tn1_rf_all = NA_real_,
                       tn1_rf_par = NA_real_,
                       tn1_rf_child = NA_real_,
                       tn1_rf_exo = NA_real_,
                       tp1_rf_all = NA_real_,
                       tp1_rf_par = NA_real_,
                       tp1_rf_child = NA_real_,
                       tp1_rf_exo = NA_real_,
                       tn2_rf_all = NA_real_,
                       tn2_rf_par = NA_real_,
                       tn2_rf_child = NA_real_,
                       tn2_rf_exo = NA_real_,
                       tp2_rf_all = NA_real_,
                       tp2_rf_par = NA_real_,
                       tp2_rf_child = NA_real_,
                       tp2_rf_exo = NA_real_,
                       tn3_rf_all = NA_real_,
                       tn3_rf_par = NA_real_,
                       tn3_rf_child = NA_real_,
                       tn3_rf_exo = NA_real_,
                       tp3_rf_all = NA_real_,
                       tp3_rf_par = NA_real_,
                       tp3_rf_child = NA_real_,
                       tp3_rf_exo = NA_real_,
                       tn4_rf_all = NA_real_,
                       tn4_rf_par = NA_real_,
                       tn4_rf_child = NA_real_,
                       tn4_rf_exo = NA_real_,
                       tp4_rf_all = NA_real_,
                       tp4_rf_par = NA_real_,
                       tp4_rf_child = NA_real_,
                       tp4_rf_exo = NA_real_,
                       tn5_rf_all = NA_real_,
                       tn5_rf_par = NA_real_,
                       tn5_rf_child = NA_real_,
                       tn5_rf_exo = NA_real_,
                       tp5_rf_all = NA_real_,
                       tp5_rf_par = NA_real_,
                       tp5_rf_child = NA_real_,
                       tp5_rf_exo = NA_real_,
                       
                          
                       balacc0_rf_all = NA_real_,
                       balacc0_rf_par = NA_real_,
                       balacc0_rf_child = NA_real_,
                       balacc0_rf_exo = NA_real_,
                       
                       balacc1_rf_all = NA_real_,
                       balacc1_rf_par = NA_real_,
                       balacc1_rf_child = NA_real_,
                       balacc1_rf_exo = NA_real_,
                       
                       balacc2_rf_all = NA_real_,
                       balacc2_rf_par = NA_real_,
                       balacc2_rf_child = NA_real_,
                       balacc2_rf_exo = NA_real_,
                       
                       balacc3_rf_all = NA_real_,
                       balacc3_rf_par = NA_real_,
                       balacc3_rf_child = NA_real_,
                       balacc3_rf_exo = NA_real_,
                       
                       balacc4_rf_all = NA_real_,
                       balacc4_rf_par = NA_real_,
                       balacc4_rf_child = NA_real_,
                       balacc4_rf_exo = NA_real_,
                       
                       balacc5_rf_all = NA_real_,
                       balacc5_rf_par = NA_real_,
                       balacc5_rf_child = NA_real_,
                       balacc5_rf_exo = NA_real_,
                       
                       auc0_rf_all = NA_real_,
                       auc0_rf_par = NA_real_,
                       auc0_rf_child = NA_real_,
                       auc0_rf_exo = NA_real_,
                       
                       auc1_rf_all = NA_real_,
                       auc1_rf_par = NA_real_,
                       auc1_rf_child = NA_real_,
                       auc1_rf_exo = NA_real_,
                       
                       auc2_rf_all = NA_real_,
                       auc2_rf_par = NA_real_,
                       auc2_rf_child = NA_real_,
                       auc2_rf_exo = NA_real_,
                       
                       auc3_rf_all = NA_real_,
                       auc3_rf_par = NA_real_,
                       auc3_rf_child = NA_real_,
                       auc3_rf_exo = NA_real_,
                       
                       auc4_rf_all = NA_real_,
                       auc4_rf_par = NA_real_,
                       auc4_rf_child = NA_real_,
                       auc4_rf_exo = NA_real_,
                       
                       auc5_rf_all = NA_real_,
                       auc5_rf_par = NA_real_,
                       auc5_rf_child = NA_real_,
                       auc5_rf_exo = NA_real_,
                        
                       Fone0_rf_all = NA_real_,
                       Fone0_rf_par = NA_real_,
                       Fone0_rf_child = NA_real_,
                       Fone0_rf_exo = NA_real_,
                       
                       Fone1_rf_all = NA_real_,
                       Fone1_rf_par = NA_real_,
                       Fone1_rf_child = NA_real_,
                       Fone1_rf_exo = NA_real_,
                       
                       Fone2_rf_all = NA_real_,
                       Fone2_rf_par = NA_real_,
                       Fone2_rf_child = NA_real_,
                       Fone2_rf_exo = NA_real_,
                       
                       Fone3_rf_all = NA_real_,
                       Fone3_rf_par = NA_real_,
                       Fone3_rf_child = NA_real_,
                       Fone3_rf_exo = NA_real_,
                       
                       Fone4_rf_all = NA_real_,
                       Fone4_rf_par = NA_real_,
                       Fone4_rf_child = NA_real_,
                       Fone4_rf_exo = NA_real_,
                       
                       Fone5_rf_all = NA_real_,
                       Fone5_rf_par = NA_real_,
                       Fone5_rf_child = NA_real_,
                       Fone5_rf_exo = NA_real_,
                       
                       tn0_gbm_all = NA_real_,
                       tn0_gbm_par = NA_real_,
                       tn0_gbm_child = NA_real_,
                       tn0_gbm_exo = NA_real_,
                       tp0_gbm_all = NA_real_,
                       tp0_gbm_par = NA_real_,
                       tp0_gbm_child = NA_real_,
                       tp0_gbm_exo = NA_real_,
                       tn1_gbm_all = NA_real_,
                       tn1_gbm_par = NA_real_,
                       tn1_gbm_child = NA_real_,
                       tn1_gbm_exo = NA_real_,
                       tp1_gbm_all = NA_real_,
                       tp1_gbm_par = NA_real_,
                       tp1_gbm_child = NA_real_,
                       tp1_gbm_exo = NA_real_,
                       tn2_gbm_all = NA_real_,
                       tn2_gbm_par = NA_real_,
                       tn2_gbm_child = NA_real_,
                       tn2_gbm_exo = NA_real_,
                       tp2_gbm_all = NA_real_,
                       tp2_gbm_par = NA_real_,
                       tp2_gbm_child = NA_real_,
                       tp2_gbm_exo = NA_real_,
                       tn3_gbm_all = NA_real_,
                       tn3_gbm_par = NA_real_,
                       tn3_gbm_child = NA_real_,
                       tn3_gbm_exo = NA_real_,
                       tp3_gbm_all = NA_real_,
                       tp3_gbm_par = NA_real_,
                       tp3_gbm_child = NA_real_,
                       tp3_gbm_exo = NA_real_,
                       tn4_gbm_all = NA_real_,
                       tn4_gbm_par = NA_real_,
                       tn4_gbm_child = NA_real_,
                       tn4_gbm_exo = NA_real_,
                       tp4_gbm_all = NA_real_,
                       tp4_gbm_par = NA_real_,
                       tp4_gbm_child = NA_real_,
                       tp4_gbm_exo = NA_real_,
                       tn5_gbm_all = NA_real_,
                       tn5_gbm_par = NA_real_,
                       tn5_gbm_child = NA_real_,
                       tn5_gbm_exo = NA_real_,
                       tp5_gbm_all = NA_real_,
                       tp5_gbm_par = NA_real_,
                       tp5_gbm_child = NA_real_,
                       tp5_gbm_exo = NA_real_,
                       
                       balacc0_gbm_all = NA_real_,
                       balacc0_gbm_par = NA_real_,
                       balacc0_gbm_child = NA_real_,
                       balacc0_gbm_exo = NA_real_,
                       
                       balacc1_gbm_all = NA_real_,
                       balacc1_gbm_par = NA_real_,
                       balacc1_gbm_child = NA_real_,
                       balacc1_gbm_exo = NA_real_,
                       
                       balacc2_gbm_all = NA_real_,
                       balacc2_gbm_par = NA_real_,
                       balacc2_gbm_child = NA_real_,
                       balacc2_gbm_exo = NA_real_,
                       
                       balacc3_gbm_all = NA_real_,
                       balacc3_gbm_par = NA_real_,
                       balacc3_gbm_child = NA_real_,
                       balacc3_gbm_exo = NA_real_,
                       
                       balacc4_gbm_all = NA_real_,
                       balacc4_gbm_par = NA_real_,
                       balacc4_gbm_child = NA_real_,
                       balacc4_gbm_exo = NA_real_,
                       
                       balacc5_gbm_all = NA_real_,
                       balacc5_gbm_par = NA_real_,
                       balacc5_gbm_child = NA_real_,
                       balacc5_gbm_exo = NA_real_,
                       
                       auc0_gbm_all = NA_real_,
                       auc0_gbm_par = NA_real_,
                       auc0_gbm_child = NA_real_,
                       auc0_gbm_exo = NA_real_,
                       
                       auc1_gbm_all = NA_real_,
                       auc1_gbm_par = NA_real_,
                       auc1_gbm_child = NA_real_,
                       auc1_gbm_exo = NA_real_,
                       
                       auc2_gbm_all = NA_real_,
                       auc2_gbm_par = NA_real_,
                       auc2_gbm_child = NA_real_,
                       auc2_gbm_exo = NA_real_,
                       
                       auc3_gbm_all = NA_real_,
                       auc3_gbm_par = NA_real_,
                       auc3_gbm_child = NA_real_,
                       auc3_gbm_exo = NA_real_,
                       
                       auc4_gbm_all = NA_real_,
                       auc4_gbm_par = NA_real_,
                       auc4_gbm_child = NA_real_,
                       auc4_gbm_exo = NA_real_,
                       
                       auc5_gbm_all = NA_real_,
                       auc5_gbm_par = NA_real_,
                       auc5_gbm_child = NA_real_,
                       auc5_gbm_exo = NA_real_,
                       
                       Fone0_gbm_all = NA_real_,
                       Fone0_gbm_par = NA_real_,
                       Fone0_gbm_child = NA_real_,
                       Fone0_gbm_exo = NA_real_,
                       
                       Fone1_gbm_all = NA_real_,
                       Fone1_gbm_par = NA_real_,
                       Fone1_gbm_child = NA_real_,
                       Fone1_gbm_exo = NA_real_,
                       
                       Fone2_gbm_all = NA_real_,
                       Fone2_gbm_par = NA_real_,
                       Fone2_gbm_child = NA_real_,
                       Fone2_gbm_exo = NA_real_,
                       
                       Fone3_gbm_all = NA_real_,
                       Fone3_gbm_par = NA_real_,
                       Fone3_gbm_child = NA_real_,
                       Fone3_gbm_exo = NA_real_,
                       
                       Fone4_gbm_all = NA_real_,
                       Fone4_gbm_par = NA_real_,
                       Fone4_gbm_child = NA_real_,
                       Fone4_gbm_exo = NA_real_,
                       
                       Fone5_gbm_all = NA_real_,
                       Fone5_gbm_par = NA_real_,
                       Fone5_gbm_child = NA_real_,
                       Fone5_gbm_exo = NA_real_,
                       
                       no_lassopred_all = NA_real_,
                       no_lassopred_par = NA_real_,
                       no_lassopred_child = NA_real_,
                       no_lassopred_exo = NA_real_,
                       
                       lambda_all = NA_real_,
                       lambda_par = NA_real_,
                       lambda_child = NA_real_,
                       lambda_exo = NA_real_,
                       
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

pred_classes = function(probabilities){
  classes = ifelse(probabilities > 0.5, 1, 0)
  
  classes_fac = factor(as.character(classes), ordered=TRUE)
  
  return(classes_fac)
  
}




get_disc = function(pred_cl, true_classes){
  
  #convert true classes to factor
  true_classes = as.factor(as.character(true_classes))
  
  conf = confusionMatrix(pred_cl, true_classes, mode = "everything")
  
  acc = round(conf$overall[['Accuracy']], digits=5)
  balacc = round(conf$byClass[['Balanced Accuracy']], digits=5)
  f1 = round(conf$byClass[['F1']], digits=5)
  
  aucCI = ci.auc(response = true_classes, predictor = pred_cl)
  auc = round(aucCI[2], digits=2)
  #auc0 = round(aucCI[1], digits=2)
  #auc1 = round(aucCI[3], digits=2)
  
  tn = as.matrix(conf)[1,1]
  fp = as.matrix(conf)[2,1]
  fn = as.matrix(conf)[2,1]
  tp = as.matrix(conf)[2,2]
  
  return(data.frame(acc, balacc, auc, f1, tn, fp, fn, tp))
  
} 


# Predict functions
predict_logistic = function(logistic_mod, test_sim_t=test_sim_t, 
                            synth_df_inta=inta_num, synth_df_inta2=inta2_num,
                            synth_df_intap=intap_num, synth_df_ints=ints_num, synth_df_inttau=inttau_num){
  #Note: It is not necessary to pass the prediction 'set' here, as in rf or gbm predict-functions. The model that is passed, includes this info already
  
  pred_0 = predict(logistic_mod, type = "response", newdata = test_sim_t)
  pred_1 = predict(logistic_mod, type = "response", newdata = synth_df_inta)
  pred_2 = predict(logistic_mod, type = "response", newdata = synth_df_intap)
  pred_3 = predict(logistic_mod, type = "response", newdata = synth_df_ints)
  pred_4 = predict(logistic_mod, type = "response", newdata = synth_df_inttau)
  pred_5 = predict(logistic_mod, type = "response", newdata = synth_df_inta2)
  
  preds = data.frame(pred_0, pred_1, pred_2, pred_3, pred_4, pred_5)
  return(preds)
}




predict_elastic = function(elastic_mod, set, test_sim_t=test_sim_t, 
                            synth_df_inta=synth_df_inta, synth_df_inta2=synth_df_inta2, 
                            synth_df_intap=synth_df_intap, synth_df_ints=synth_df_ints, 
                            synth_df_inttau=inttau_num){
  s = elastic_mod$lambda
  
  pred_0 = c(predict(elastic_mod, type = "response", s=s, newx = data.matrix(test_sim_t[, set])))
  pred_1 = c(predict(elastic_mod, type = "response", s=s, newx = data.matrix(synth_df_inta[, set])))
  pred_2 = c(predict(elastic_mod, type = "response", s=s, newx = data.matrix(synth_df_intap[, set])))
  pred_3 = c(predict(elastic_mod, type = "response", s=s, newx = data.matrix(synth_df_ints[, set])))
  pred_4 = c(predict(elastic_mod, type = "response", s=s, newx = data.matrix(synth_df_inttau[, set])))
  pred_5 = c(predict(elastic_mod, type = "response", s=s, newx = data.matrix(synth_df_inta2[, set])))
  
  preds = data.frame(pred_0, pred_1, pred_2, pred_3, pred_4, pred_5)
  return(preds)
}

predict_rf = function(rf_mod, set, test_data_predictors=test_num, 
                      inta_predictors=inta_num, inta2_predictors=inta2_num, 
                      intap_predictors=intap_num, ints_predictors=ints_num, 
                      inttau_predictors=inttau_num){
  
  #pred_0 = rf_mod$test$votes[, 2]
  pred_0 = predict(rf_mod, test_data_predictors[, set], type='prob', predict.all=FALSE, norm.votes=TRUE)[, 2]
  pred_1 = predict(rf_mod, inta_predictors[, set], type='prob', predict.all=FALSE, norm.votes=TRUE)[, 2]
  pred_2 = predict(rf_mod, intap_predictors[, set], type='prob', predict.all=FALSE, norm.votes=TRUE)[, 2]
  pred_3 = predict(rf_mod, ints_predictors[, set], type='prob', predict.all=FALSE, norm.votes=TRUE)[, 2]
  pred_4 = predict(rf_mod, inttau_predictors[, set], type='prob', predict.all=FALSE, norm.votes=TRUE)[, 2]
  pred_5 = predict(rf_mod, inta2_predictors[, set], type='prob', predict.all=FALSE, norm.votes=TRUE)[, 2]
  
  preds = data.frame(pred_0, pred_1, pred_2, pred_3, pred_4, pred_5)
  return(preds)
}

predict_gbm = function(gbm_mod, set, 
                       train_data_predictors=dev_sim_t,
                       test_data_predictors=test_sim_t, 
                       inta_predictors=synth_df_inta, 
                       inta2_predictors=synth_df_inta2, 
                       intap_predictors=synth_df_intap, ints_predictors=synth_df_ints, 
                       inttau_predictors=synth_df_inttau){
  
  pred_0 = predict.gbm(gbm_mod, test_data_predictors[, set], type='response') 
  pred_1 = predict.gbm(gbm_mod, inta_predictors[, set], type='response')
  pred_2 = predict.gbm(gbm_mod, intap_predictors[, set], type='response')
  pred_3 = predict.gbm(gbm_mod, ints_predictors[, set], type='response')
  pred_4 = predict.gbm(gbm_mod, inttau_predictors[, set], type='response')
  pred_5 = predict.gbm(gbm_mod, inta2_predictors[, set], type='response')
  
  preds = data.frame(pred_0, pred_1, pred_2, pred_3, pred_4, pred_5)
  return(preds)
}


