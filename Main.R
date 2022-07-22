library(semTools)
library(lavaan)
library(doParallel)
library(plyr)
library(dplyr)
library(stringr)
library(gbm)
library(randomForest)
library(glmnet)
library(caret)
library(SpecsVerification)
library(arsenal)


#Specify paths
output_path = 'XXX'
local_path = 'XXXX'



source(paste0(local_path,'utility_func.R'))
source(paste0(local_path,'Impute.R'))

#Specify parameters for Simulating data
nCluster = 24
n_new_samples = 10000

#Fit prediction models
n_sim = 30000

##################
# 1. IMPUTE DATA #
##################

#Impute data and revalue dxL. Save means and sds of variables.
imputed_data = impute(imp_df)
imputed_df = imputed_data[[1]] #transformed values
means_imputed_df = imputed_data[[2]] #means
sd_imputed_df = imputed_data[[3]]  #sd

#Save means of age (needed later for simulation)
mean_a = means_imputed_df[['age']]
sd_a = sd_imputed_df[['age']]


#Simulate data
synth_df_lin = data.frame(a = imputed_df$age, ap = factor(imputed_df$apoe4), 
                          s = factor(imputed_df$sex), dx = factor(imputed_df$dx_bl, ordered=TRUE), 
                          educ = imputed_df$educ, bmi = imputed_df$bmi, 
                          hyp = factor(imputed_df$hyp_ind), 
                          alc = factor(imputed_df$alc_cat), smo = factor(imputed_df$smok_cat),
                          hi = imputed_df$hippocampus, cardio = imputed_df$cardio,
                          ab = imputed_df$abeta, tau = imputed_df$tau, 
                          ve = imputed_df$ventricles, icv = imputed_df$icv, 
                          fdg = imputed_df$fdg,
                          me = imputed_df$mmse)


model=' 
  educ ~ s + a 
  bmi ~ a + s + educ
  alc ~ a + s + educ
  smo ~ a + s + educ
  hyp ~ a + s + bmi + alc + smo + educ
  cardio ~ a + s + bmi + educ + alc + smo + hyp 
  ab ~ a + ap 
  tau ~ a + ap + alc + hyp
  dx ~ a + ap + s + educ + ab + tau + cardio + bmi
  fdg ~ a + s + dx + ab + tau + alc + smo + educ + hyp + cardio
  hi ~ a + ap + s + bmi + dx + ab + alc + tau
  ve ~ dx + ab + alc
  icv ~ a + s + dx + ab + tau
  me ~ dx + educ + a + ab + tau + fdg + ve + hi + icv
'

#Specify variable lists
numeric_cols = c('bmi', 'hi', 've', 'icv', 'fdg', 'ab', 'tau', 'educ', 'me') 
outcome = 'dx'
all = c('a', 'ap', 's', 'educ', 'bmi', 'hyp', 'alc', 'smo', 'cardio',
        'tau', 'ab', 'hi', 've', 'icv', 'fdg', 'me') 
parents = c('a', 'ap', 's','educ', 'ab', 'tau', 'cardio', 'bmi')
children = c('me', 'hi', 've', 'icv', 'fdg') 

##############
# 2. FIT SEM #
##############
fit = sem(model, data=synth_df_lin)  
genData_df = extract_semparams(parTable(fit))[[1]]
ths_df = extract_semparams(parTable(fit))[[2]]

#Create copy for tau intervention
genData_df_tau = genData_df

#Write outputfiles
write.csv(summary(fit), file = paste0(output_path, 'sem_fit', Sys.Date(), '.csv'))
write.table(genData_df, file=paste0(output_path, 'Datageneration_equations', Sys.Date(),'.csv'), sep=',')
write.table(ths_df, file=paste0(output_path, 'Datageneration_thresholds', Sys.Date(),'.csv'), sep=',')


###############################################
# 3. SIMULATE Data, train and validate models #
###############################################
#Initialize report and output files
report = init_report()

print('Running simulations and models....')
cl = parallel::makeCluster(nCluster)
registerDoParallel(cl)

report<- foreach(sim=1:n_sim, .combine=rbind)%dopar%{
  library(gbm)
  library(randomForest)
  library(DescTools)
  library(caret)
  library(SpecsVerification)
  library(dplyr)
  
  #Simulate data for training
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('a', 'ap', 's')]
  s_a = exo_boot$a
  s_ap = as.numeric(as.character(exo_boot$ap))
  s_s = as.numeric(as.character(exo_boot$s))
  dev_sim_t = generate_endogenous(genData_df, ths_df = ths_df, 
                                           a = s_a, ap = s_ap, s = s_s,  
                                           n_new_samples, numeric_cols = numeric_cols)
  
  #Simulate internal validation data
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('a', 'ap', 's')]
  s_a = exo_boot$a
  s_ap = as.numeric(as.character(exo_boot$ap))
  s_s = as.numeric(as.character(exo_boot$s))
  
  test_sim_t = generate_endogenous(genData_df, ths_df = ths_df, 
                                            a = s_a, ap = s_ap, s = s_s, 
                                            n_new_samples, numeric_cols = numeric_cols)
  
  #Simulate age-intervention data
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('ap', 's')]
  s_ap = as.numeric(as.character(exo_boot$ap))
  s_s = as.numeric(as.character(exo_boot$s))
  s_a = scale(rnorm(n_new_samples, mean=35, sd=10), center=mean_a, scale = sd_a)  
  synth_df_inta = generate_endogenous(genData_df, ths_df = ths_df,  
                                               a = s_a, ap = s_ap, s = s_s, 
                                               n_new_samples, numeric_cols = numeric_cols)
  
  #Simulate apoe-intervention data
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('a', 's')]
  s_a = exo_boot$a
  s_s = as.numeric(as.character(exo_boot$s))
  s_ap = rbinom(n_new_samples, size=1, prob=0.1)
  synth_df_intap = generate_endogenous(genData_df, ths_df = ths_df,  
                                                a = s_a, ap = s_ap, s = s_s, 
                                                n_new_samples, numeric_cols = numeric_cols)
  
  #Simulate sex-intervention data (not published)
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('a', 'ap')]
  s_a = exo_boot$a
  s_ap = as.numeric(as.character(exo_boot$ap))
  s_s = rbinom(n_new_samples, size=1, prob=0.1)
  synth_df_ints = generate_endogenous(genData_df, ths_df = ths_df, 
                                               a = s_a, ap = s_ap, s = s_s, 
                                               n_new_samples, numeric_cols = numeric_cols)
  
  #Simulate tau-intervention
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('a', 'ap', 's')]
  s_a = exo_boot$a
  s_ap = as.numeric(as.character(exo_boot$ap))
  s_s = as.numeric(as.character(exo_boot$s))
  
  genData_df_tau$model_strs[genData_df_tau$endov=='tau'] = '0.7 * a + 0.8 * ap + 0.001 * alc + 0.3 * hyp ' #+ ( -0.5 )
  
  synth_df_inttau = generate_endogenous(genData_df_tau, ths_df = ths_df, 
                                                 a = s_a, ap = s_ap, s = s_s, 
                                                 n_new_samples, numeric_cols = numeric_cols)
 
  #Copy of datasets, to convert factors to numeric --> Logistic and GBM can only run with numeric.  
  dev_num  = data.frame(sapply(dev_sim_t, as.character))
  dev_num  = data.frame(sapply(dev_num, as.numeric))
  test_num  = data.frame(sapply(test_sim_t, as.character))
  test_num  = data.frame(sapply(test_num, as.numeric))
  
  inta_num = data.frame(sapply(synth_df_inta, as.character))
  inta_num = data.frame(sapply(inta_num, as.numeric))
  
  intap_num = data.frame(sapply(synth_df_intap, as.character))
  intap_num = data.frame(sapply(intap_num, as.numeric))
  
  ints_num = data.frame(sapply(synth_df_ints, as.character))
  ints_num = data.frame(sapply(ints_num, as.numeric))
  
  inttau_num = data.frame(sapply(synth_df_inttau, as.character))
  inttau_num = data.frame(sapply(inttau_num, as.numeric))
  
  #Save summary statistics as string in report
  #For Numeric: Save mean of each variable. 
  report$sum_dev_num = paste0(round(sapply(select(dev_sim_t, all_of(c('a', numeric_cols))), mean), digits=3), collapse = ",")
  report$sum_dev_cat = paste0(round(sapply(select(dev_num, !any_of(c('a', numeric_cols))), sum)/n_new_samples*100, digits=2), collapse = ",")
  report$sum_test_num = paste0(round(sapply(select(test_sim_t, all_of(c('a', numeric_cols))), mean), digits=3), collapse = ",")
  report$sum_test_cat = paste0(round(sapply(select(test_num, !any_of(c('a', numeric_cols))), sum)/n_new_samples*100, digits=2), collapse = ",")
  report$sum_inta_num = paste0(round(sapply(select(synth_df_inta, all_of(c('a', numeric_cols))), mean), digits=3), collapse = ",")
  report$sum_inta_cat = paste0(round(sapply(select(inta_num, !any_of(c('a', numeric_cols))), sum)/n_new_samples*100, digits=2), collapse = ",")
  report$sum_intap_num = paste0(round(sapply(select(synth_df_intap, all_of(c('a', numeric_cols))), mean), digits=3), collapse = ",")
  report$sum_intap_cat = paste0(round(sapply(select(intap_num, !any_of(c('a', numeric_cols))), sum)/n_new_samples*100, digits=2), collapse = ",")
  report$sum_ints_num = paste0(round(sapply(select(synth_df_ints, all_of(c('a', numeric_cols))), mean), digits=3), collapse = ",")
  report$sum_ints_cat = paste0(round(sapply(select(ints_num, !any_of(c('a', numeric_cols))), sum)/n_new_samples*100, digits=2), collapse = ",")
  report$sum_inttau_num = paste0(round(sapply(select(synth_df_inttau, all_of(c('a', numeric_cols))), mean), digits=3), collapse = ",")
  report$sum_inttau_cat = paste0(round(sapply(select(inttau_num, !any_of(c('a', numeric_cols))), sum)/n_new_samples*100, digits=2), collapse = ",")
  
  
  ########
  #### MODELS
  ########
  # Models
  report$sim = sim
  
  # Formulas
  formula_all <- as.formula('dx ~ a + ap + s + educ + bmi + hyp + alc + smo + tau + ab + cardio + hi + ve + icv + fdg + me')
  formula_par <-as.formula('dx ~ a + ap + s + educ + ab + tau + cardio + bmi')
  formula_child <- as.formula('dx ~ me + hi + ve + icv + fdg')
  
  ###
  # Logistic Regression
  ###
  logistic_all <- glm(formula_all, data= dev_num, family = binomial(link = "logit"))  
  logistic_par <- glm(formula_par, data= dev_num, family = binomial(link = "logit")) 
  logistic_child <- glm(formula_child, data= dev_num, family = binomial(link = "logit"))
  
  #Predicting with logistic regression
  pred_logistic_all = predict_logistic(logistic_all, test_sim_t=test_num,
                                       synth_df_inta=inta_num, synth_df_intap=intap_num, 
                                       synth_df_ints=ints_num, synth_df_inttau=inttau_num)
  pred_logistic_par = predict_logistic(logistic_par, test_sim_t=test_num, synth_df_inta=inta_num, 
                                       synth_df_intap=intap_num, synth_df_ints=ints_num, 
                                       synth_df_inttau=inttau_num)
  pred_logistic_child = predict_logistic(logistic_child, test_sim_t=test_num, 
                                         synth_df_inta=inta_num, synth_df_intap=intap_num, 
                                         synth_df_ints=ints_num, synth_df_inttau=inttau_num)
  
  #Evaluation
  report$ICI0_logistic_all <- ICI(outcome = test_num[, outcome], prediction = pred_logistic_all$pred_0)
  report$ICI1_logistic_all <- ICI(outcome = inta_num[, outcome], prediction = pred_logistic_all$pred_1)
  report$ICI2_logistic_all <- ICI(outcome = intap_num[, outcome], prediction = pred_logistic_all$pred_2)
  report$ICI3_logistic_all <- ICI(outcome = ints_num[, outcome], prediction = pred_logistic_all$pred_3)
  report$ICI4_logistic_all <- ICI(outcome = inttau_num[, outcome], prediction = pred_logistic_all$pred_4)
  report$Brier0_logistic_all = brierscore(test_num[, outcome], pred_logistic_all$pred_0)
  report$Brier1_logistic_all = brierscore(inta_num[, outcome], pred_logistic_all$pred_1)
  report$Brier2_logistic_all = brierscore(intap_num[, outcome], pred_logistic_all$pred_2)
  report$Brier3_logistic_all = brierscore(ints_num[, outcome], pred_logistic_all$pred_3)
  report$Brier4_logistic_all = brierscore(inttau_num[, outcome], pred_logistic_all$pred_4)
  
  report$ICI0_logistic_par <-ICI(outcome = test_num[, outcome], prediction = pred_logistic_par$pred_0)
  report$ICI1_logistic_par <-ICI(outcome = inta_num[, outcome], prediction = pred_logistic_par$pred_1)
  report$ICI2_logistic_par <-ICI(outcome = intap_num[, outcome], prediction = pred_logistic_par$pred_2)
  report$ICI3_logistic_par <- ICI(outcome = ints_num[, outcome], prediction = pred_logistic_par$pred_3)
  report$ICI4_logistic_par <- ICI(outcome = inttau_num[, outcome], prediction = pred_logistic_par$pred_4)
  report$Brier0_logistic_par = brierscore(test_num[, outcome], pred=pred_logistic_par$pred_0)
  report$Brier1_logistic_par = brierscore(inta_num[, outcome], pred=pred_logistic_par$pred_1)
  report$Brier2_logistic_par = brierscore(intap_num[, outcome], pred=pred_logistic_par$pred_2)
  report$Brier3_logistic_par = brierscore(ints_num[, outcome], pred=pred_logistic_par$pred_3)
  report$Brier4_logistic_par = brierscore(inttau_num[, outcome], pred_logistic_par$pred_4)
  
  report$ICI0_logistic_child <-ICI(outcome = test_num[, outcome], prediction = pred_logistic_child$pred_0)
  report$ICI1_logistic_child <-ICI(outcome = inta_num[, outcome], prediction = pred_logistic_child$pred_1)
  report$ICI2_logistic_child <-ICI(outcome = intap_num[, outcome], prediction = pred_logistic_child$pred_2)
  report$ICI3_logistic_child <- ICI(outcome = ints_num[, outcome], prediction = pred_logistic_child$pred_3)
  report$ICI4_logistic_child = ICI(outcome = inttau_num[, outcome], prediction = pred_logistic_child$pred_4)
  
  report$Brier0_logistic_child = brierscore(test_num[, outcome], pred=pred_logistic_child$pred_0)
  report$Brier1_logistic_child = brierscore(inta_num[, outcome], pred=pred_logistic_child$pred_1)
  report$Brier2_logistic_child = brierscore(intap_num[, outcome], pred=pred_logistic_child$pred_2)
  report$Brier3_logistic_child = brierscore(ints_num[, outcome], pred=pred_logistic_child$pred_3)
  report$Brier4_logistic_child = brierscore(inttau_num[, outcome], pred_logistic_child$pred_4)
  
  ###
  # Elasticnet (Lasso with alpha=1)
  ###
  #Gridsearch for lambda
  lambda.1se_all = grid_search_elastnet(dev_num)
  lambda.1se_par = grid_search_elastnet(dev_num[,c(parents, outcome)])
  lambda.1se_child = grid_search_elastnet(dev_num[,c(children, outcome)])
  
  #Models
  elast_all <- glmnet::glmnet(x = as.matrix(dev_num[, all]), y = dev_num[, outcome],
                              alpha = 1, lambda = lambda.1se_all, family = "binomial") 
  elast_par <- glmnet::glmnet(x = as.matrix(dev_num[,parents]), y = dev_num[, outcome],
                              alpha = 1, lambda = lambda.1se_par, family = "binomial") 
  elast_child <- glmnet::glmnet(x = as.matrix(dev_num[,children]), y = dev_num[, outcome],
                                alpha = 1, lambda = lambda.1se_child, family = "binomial") 
  #Predictions
  pred_elastic_all = predict_elastic(elast_all, set=all, test_sim_t=test_num,
                                     synth_df_inta=inta_num, synth_df_intap=intap_num, synth_df_ints=ints_num,
                                     synth_df_inttau = inttau_num)
  pred_elastic_par = predict_elastic(elast_par, set=parents, test_sim_t=test_num, synth_df_inta=inta_num, 
                                     synth_df_intap=intap_num, synth_df_ints=ints_num,
                                     synth_df_inttau = inttau_num)
  pred_elastic_child = predict_elastic(elast_child, set=children, test_sim_t=test_num, 
                                       synth_df_inta=inta_num, synth_df_intap=intap_num, 
                                       synth_df_ints=ints_num, synth_df_inttau = inttau_num)
  
  report$no_lassopred_all <- elast_all$df
  report$no_lassopred_par <- elast_par$df
  report$no_lassopred_child <- elast_child$df
  report$lambda_all <- lambda.1se_all
  report$lambda_par <- lambda.1se_par
  report$lambda_child <-lambda.1se_child
  
  report$elast_all_a = elast_all$beta['a',]
  report$elast_all_ap = elast_all$beta['ap',]
  report$elast_all_s = elast_all$beta['s',]
  report$elast_all_educ = elast_all$beta['educ',]
  report$elast_all_bmi = elast_all$beta['bmi',]
  report$elast_all_hyp = elast_all$beta['hyp',]
  report$elast_all_alc = elast_all$beta['alc',]
  report$elast_all_smo = elast_all$beta['smo',]
  report$elast_all_cardio = elast_all$beta['cardio',]
  report$elast_all_tau = elast_all$beta['tau',]
  report$elast_all_ab = elast_all$beta['ab',]
  report$elast_all_hi = elast_all$beta['hi',]
  report$elast_all_ve = elast_all$beta['ve',]
  report$elast_all_icv = elast_all$beta['icv',]
  report$elast_all_fdg = elast_all$beta['fdg',]
  report$elast_all_me = elast_all$beta['me',]
  
  report$ICI0_elastnet_all <-ICI(outcome = test_num[, outcome], prediction = pred_elastic_all$pred_0)
  report$ICI1_elastnet_all <-ICI(outcome = inta_num[, outcome], prediction = pred_elastic_all$pred_1)
  report$ICI2_elastnet_all <-ICI(outcome = intap_num[, outcome], prediction = pred_elastic_all$pred_2)
  report$ICI3_elastnet_all <-ICI(outcome = ints_num[, outcome], prediction = pred_elastic_all$pred_3)
  report$ICI4_elastnet_all <-ICI(outcome = inttau_num[, outcome], prediction = pred_elastic_all$pred_4)
  report$Brier0_elastnet_all <-brierscore(out = test_num[, outcome], pred = pred_elastic_all$pred_0)
  report$Brier1_elastnet_all <-brierscore(out = inta_num[, outcome], pred = pred_elastic_all$pred_1)
  report$Brier2_elastnet_all <-brierscore(out = intap_num[, outcome], pred = pred_elastic_all$pred_2)
  report$Brier3_elastnet_all <-brierscore(out = ints_num[, outcome], pred = pred_elastic_all$pred_3)
  report$Brier4_elastnet_all <-brierscore(out = inttau_num[, outcome], pred = pred_elastic_all$pred_4)
  
  
  report$ICI0_elastnet_par <-ICI(outcome = test_num[, outcome], prediction = pred_elastic_par$pred_0)
  report$ICI1_elastnet_par <-ICI(outcome = inta_num[, outcome], prediction = pred_elastic_par$pred_1)
  report$ICI2_elastnet_par <-ICI(outcome = intap_num[, outcome], prediction = pred_elastic_par$pred_2)
  report$ICI3_elastnet_par <-ICI(outcome = ints_num[, outcome], prediction = pred_elastic_par$pred_3)
  report$ICI4_elastnet_par <-ICI(outcome = inttau_num[, outcome], prediction = pred_elastic_par$pred_4)
  
  report$Brier0_elastnet_par <-brierscore(out = test_num[, outcome], pred = pred_elastic_par$pred_0)
  report$Brier1_elastnet_par <-brierscore(out = inta_num[, outcome], pred = pred_elastic_par$pred_1)
  report$Brier2_elastnet_par <-brierscore(out = intap_num[, outcome], pred = pred_elastic_par$pred_2)
  report$Brier3_elastnet_par <-brierscore(out = ints_num[, outcome], pred = pred_elastic_par$pred_3)
  report$Brier4_elastnet_par <-brierscore(out = inttau_num[, outcome], pred = pred_elastic_par$pred_4)
  
  
  report$ICI0_elastnet_child <-ICI(outcome = test_num[, outcome], prediction = pred_elastic_child$pred_0)
  report$ICI1_elastnet_child <-ICI(outcome = inta_num[, outcome], prediction = pred_elastic_child$pred_1)
  report$ICI2_elastnet_child <-ICI(outcome = intap_num[, outcome], prediction = pred_elastic_child$pred_2)
  report$ICI3_elastnet_child <-ICI(outcome = ints_num[, outcome], prediction = pred_elastic_child$pred_3)
  report$ICI4_elastnet_child <-ICI(outcome = inttau_num[, outcome], prediction = pred_elastic_child$pred_4)
  
  report$Brier0_elastnet_child <-brierscore(out = test_num[, outcome], pred = pred_elastic_child$pred_0)
  report$Brier1_elastnet_child <-brierscore(out = inta_num[, outcome], pred = pred_elastic_child$pred_1)
  report$Brier2_elastnet_child <-brierscore(out = intap_num[, outcome], pred = pred_elastic_child$pred_2)
  report$Brier3_elastnet_child <-brierscore(out = ints_num[, outcome], pred = pred_elastic_child$pred_3)
  report$Brier4_elastnet_child <-brierscore(out = inttau_num[, outcome], pred = pred_elastic_child$pred_4)
  
  ###
  # Random Forest
  ###
  # Models
  rf_all = randomForest::randomForest(formula = formula_all, data = dev_sim_t, keep.forest=TRUE)
  rf_par = randomForest::randomForest(formula = formula_par, data = dev_sim_t, keep.forest=TRUE)
  rf_child = randomForest::randomForest(formula = formula_child, data = dev_sim_t, keep.forest=TRUE)
  
  #for train() to work we have to relevel the dx variable
  levels(dev_sim_t$dx) <- c("zero","one")
  
    
  #Predictions
  pred_rf_all = predict_rf(rf_all, set=all, test_data_predictors=test_sim_t, 
                           inta_predictors=synth_df_inta, 
                           intap_predictors=synth_df_intap, ints_predictors=synth_df_ints,
                           inttau_predictors=synth_df_inttau)
  pred_rf_par = predict_rf(rf_par, set=parents, test_data_predictors=test_sim_t, 
                           inta_predictors=synth_df_inta, 
                           intap_predictors=synth_df_intap, ints_predictors=synth_df_ints,
                           inttau_predictors=synth_df_inttau)
  pred_rf_child = predict_rf(rf_child, set=children, test_data_predictors=test_sim_t, 
                             inta_predictors=synth_df_inta, 
                             intap_predictors=synth_df_intap, ints_predictors=synth_df_ints,
                             inttau_predictors=synth_df_inttau)
  
  #Evaluation
  report$ICI0_rf_all <- ICI(outcome = test_num[, outcome], prediction = pred_rf_all$pred_0)
  report$ICI1_rf_all <- ICI(outcome = inta_num[, outcome], prediction = pred_rf_all$pred_1)
  report$ICI2_rf_all <- ICI(outcome = intap_num[, outcome], prediction = pred_rf_all$pred_2)
  report$ICI3_rf_all <- ICI(outcome = ints_num[, outcome], prediction = pred_rf_all$pred_3)
  report$ICI4_rf_all <- ICI(outcome = inttau_num[, outcome], prediction = pred_rf_all$pred_4)
  report$Brier0_rf_all = brierscore(test_num[, outcome], pred=pred_rf_all$pred_0)
  report$Brier1_rf_all = brierscore(inta_num[, outcome], pred=pred_rf_all$pred_1)
  report$Brier2_rf_all = brierscore(intap_num[, outcome], pred=pred_rf_all$pred_2)
  report$Brier3_rf_all = brierscore(ints_num[, outcome], pred=pred_rf_all$pred_3)
  report$Brier4_rf_all = brierscore(inttau_num[, outcome], pred=pred_rf_all$pred_4)
  
  report$ICI0_rf_par <- ICI(outcome = test_num[, outcome], prediction = pred_rf_par$pred_0)
  report$ICI1_rf_par <- ICI(outcome = inta_num[, outcome], prediction = pred_rf_par$pred_1)
  report$ICI2_rf_par <- ICI(outcome = intap_num[, outcome], prediction = pred_rf_par$pred_2)
  report$ICI3_rf_par <- ICI(outcome = ints_num[, outcome], prediction = pred_rf_par$pred_3)
  report$ICI4_rf_par <- ICI(outcome = inttau_num[, outcome], prediction = pred_rf_par$pred_4)
  
  report$Brier0_rf_par = brierscore(test_num[, outcome], pred=pred_rf_par$pred_0)
  report$Brier1_rf_par = brierscore(inta_num[, outcome], pred=pred_rf_par$pred_1)
  report$Brier2_rf_par = brierscore(intap_num[, outcome], pred=pred_rf_par$pred_2)
  report$Brier3_rf_par = brierscore(ints_num[, outcome], pred=pred_rf_par$pred_3)
  report$Brier4_rf_par = brierscore(inttau_num[, outcome], pred=pred_rf_par$pred_4)
  
  report$ICI0_rf_child <- ICI(outcome = test_num[, outcome], prediction = pred_rf_child$pred_0)
  report$ICI1_rf_child <- ICI(outcome = inta_num[, outcome], prediction = pred_rf_child$pred_1)
  report$ICI2_rf_child <- ICI(outcome = intap_num[, outcome], prediction = pred_rf_child$pred_2)
  report$ICI3_rf_child <- ICI(outcome = ints_num[, outcome], prediction = pred_rf_child$pred_3)
  report$ICI4_rf_child <- ICI(outcome = inttau_num[, outcome], prediction = pred_rf_child$pred_4)
  
  report$Brier0_rf_child = brierscore(test_num[, outcome], pred=pred_rf_child$pred_0)
  report$Brier1_rf_child = brierscore(inta_num[, outcome], pred=pred_rf_child$pred_1)
  report$Brier2_rf_child = brierscore(intap_num[, outcome], pred=pred_rf_child$pred_2)
  report$Brier3_rf_child = brierscore(ints_num[, outcome], pred=pred_rf_child$pred_3)
  report$Brier4_rf_child = brierscore(inttau_num[, outcome], pred=pred_rf_child$pred_4)
  
  ###
  # GBM
  ###
  # Models
  gbmMod_all <- gbm(formula_all, data=dev_num, distribution = "bernoulli", bag.fraction=1)
  gbmMod_par <- gbm(formula_par, data=dev_num, distribution = "bernoulli", bag.fraction=1)
  gbmMod_child <- gbm(formula_child, data=dev_num, distribution = "bernoulli", bag.fraction=1)
   
  #Predictions
  pred_gbm_all = predict_gbm(gbm_mod=gbmMod_all, set=all, 
                              train_data_predictors=dev_num,
                              test_data_predictors=test_num, 
                              inta_predictors=inta_num, 
                              intap_predictors=intap_num, ints_predictors=ints_num, 
                              inttau_predictors=inttau_num)
   
  pred_gbm_parents = predict_gbm(gbm_mod=gbmMod_par, set=parents, 
                                 train_data_predictors=dev_num,
                                 test_data_predictors=test_num, 
                                  inta_predictors=inta_num, 
                                  intap_predictors=intap_num, ints_predictors=ints_num, 
                                 inttau_predictors=inttau_num)
   
  pred_gbm_child = predict_gbm(gbm_mod=gbmMod_child, set=children, 
                                train_data_predictors=dev_num,
                                test_data_predictors=test_num, 
                                inta_predictors=inta_num, 
                                intap_predictors=intap_num, ints_predictors=ints_num, 
                               inttau_predictors=inttau_num)
  
  #Validation on internal and external data
  report$ICI0_gbm_all <-ICI(outcome = test_num[, outcome], prediction = pred_gbm_all$pred_0)
  report$ICI1_gbm_all <-ICI(outcome = inta_num[, outcome], prediction = pred_gbm_all$pred_1)
  report$ICI2_gbm_all <-ICI(outcome = intap_num[, outcome], prediction = pred_gbm_all$pred_2)
  report$ICI3_gbm_all <-ICI(outcome = ints_num[, outcome], prediction = pred_gbm_all$pred_3)
  report$ICI4_gbm_all <-ICI(outcome = inttau_num[, outcome], prediction = pred_gbm_all$pred_4)
  report$Brier0_gbm_all = brierscore(test_num[, outcome], pred=pred_gbm_all$pred_0)
  report$Brier1_gbm_all = brierscore(inta_num[, outcome], pred=pred_gbm_all$pred_1)
  report$Brier2_gbm_all = brierscore(intap_num[, outcome], pred=pred_gbm_all$pred_2)
  report$Brier3_gbm_all = brierscore(ints_num[, outcome], pred=pred_gbm_all$pred_3)
  report$Brier4_gbm_all = brierscore(inttau_num[, outcome], pred=pred_gbm_all$pred_4)
  
  report$ICI0_gbm_par <-ICI(outcome = test_num[, outcome], prediction = pred_gbm_parents$pred_0)
  report$ICI1_gbm_par <-ICI(outcome = inta_num[, outcome], prediction = pred_gbm_parents$pred_1)
  report$ICI2_gbm_par <-ICI(outcome = intap_num[, outcome], prediction = pred_gbm_parents$pred_2)
  report$ICI3_gbm_par <-ICI(outcome = ints_num[, outcome], prediction = pred_gbm_parents$pred_3)
  report$ICI4_gbm_par <-ICI(outcome = inttau_num[, outcome], prediction = pred_gbm_parents$pred_4)
  report$Brier0_gbm_par = brierscore(test_num[, outcome], pred=pred_gbm_parents$pred_0)
  report$Brier1_gbm_par = brierscore(inta_num[, outcome], pred=pred_gbm_parents$pred_1)
  report$Brier2_gbm_par = brierscore(intap_num[, outcome], pred=pred_gbm_parents$pred_2)
  report$Brier3_gbm_par = brierscore(ints_num[, outcome], pred=pred_gbm_parents$pred_3)
  report$Brier4_gbm_par = brierscore(inttau_num[, outcome], pred=pred_gbm_parents$pred_4)
   
  report$ICI0_gbm_child <-ICI(outcome = test_num[, outcome], prediction = pred_gbm_child$pred_0)
  report$ICI1_gbm_child <-ICI(outcome = inta_num[, outcome], prediction = pred_gbm_child$pred_1)
  report$ICI2_gbm_child <-ICI(outcome = intap_num[, outcome], prediction = pred_gbm_child$pred_2)
  report$ICI3_gbm_child <-ICI(outcome = ints_num[, outcome], prediction = pred_gbm_child$pred_3)
  report$ICI4_gbm_child <-ICI(outcome = inttau_num[, outcome], prediction = pred_gbm_child$pred_4)
  report$Brier0_gbm_child = brierscore(test_num[, outcome], pred=pred_gbm_child$pred_0)
  report$Brier1_gbm_child = brierscore(inta_num[, outcome], pred=pred_gbm_child$pred_1)
  report$Brier2_gbm_child = brierscore(intap_num[, outcome], pred=pred_gbm_child$pred_2)
  report$Brier3_gbm_child = brierscore(ints_num[, outcome], pred=pred_gbm_child$pred_3)
  report$Brier4_gbm_child = brierscore(inttau_num[, outcome], pred=pred_gbm_child$pred_4)
   
  return(report)
} 
stopCluster(cl)


write.csv(report, file = paste0(output_path, 'report_', Sys.Date(), Sys.time(), '.csv'))

print('Done...\n')

