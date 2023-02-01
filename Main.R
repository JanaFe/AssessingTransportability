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
#output_path = 'XXX'
#local_path = 'XXXX'


source(paste0(local_path,'utility_func.R'))


#Specify seed from argument. Seed should be the first argument
args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  seed=42 #default
} else  {
  seed=args[1] 
}


#Specify parameters for Simulating data
nCluster = 64
n_new_samples = 10000

#Fit prediction models
n_sim = 10000

#Load processed ADNI dataset
imp_df = read.csv(paste0(local_path,'/processed_df.csv'), header=TRUE, sep=',')


##################
# 1. IMPUTE DATA #
##################

#Impute data and revalue dxL. Save means and sds of variables.
imputed_data = impute(imp_df, seed=seed)
imputed_df = imputed_data[[1]] #transformed values
means_imputed_df = imputed_data[[2]] #means
sd_imputed_df = imputed_data[[3]]  #sd

###Save summary of imputed data for categorical variables
  #Convert to numeric
impnum = data.frame(sapply(imputed_df[,c('dx_bl', 'alc_cat', 'smok_cat', 'hyp_ind', 'cardio', 'sex', 'apoe4')], as.character))
impnum = data.frame(sapply(impnum[,c('dx_bl', 'alc_cat', 'smok_cat', 'hyp_ind', 'cardio', 'sex', 'apoe4')], as.numeric))
  # Calculate prevalence of categorical data
sumdf = apply(impnum, 2, sum)
cat_imp_sum = round((sumdf/nrow(imp_df))*100, digits=2)


#Save means of age (needed later for simulation)
mean_a = means_imputed_df[['age']]
sd_a = sd_imputed_df[['age']]


#Simulate data
synth_df_lin = data.frame(a = imputed_df$age, ap = as.numeric(as.character(imputed_df$apoe4)), 
                          s = as.numeric(as.character(imputed_df$sex)), dx = factor(imputed_df$dx_bl, ordered=TRUE), 
                          educ = imputed_df$educ, bmi = imputed_df$bmi, 
                          hyp = factor(imputed_df$hyp_ind, ordered=TRUE), 
                          alc = factor(imputed_df$alc_cat, ordered=TRUE), smo = factor(imputed_df$smok_cat, ordered=TRUE),
                          hi = imputed_df$hippocampus, cardio = factor(imputed_df$cardio, ordered=TRUE),
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
exo = c('a', 'ap', 's')

##############
# 2. FIT SEM #
##############
fit = sem(model, data=synth_df_lin)  
genData_df = extract_semparams(parTable(fit))[[1]]
ths_df = extract_semparams(parTable(fit))[[2]]

#Create copy for tau intervention
genData_df_tau = genData_df

#Write outputfiles
write.csv(summary(fit), file = paste0(output_path, 'sem_fit_', seed,'_', Sys.Date(), '.csv'))
write.table(genData_df, file=paste0(output_path, 'Datageneration_equations_', seed, '_', Sys.Date(),'.csv'), sep=',')
write.table(ths_df, file=paste0(output_path, 'Datageneration_thresholds_', seed, '_', Sys.Date(),'.csv'), sep=',')
write.table(data.frame(means_imputed_df), file=paste0(output_path, 'Means_imputed_df_', seed, '_', Sys.Date(),'.csv'), sep=',')
write.table(data.frame(sd_imputed_df), file=paste0(output_path, 'sd_imputed_df_', seed, '_', Sys.Date(),'.csv'), sep=',')
write.table(cat_imp_sum , file=paste0(output_path, 'cat_imp_sum_', seed, '_', Sys.Date(),'.csv'), sep=',')


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
  library(pROC)
  
  #Simulate data for training
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('a', 'ap', 's')]
  s_a = exo_boot$a
  s_ap = exo_boot$ap
  s_s = exo_boot$s
  dev_sim_t = generate_endogenous(genData_df, ths_df = ths_df, 
                                           a = s_a, ap = s_ap, s = s_s,  
                                           n_new_samples, numeric_cols = numeric_cols)
  
  #Simulate internal validation data
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('a', 'ap', 's')]
  s_a = exo_boot$a
  s_ap = exo_boot$ap
  s_s = exo_boot$s
  
  test_sim_t = generate_endogenous(genData_df, ths_df = ths_df, 
                                            a = s_a, ap = s_ap, s = s_s, 
                                            n_new_samples, numeric_cols = numeric_cols)
  
  #Simulate age-intervention data
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('ap', 's')]
  s_ap = exo_boot$ap
  s_s = exo_boot$s
  s_a = scale(rnorm(n_new_samples, mean=35, sd=10), center=mean_a, scale = sd_a)  
  synth_df_inta = generate_endogenous(genData_df, ths_df = ths_df,  
                                               a = s_a, ap = s_ap, s = s_s, 
                                               n_new_samples, numeric_cols = numeric_cols)
  
  #Simulate second age-intervention data
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('ap', 's')]
  s_ap = exo_boot$ap
  s_s = exo_boot$s
  s_a = scale(rnorm(n_new_samples, mean=65, sd=10), center=mean_a, scale = sd_a)  
  synth_df_inta2 = generate_endogenous(genData_df, ths_df = ths_df,  
                                      a = s_a, ap = s_ap, s = s_s, 
                                      n_new_samples, numeric_cols = numeric_cols)
  
  #Simulate apoe-intervention data
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('a', 's')]
  s_a = exo_boot$a
  s_s = as.numeric(as.character(exo_boot$s))
  s_ap = rbinom(n_new_samples, size=1, prob=0.05) #Reduce prevalence to 5%
  synth_df_intap = generate_endogenous(genData_df, ths_df = ths_df,  
                                                a = s_a, ap = s_ap, s = s_s, 
                                                n_new_samples, numeric_cols = numeric_cols)
  
  #Simulate sex-intervention data (not published)
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('a', 'ap')]
  s_a = exo_boot$a
  s_ap = exo_boot$ap
  s_s = rbinom(n_new_samples, size=1, prob=0.1)
  synth_df_ints = generate_endogenous(genData_df, ths_df = ths_df, 
                                               a = s_a, ap = s_ap, s = s_s, 
                                               n_new_samples, numeric_cols = numeric_cols)
  
  #Simulate tau-intervention
  exo_boot = synth_df_lin[sample(1:nrow(synth_df_lin), size = n_new_samples, replace=TRUE), c('a', 'ap', 's')]
  s_a = exo_boot$a
  s_ap = as.numeric(as.character(exo_boot$ap))
  s_s = as.numeric(as.character(exo_boot$s))
  
  #genData_df_tau$model_strs[genData_df_tau$endov=='tau'] = '0.7 * a + 0.8 * ap + 0.001 * alc + 0.3 * hyp ' #+ ( -0.5 )
  genData_df_tau$model_strs[genData_df_tau$endov=='tau'] = '0.9 * a + 0.8 * ap + 0.001 * alc + 0.9 * hyp +  ( 0.9 )'
  
  
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
  
  inta2_num = data.frame(sapply(synth_df_inta2, as.character))
  inta2_num = data.frame(sapply(inta2_num, as.numeric))
  
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
  report$sum_inta2_num = paste0(round(sapply(select(synth_df_inta2, all_of(c('a', numeric_cols))), mean), digits=3), collapse = ",")
  
  report$sum_inta_cat = paste0(round(sapply(select(inta_num, !any_of(c('a', numeric_cols))), sum)/n_new_samples*100, digits=2), collapse = ",")
  report$sum_inta2_cat = paste0(round(sapply(select(inta2_num, !any_of(c('a', numeric_cols))), sum)/n_new_samples*100, digits=2), collapse = ",")
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
  formula_exo <- as.formula('dx ~ a + ap + s')
  
  ###
  # Logistic Regression
  ###
  logistic_all <- glm(formula_all, data= dev_num, family = binomial(link = "logit"))  
  logistic_par <- glm(formula_par, data= dev_num, family = binomial(link = "logit")) 
  logistic_child <- glm(formula_child, data= dev_num, family = binomial(link = "logit"))
  logistic_exo <- glm(formula_exo, data= dev_num, family = binomial(link = "logit"))
  
  #Predicting with logistic regression
  pred_logistic_all = predict_logistic(logistic_all, test_sim_t=test_num,
                                       synth_df_inta=inta_num, synth_df_inta2=inta2_num, 
                                       synth_df_intap=intap_num, 
                                       synth_df_ints=ints_num, synth_df_inttau=inttau_num)
  pred_logistic_par = predict_logistic(logistic_par, test_sim_t=test_num, 
                                       synth_df_inta=inta_num, synth_df_inta2=inta2_num,
                                       synth_df_intap=intap_num, synth_df_ints=ints_num, 
                                       synth_df_inttau=inttau_num)
  pred_logistic_child = predict_logistic(logistic_child, test_sim_t=test_num, 
                                         synth_df_inta=inta_num, synth_df_inta2=inta2_num, 
                                         synth_df_intap=intap_num, 
                                         synth_df_ints=ints_num, synth_df_inttau=inttau_num)
  pred_logistic_exo = predict_logistic(logistic_exo, test_sim_t=test_num, 
                                         synth_df_inta=inta_num, synth_df_inta2=inta2_num,
                                         synth_df_intap=intap_num, 
                                         synth_df_ints=ints_num, synth_df_inttau=inttau_num)
  
  #Evaluation
  #Calibration
  report$ICI0_logistic_all <- ICI(outcome = test_num[, outcome], prediction = pred_logistic_all$pred_0)
  report$ICI1_logistic_all <- ICI(outcome = inta_num[, outcome], prediction = pred_logistic_all$pred_1)
  report$ICI2_logistic_all <- ICI(outcome = intap_num[, outcome], prediction = pred_logistic_all$pred_2)
  report$ICI3_logistic_all <- ICI(outcome = ints_num[, outcome], prediction = pred_logistic_all$pred_3)
  report$ICI4_logistic_all <- ICI(outcome = inttau_num[, outcome], prediction = pred_logistic_all$pred_4)
  report$ICI5_logistic_all <- ICI(outcome = inta2_num[, outcome], prediction = pred_logistic_all$pred_5)
  report$Brier0_logistic_all = brierscore(test_num[, outcome], pred_logistic_all$pred_0)
  report$Brier1_logistic_all = brierscore(inta_num[, outcome], pred_logistic_all$pred_1)
  report$Brier2_logistic_all = brierscore(intap_num[, outcome], pred_logistic_all$pred_2)
  report$Brier3_logistic_all = brierscore(ints_num[, outcome], pred_logistic_all$pred_3)
  report$Brier4_logistic_all = brierscore(inttau_num[, outcome], pred_logistic_all$pred_4)
  report$Brier5_logistic_all = brierscore(inta2_num[, outcome], pred_logistic_all$pred_5)
  
  report$ICI0_logistic_par <-ICI(outcome = test_num[, outcome], prediction = pred_logistic_par$pred_0)
  report$ICI1_logistic_par <-ICI(outcome = inta_num[, outcome], prediction = pred_logistic_par$pred_1)
  report$ICI2_logistic_par <-ICI(outcome = intap_num[, outcome], prediction = pred_logistic_par$pred_2)
  report$ICI3_logistic_par <- ICI(outcome = ints_num[, outcome], prediction = pred_logistic_par$pred_3)
  report$ICI4_logistic_par <- ICI(outcome = inttau_num[, outcome], prediction = pred_logistic_par$pred_4)
  report$ICI5_logistic_par <-ICI(outcome = inta2_num[, outcome], prediction = pred_logistic_par$pred_5)
  report$Brier0_logistic_par = brierscore(test_num[, outcome], pred=pred_logistic_par$pred_0)
  report$Brier1_logistic_par = brierscore(inta_num[, outcome], pred=pred_logistic_par$pred_1)
  report$Brier2_logistic_par = brierscore(intap_num[, outcome], pred=pred_logistic_par$pred_2)
  report$Brier3_logistic_par = brierscore(ints_num[, outcome], pred=pred_logistic_par$pred_3)
  report$Brier4_logistic_par = brierscore(inttau_num[, outcome], pred_logistic_par$pred_4)
  report$Brier5_logistic_par = brierscore(inta2_num[, outcome], pred_logistic_par$pred_5)
  
  report$ICI0_logistic_child <-ICI(outcome = test_num[, outcome], prediction = pred_logistic_child$pred_0)
  report$ICI1_logistic_child <-ICI(outcome = inta_num[, outcome], prediction = pred_logistic_child$pred_1)
  report$ICI2_logistic_child <-ICI(outcome = intap_num[, outcome], prediction = pred_logistic_child$pred_2)
  report$ICI3_logistic_child <- ICI(outcome = ints_num[, outcome], prediction = pred_logistic_child$pred_3)
  report$ICI4_logistic_child <- ICI(outcome = inttau_num[, outcome], prediction = pred_logistic_child$pred_4)
  report$ICI5_logistic_child <-ICI(outcome = inta2_num[, outcome], prediction = pred_logistic_child$pred_5)
  report$Brier0_logistic_child = brierscore(test_num[, outcome], pred=pred_logistic_child$pred_0)
  report$Brier1_logistic_child = brierscore(inta_num[, outcome], pred=pred_logistic_child$pred_1)
  report$Brier2_logistic_child = brierscore(intap_num[, outcome], pred=pred_logistic_child$pred_2)
  report$Brier3_logistic_child = brierscore(ints_num[, outcome], pred=pred_logistic_child$pred_3)
  report$Brier4_logistic_child = brierscore(inttau_num[, outcome], pred_logistic_child$pred_4)
  report$Brier5_logistic_child = brierscore(inta2_num[, outcome], pred=pred_logistic_child$pred_5)
  
  report$ICI0_logistic_exo <-ICI(outcome = test_num[, outcome], prediction = pred_logistic_exo$pred_0)
  report$ICI1_logistic_exo <-ICI(outcome = inta_num[, outcome], prediction = pred_logistic_exo$pred_1)
  report$ICI2_logistic_exo <-ICI(outcome = intap_num[, outcome], prediction = pred_logistic_exo$pred_2)
  report$ICI3_logistic_exo <- ICI(outcome = ints_num[, outcome], prediction = pred_logistic_exo$pred_3)
  report$ICI4_logistic_exo = ICI(outcome = inttau_num[, outcome], prediction = pred_logistic_exo$pred_4)
  report$ICI5_logistic_exo = ICI(outcome = inta2_num[, outcome], prediction = pred_logistic_exo$pred_5)
  report$Brier0_logistic_exo = brierscore(test_num[, outcome], pred=pred_logistic_exo$pred_0)
  report$Brier1_logistic_exo = brierscore(inta_num[, outcome], pred=pred_logistic_exo$pred_1)
  report$Brier2_logistic_exo = brierscore(intap_num[, outcome], pred=pred_logistic_exo$pred_2)
  report$Brier3_logistic_exo = brierscore(ints_num[, outcome], pred=pred_logistic_exo$pred_3)
  report$Brier4_logistic_exo = brierscore(inttau_num[, outcome], pred_logistic_exo$pred_4)
  report$Brier5_logistic_exo = brierscore(inta2_num[, outcome], pred=pred_logistic_exo$pred_5)
  
  
  #Discrimination metrics
  disc0_logistic_all = get_disc(pred_cl = pred_classes(pred_logistic_all$pred_0), true_classes=test_num[, outcome])
  disc1_logistic_all = get_disc(pred_cl = pred_classes(pred_logistic_all$pred_1), true_classes=inta_num[, outcome])
  disc2_logistic_all = get_disc(pred_cl = pred_classes(pred_logistic_all$pred_2), true_classes=intap_num[, outcome])
  disc3_logistic_all = get_disc(pred_cl = pred_classes(pred_logistic_all$pred_3), true_classes=ints_num[, outcome])
  disc4_logistic_all = get_disc(pred_cl = pred_classes(pred_logistic_all$pred_4), true_classes=inttau_num[, outcome])
  disc5_logistic_all = get_disc(pred_cl = pred_classes(pred_logistic_all$pred_5), true_classes=inta2_num[, outcome])
  disc0_logistic_par = get_disc(pred_cl = pred_classes(pred_logistic_par$pred_0), true_classes=test_num[, outcome])
  disc1_logistic_par = get_disc(pred_cl = pred_classes(pred_logistic_par$pred_1), true_classes=inta_num[, outcome])
  disc2_logistic_par = get_disc(pred_cl = pred_classes(pred_logistic_par$pred_2), true_classes=intap_num[, outcome])
  disc3_logistic_par = get_disc(pred_cl = pred_classes(pred_logistic_par$pred_3), true_classes=ints_num[, outcome])
  disc4_logistic_par = get_disc(pred_cl = pred_classes(pred_logistic_par$pred_4), true_classes=inttau_num[, outcome])
  disc5_logistic_par = get_disc(pred_cl = pred_classes(pred_logistic_par$pred_5), true_classes=inta2_num[, outcome])
  disc0_logistic_child = get_disc(pred_cl = pred_classes(pred_logistic_child$pred_0), true_classes=test_num[, outcome])
  disc1_logistic_child = get_disc(pred_cl = pred_classes(pred_logistic_child$pred_1), true_classes=inta_num[, outcome])
  disc2_logistic_child = get_disc(pred_cl = pred_classes(pred_logistic_child$pred_2), true_classes=intap_num[, outcome])
  disc3_logistic_child = get_disc(pred_cl = pred_classes(pred_logistic_child$pred_3), true_classes=ints_num[, outcome])
  disc4_logistic_child = get_disc(pred_cl = pred_classes(pred_logistic_child$pred_4), true_classes=inttau_num[, outcome])
  disc5_logistic_child = get_disc(pred_cl = pred_classes(pred_logistic_child$pred_5), true_classes=inta2_num[, outcome])
  disc0_logistic_exo = get_disc(pred_cl = pred_classes(pred_logistic_exo$pred_0), true_classes=test_num[, outcome])
  disc1_logistic_exo = get_disc(pred_cl = pred_classes(pred_logistic_exo$pred_1), true_classes=inta_num[, outcome])
  disc2_logistic_exo = get_disc(pred_cl = pred_classes(pred_logistic_exo$pred_2), true_classes=intap_num[, outcome])
  disc3_logistic_exo = get_disc(pred_cl = pred_classes(pred_logistic_exo$pred_3), true_classes=ints_num[, outcome])
  disc4_logistic_exo = get_disc(pred_cl = pred_classes(pred_logistic_exo$pred_4), true_classes=inttau_num[, outcome])
  disc5_logistic_exo = get_disc(pred_cl = pred_classes(pred_logistic_exo$pred_5), true_classes=inta2_num[, outcome])
  
  
  report$tn0_logistic_all = disc0_logistic_all$tn #internal val.
  report$tn1_logistic_all = disc1_logistic_all$tn 
  report$tn2_logistic_all = disc2_logistic_all$tn
  report$tn3_logistic_all = disc3_logistic_all$tn
  report$tn4_logistic_all = disc4_logistic_all$tn
  report$tn5_logistic_all = disc5_logistic_all$tn
  report$tp0_logistic_all = disc0_logistic_all$tp #internal val.
  report$tp1_logistic_all = disc1_logistic_all$tp 
  report$tp2_logistic_all = disc2_logistic_all$tp
  report$tp3_logistic_all = disc3_logistic_all$tp
  report$tp4_logistic_all = disc4_logistic_all$tp
  report$tp5_logistic_all = disc5_logistic_all$tp
  
  report$tn0_logistic_par = disc0_logistic_par$tn #internal val.
  report$tn1_logistic_par = disc1_logistic_par$tn 
  report$tn2_logistic_par = disc2_logistic_par$tn
  report$tn3_logistic_par = disc3_logistic_par$tn
  report$tn4_logistic_par = disc4_logistic_par$tn
  report$tn5_logistic_par = disc5_logistic_par$tn
  report$tp0_logistic_par = disc0_logistic_par$tp #internal val.
  report$tp1_logistic_par = disc1_logistic_par$tp 
  report$tp2_logistic_par = disc2_logistic_par$tp
  report$tp3_logistic_par = disc3_logistic_par$tp
  report$tp4_logistic_par = disc4_logistic_par$tp
  report$tp5_logistic_par = disc5_logistic_par$tp
  
  report$tn0_logistic_child = disc0_logistic_child$tn #internal val.
  report$tn1_logistic_child = disc1_logistic_child$tn 
  report$tn2_logistic_child = disc2_logistic_child$tn
  report$tn3_logistic_child = disc3_logistic_child$tn
  report$tn4_logistic_child = disc4_logistic_child$tn
  report$tn5_logistic_child = disc5_logistic_child$tn
  report$tp0_logistic_child = disc0_logistic_child$tp #internal val.
  report$tp1_logistic_child = disc1_logistic_child$tp 
  report$tp2_logistic_child = disc2_logistic_child$tp
  report$tp3_logistic_child = disc3_logistic_child$tp
  report$tp4_logistic_child = disc4_logistic_child$tp
  report$tp5_logistic_child = disc5_logistic_child$tp
  
  report$tn0_logistic_exo = disc0_logistic_exo$tn #internal val.
  report$tn1_logistic_exo = disc1_logistic_exo$tn 
  report$tn2_logistic_exo = disc2_logistic_exo$tn
  report$tn3_logistic_exo = disc3_logistic_exo$tn
  report$tn4_logistic_exo = disc4_logistic_exo$tn
  report$tn5_logistic_exo = disc5_logistic_exo$tn
  report$tp0_logistic_exo = disc0_logistic_exo$tp #internal val.
  report$tp1_logistic_exo = disc1_logistic_exo$tp 
  report$tp2_logistic_exo = disc2_logistic_exo$tp
  report$tp3_logistic_exo = disc3_logistic_exo$tp
  report$tp4_logistic_exo = disc4_logistic_exo$tp
  report$tp5_logistic_exo = disc5_logistic_exo$tp
  
  #AUC
  report$auc0_logistic_all = disc0_logistic_all$auc #internal val.
  report$auc1_logistic_all = disc1_logistic_all$auc 
  report$auc2_logistic_all = disc2_logistic_all$auc
  report$auc3_logistic_all = disc3_logistic_all$auc
  report$auc4_logistic_all = disc4_logistic_all$auc
  report$auc5_logistic_all = disc5_logistic_all$auc
  report$auc0_logistic_par = disc0_logistic_par$auc #internal val.
  report$auc1_logistic_par = disc1_logistic_par$auc 
  report$auc2_logistic_par = disc2_logistic_par$auc
  report$auc3_logistic_par = disc3_logistic_par$auc
  report$auc4_logistic_par = disc4_logistic_par$auc
  report$auc5_logistic_par = disc5_logistic_par$auc
  report$auc0_logistic_child = disc0_logistic_child$auc #internal val.
  report$auc1_logistic_child = disc1_logistic_child$auc 
  report$auc2_logistic_child = disc2_logistic_child$auc
  report$auc3_logistic_child = disc3_logistic_child$auc
  report$auc4_logistic_child = disc4_logistic_child$auc
  report$auc5_logistic_child = disc5_logistic_child$auc
  report$auc0_logistic_exo = disc0_logistic_exo$auc #internal val.
  report$auc1_logistic_exo = disc1_logistic_exo$auc 
  report$auc2_logistic_exo = disc2_logistic_exo$auc
  report$auc3_logistic_exo = disc3_logistic_exo$auc
  report$auc4_logistic_exo = disc4_logistic_exo$auc
  report$auc5_logistic_exo = disc5_logistic_exo$auc
  
  #Balanced Acc
  report$balacc0_logistic_all = disc0_logistic_all$balacc #internal val.
  report$balacc1_logistic_all = disc1_logistic_all$balacc 
  report$balacc2_logistic_all = disc2_logistic_all$balacc
  report$balacc3_logistic_all = disc3_logistic_all$balacc
  report$balacc4_logistic_all = disc4_logistic_all$balacc
  report$balacc5_logistic_all = disc5_logistic_all$balacc
  report$balacc0_logistic_par = disc0_logistic_par$balacc #internal val.
  report$balacc1_logistic_par = disc1_logistic_par$balacc 
  report$balacc2_logistic_par = disc2_logistic_par$balacc
  report$balacc3_logistic_par = disc3_logistic_par$balacc
  report$balacc4_logistic_par = disc4_logistic_par$balacc
  report$balacc5_logistic_par = disc5_logistic_par$balacc
  report$balacc0_logistic_child = disc0_logistic_child$balacc #internal val.
  report$balacc1_logistic_child = disc1_logistic_child$balacc 
  report$balacc2_logistic_child = disc2_logistic_child$balacc
  report$balacc3_logistic_child = disc3_logistic_child$balacc
  report$balacc4_logistic_child = disc4_logistic_child$balacc
  report$balacc5_logistic_child = disc5_logistic_child$balacc
  report$balacc0_logistic_exo = disc0_logistic_exo$balacc #internal val.
  report$balacc1_logistic_exo = disc1_logistic_exo$balacc 
  report$balacc2_logistic_exo = disc2_logistic_exo$balacc
  report$balacc3_logistic_exo = disc3_logistic_exo$balacc
  report$balacc4_logistic_exo = disc4_logistic_exo$balacc
  report$balacc5_logistic_exo = disc5_logistic_exo$balacc
  
  #F1
  report$Fone0_logistic_all = disc0_logistic_all$f1 #internal val.
  report$Fone1_logistic_all = disc1_logistic_all$f1 
  report$Fone2_logistic_all = disc2_logistic_all$f1
  report$Fone3_logistic_all = disc3_logistic_all$f1
  report$Fone4_logistic_all = disc4_logistic_all$f1
  report$Fone5_logistic_all = disc5_logistic_all$f1
  report$Fone0_logistic_par = disc0_logistic_par$f1 #internal val.
  report$Fone1_logistic_par = disc1_logistic_par$f1 
  report$Fone2_logistic_par = disc2_logistic_par$f1
  report$Fone3_logistic_par = disc3_logistic_par$f1
  report$Fone4_logistic_par = disc4_logistic_par$f1
  report$Fone5_logistic_par = disc5_logistic_par$f1
  report$Fone0_logistic_child = disc0_logistic_child$f1 #internal val.
  report$Fone1_logistic_child = disc1_logistic_child$f1 
  report$Fone2_logistic_child = disc2_logistic_child$f1
  report$Fone3_logistic_child = disc3_logistic_child$f1
  report$Fone4_logistic_child = disc4_logistic_child$f1
  report$Fone5_logistic_child = disc5_logistic_child$f1
  report$Fone0_logistic_exo = disc0_logistic_exo$f1 #internal val.
  report$Fone1_logistic_exo = disc1_logistic_exo$f1 
  report$Fone2_logistic_exo = disc2_logistic_exo$f1
  report$Fone3_logistic_exo = disc3_logistic_exo$f1
  report$Fone4_logistic_exo = disc4_logistic_exo$f1
  report$Fone5_logistic_exo = disc5_logistic_exo$f1
  
  
  ###
  # Elasticnet (Lasso with alpha=1)
  ###
  #Gridsearch for lambda
  lambda.1se_all = grid_search_elastnet(dev_num)
  lambda.1se_par = grid_search_elastnet(dev_num[,c(parents, outcome)])
  lambda.1se_child = grid_search_elastnet(dev_num[,c(children, outcome)])
  lambda.1se_exo = grid_search_elastnet(dev_num[,c(exo, outcome)])
  
  #Models
  elast_all <- glmnet::glmnet(x = as.matrix(dev_num[, all]), y = dev_num[, outcome],
                              alpha = 1, lambda = lambda.1se_all, family = "binomial") 
  elast_par <- glmnet::glmnet(x = as.matrix(dev_num[,parents]), y = dev_num[, outcome],
                              alpha = 1, lambda = lambda.1se_par, family = "binomial") 
  elast_child <- glmnet::glmnet(x = as.matrix(dev_num[,children]), y = dev_num[, outcome],
                                alpha = 1, lambda = lambda.1se_child, family = "binomial")
  elast_exo <- glmnet::glmnet(x = as.matrix(dev_num[,exo]), y = dev_num[, outcome],
                                alpha = 1, lambda = lambda.1se_exo, family = "binomial")
  
  #Predictions
  pred_elastic_all = predict_elastic(elast_all, set=all, test_sim_t=test_num,
                                     synth_df_inta=inta_num, synth_df_inta2=inta2_num,
                                     synth_df_intap=intap_num, synth_df_ints=ints_num,
                                     synth_df_inttau = inttau_num)
  pred_elastic_par = predict_elastic(elast_par, set=parents, test_sim_t=test_num, 
                                     synth_df_inta=inta_num, synth_df_inta2=inta2_num,
                                     synth_df_intap=intap_num, synth_df_ints=ints_num,
                                     synth_df_inttau = inttau_num)
  pred_elastic_child = predict_elastic(elast_child, set=children, test_sim_t=test_num, 
                                       synth_df_inta=inta_num, synth_df_inta2=inta2_num,
                                       synth_df_intap=intap_num, 
                                       synth_df_ints=ints_num, synth_df_inttau = inttau_num)
  pred_elastic_exo = predict_elastic(elast_exo, set=exo, test_sim_t=test_num, 
                                     synth_df_inta=inta_num, synth_df_inta2=inta2_num,
                                     synth_df_intap=intap_num, 
                                     synth_df_ints=ints_num, synth_df_inttau = inttau_num)
  
  report$no_lassopred_all <- elast_all$df
  report$no_lassopred_par <- elast_par$df
  report$no_lassopred_child <- elast_child$df
  report$no_lassopred_exo <- elast_exo$df
  report$lambda_all <- lambda.1se_all
  report$lambda_par <- lambda.1se_par
  report$lambda_child <-lambda.1se_child
  report$lambda_exo <-lambda.1se_exo
  
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
  report$ICI5_elastnet_all <-ICI(outcome = inta2_num[, outcome], prediction = pred_elastic_all$pred_5)
  report$ICI0_elastnet_par <-ICI(outcome = test_num[, outcome], prediction = pred_elastic_par$pred_0)
  report$ICI1_elastnet_par <-ICI(outcome = inta_num[, outcome], prediction = pred_elastic_par$pred_1)
  report$ICI2_elastnet_par <-ICI(outcome = intap_num[, outcome], prediction = pred_elastic_par$pred_2)
  report$ICI3_elastnet_par <-ICI(outcome = ints_num[, outcome], prediction = pred_elastic_par$pred_3)
  report$ICI4_elastnet_par <-ICI(outcome = inttau_num[, outcome], prediction = pred_elastic_par$pred_4)
  report$ICI5_elastnet_par <-ICI(outcome = inta2_num[, outcome], prediction = pred_elastic_par$pred_5)
  report$ICI0_elastnet_child <-ICI(outcome = test_num[, outcome], prediction = pred_elastic_child$pred_0)
  report$ICI1_elastnet_child <-ICI(outcome = inta_num[, outcome], prediction = pred_elastic_child$pred_1)
  report$ICI2_elastnet_child <-ICI(outcome = intap_num[, outcome], prediction = pred_elastic_child$pred_2)
  report$ICI3_elastnet_child <-ICI(outcome = ints_num[, outcome], prediction = pred_elastic_child$pred_3)
  report$ICI4_elastnet_child <-ICI(outcome = inttau_num[, outcome], prediction = pred_elastic_child$pred_4)
  report$ICI5_elastnet_child <-ICI(outcome = inta2_num[, outcome], prediction = pred_elastic_child$pred_5)
  report$ICI0_elastnet_exo <-ICI(outcome = test_num[, outcome], prediction = pred_elastic_exo$pred_0)
  report$ICI1_elastnet_exo <-ICI(outcome = inta_num[, outcome], prediction = pred_elastic_exo$pred_1)
  report$ICI2_elastnet_exo <-ICI(outcome = intap_num[, outcome], prediction = pred_elastic_exo$pred_2)
  report$ICI3_elastnet_exo <-ICI(outcome = ints_num[, outcome], prediction = pred_elastic_exo$pred_3)
  report$ICI4_elastnet_exo <-ICI(outcome = inttau_num[, outcome], prediction = pred_elastic_exo$pred_4)
  report$ICI5_elastnet_exo <-ICI(outcome = inta2_num[, outcome], prediction = pred_elastic_exo$pred_5)
  
  report$Brier0_elastnet_all <-brierscore(out = test_num[, outcome], pred = pred_elastic_all$pred_0)
  report$Brier1_elastnet_all <-brierscore(out = inta_num[, outcome], pred = pred_elastic_all$pred_1)
  report$Brier2_elastnet_all <-brierscore(out = intap_num[, outcome], pred = pred_elastic_all$pred_2)
  report$Brier3_elastnet_all <-brierscore(out = ints_num[, outcome], pred = pred_elastic_all$pred_3)
  report$Brier4_elastnet_all <-brierscore(out = inttau_num[, outcome], pred = pred_elastic_all$pred_4)
  report$Brier5_elastnet_all <-brierscore(out = inta2_num[, outcome], pred = pred_elastic_all$pred_5)
  report$Brier0_elastnet_par <-brierscore(out = test_num[, outcome], pred = pred_elastic_par$pred_0)
  report$Brier1_elastnet_par <-brierscore(out = inta_num[, outcome], pred = pred_elastic_par$pred_1)
  report$Brier2_elastnet_par <-brierscore(out = intap_num[, outcome], pred = pred_elastic_par$pred_2)
  report$Brier3_elastnet_par <-brierscore(out = ints_num[, outcome], pred = pred_elastic_par$pred_3)
  report$Brier4_elastnet_par <-brierscore(out = inttau_num[, outcome], pred = pred_elastic_par$pred_4)
  report$Brier5_elastnet_par <-brierscore(out = inta2_num[, outcome], pred = pred_elastic_par$pred_5)
  report$Brier0_elastnet_child <-brierscore(out = test_num[, outcome], pred = pred_elastic_child$pred_0)
  report$Brier1_elastnet_child <-brierscore(out = inta_num[, outcome], pred = pred_elastic_child$pred_1)
  report$Brier2_elastnet_child <-brierscore(out = intap_num[, outcome], pred = pred_elastic_child$pred_2)
  report$Brier3_elastnet_child <-brierscore(out = ints_num[, outcome], pred = pred_elastic_child$pred_3)
  report$Brier4_elastnet_child <-brierscore(out = inttau_num[, outcome], pred = pred_elastic_child$pred_4)
  report$Brier5_elastnet_child <-brierscore(out = inta2_num[, outcome], pred = pred_elastic_child$pred_5)
  report$Brier0_elastnet_exo <-brierscore(out = test_num[, outcome], pred = pred_elastic_exo$pred_0)
  report$Brier1_elastnet_exo <-brierscore(out = inta_num[, outcome], pred = pred_elastic_exo$pred_1)
  report$Brier2_elastnet_exo <-brierscore(out = intap_num[, outcome], pred = pred_elastic_exo$pred_2)
  report$Brier3_elastnet_exo <-brierscore(out = ints_num[, outcome], pred = pred_elastic_exo$pred_3)
  report$Brier4_elastnet_exo <-brierscore(out = inttau_num[, outcome], pred = pred_elastic_exo$pred_4)
  report$Brier5_elastnet_exo <-brierscore(out = inta2_num[, outcome], pred = pred_elastic_exo$pred_5)
  
  disc0_elastnet_all = get_disc(pred_cl = pred_classes(pred_elastic_all$pred_0), true_classes=test_num[, outcome])
  disc1_elastnet_all = get_disc(pred_cl = pred_classes(pred_elastic_all$pred_1), true_classes=inta_num[, outcome])
  disc2_elastnet_all = get_disc(pred_cl = pred_classes(pred_elastic_all$pred_2), true_classes=intap_num[, outcome])
  disc3_elastnet_all = get_disc(pred_cl = pred_classes(pred_elastic_all$pred_3), true_classes=ints_num[, outcome])
  disc4_elastnet_all = get_disc(pred_cl = pred_classes(pred_elastic_all$pred_4), true_classes=inttau_num[, outcome])
  disc5_elastnet_all = get_disc(pred_cl = pred_classes(pred_elastic_all$pred_5), true_classes=inta2_num[, outcome])
  disc0_elastnet_par = get_disc(pred_cl = pred_classes(pred_elastic_par$pred_0), true_classes=test_num[, outcome])
  disc1_elastnet_par = get_disc(pred_cl = pred_classes(pred_elastic_par$pred_1), true_classes=inta_num[, outcome])
  disc2_elastnet_par = get_disc(pred_cl = pred_classes(pred_elastic_par$pred_2), true_classes=intap_num[, outcome])
  disc3_elastnet_par = get_disc(pred_cl = pred_classes(pred_elastic_par$pred_3), true_classes=ints_num[, outcome])
  disc4_elastnet_par = get_disc(pred_cl = pred_classes(pred_elastic_par$pred_4), true_classes=inttau_num[, outcome])
  disc5_elastnet_par = get_disc(pred_cl = pred_classes(pred_elastic_par$pred_5), true_classes=inta2_num[, outcome])
  disc0_elastnet_child = get_disc(pred_cl = pred_classes(pred_elastic_child$pred_0), true_classes=test_num[, outcome])
  disc1_elastnet_child = get_disc(pred_cl = pred_classes(pred_elastic_child$pred_1), true_classes=inta_num[, outcome])
  disc2_elastnet_child = get_disc(pred_cl = pred_classes(pred_elastic_child$pred_2), true_classes=intap_num[, outcome])
  disc3_elastnet_child = get_disc(pred_cl = pred_classes(pred_elastic_child$pred_3), true_classes=ints_num[, outcome])
  disc4_elastnet_child = get_disc(pred_cl = pred_classes(pred_elastic_child$pred_4), true_classes=inttau_num[, outcome])
  disc5_elastnet_child = get_disc(pred_cl = pred_classes(pred_elastic_child$pred_5), true_classes=inta2_num[, outcome])
  disc0_elastnet_exo = get_disc(pred_cl = pred_classes(pred_elastic_exo$pred_0), true_classes=test_num[, outcome])
  disc1_elastnet_exo = get_disc(pred_cl = pred_classes(pred_elastic_exo$pred_1), true_classes=inta_num[, outcome])
  disc2_elastnet_exo = get_disc(pred_cl = pred_classes(pred_elastic_exo$pred_2), true_classes=intap_num[, outcome])
  disc3_elastnet_exo = get_disc(pred_cl = pred_classes(pred_elastic_exo$pred_3), true_classes=ints_num[, outcome])
  disc4_elastnet_exo = get_disc(pred_cl = pred_classes(pred_elastic_exo$pred_4), true_classes=inttau_num[, outcome])
  disc5_elastnet_exo = get_disc(pred_cl = pred_classes(pred_elastic_exo$pred_5), true_classes=inta2_num[, outcome])
  
  #AUC
  report$auc0_elastnet_all = disc0_elastnet_all$auc #internal val.
  report$auc1_elastnet_all = disc1_elastnet_all$auc 
  report$auc2_elastnet_all = disc2_elastnet_all$auc
  report$auc3_elastnet_all = disc3_elastnet_all$auc
  report$auc4_elastnet_all = disc4_elastnet_all$auc
  report$auc5_elastnet_all = disc5_elastnet_all$auc
  report$auc0_elastnet_par = disc0_elastnet_par$auc #internal val.
  report$auc1_elastnet_par = disc1_elastnet_par$auc 
  report$auc2_elastnet_par = disc2_elastnet_par$auc
  report$auc3_elastnet_par = disc3_elastnet_par$auc
  report$auc4_elastnet_par = disc4_elastnet_par$auc
  report$auc5_elastnet_par = disc5_elastnet_par$auc
  report$auc0_elastnet_child = disc0_elastnet_child$auc #internal val.
  report$auc1_elastnet_child = disc1_elastnet_child$auc 
  report$auc2_elastnet_child = disc2_elastnet_child$auc
  report$auc3_elastnet_child = disc3_elastnet_child$auc
  report$auc4_elastnet_child = disc4_elastnet_child$auc
  report$auc5_elastnet_child = disc5_elastnet_child$auc
  report$auc0_elastnet_exo = disc0_elastnet_exo$auc #internal val.
  report$auc1_elastnet_exo = disc1_elastnet_exo$auc 
  report$auc2_elastnet_exo = disc2_elastnet_exo$auc
  report$auc3_elastnet_exo = disc3_elastnet_exo$auc
  report$auc4_elastnet_exo = disc4_elastnet_exo$auc
  report$auc5_elastnet_exo = disc5_elastnet_exo$auc
  
  #Balacc & F1
  report$balacc0_elastnet_all = disc0_elastnet_all$balacc #internal val.
  report$balacc1_elastnet_all = disc1_elastnet_all$balacc 
  report$balacc2_elastnet_all = disc2_elastnet_all$balacc
  report$balacc3_elastnet_all = disc3_elastnet_all$balacc
  report$balacc4_elastnet_all = disc4_elastnet_all$balacc
  report$balacc5_elastnet_all = disc5_elastnet_all$balacc
  report$Fone0_elastnet_all = disc0_elastnet_all$f1 #internal val.
  report$Fone1_elastnet_all = disc1_elastnet_all$f1 
  report$Fone2_elastnet_all = disc2_elastnet_all$f1
  report$Fone3_elastnet_all = disc3_elastnet_all$f1
  report$Fone4_elastnet_all = disc4_elastnet_all$f1
  report$Fone5_elastnet_all = disc5_elastnet_all$f1
  
  report$balacc0_elastnet_par = disc0_elastnet_par$balacc #internal val.
  report$balacc1_elastnet_par = disc1_elastnet_par$balacc 
  report$balacc2_elastnet_par = disc2_elastnet_par$balacc
  report$balacc3_elastnet_par = disc3_elastnet_par$balacc
  report$balacc4_elastnet_par = disc4_elastnet_par$balacc
  report$balacc5_elastnet_par = disc5_elastnet_par$balacc
  report$Fone0_elastnet_par = disc0_elastnet_par$f1 #internal val.
  report$Fone1_elastnet_par = disc1_elastnet_par$f1 
  report$Fone2_elastnet_par = disc2_elastnet_par$f1
  report$Fone3_elastnet_par = disc3_elastnet_par$f1
  report$Fone4_elastnet_par = disc4_elastnet_par$f1
  report$Fone5_elastnet_par = disc5_elastnet_par$f1
  
  report$balacc0_elastnet_child = disc0_elastnet_child$balacc #internal val.
  report$balacc1_elastnet_child = disc1_elastnet_child$balacc 
  report$balacc2_elastnet_child = disc2_elastnet_child$balacc
  report$balacc3_elastnet_child = disc3_elastnet_child$balacc
  report$balacc4_elastnet_child = disc4_elastnet_child$balacc
  report$balacc5_elastnet_child = disc5_elastnet_child$balacc
  report$Fone0_elastnet_child = disc0_elastnet_child$f1 #internal val.
  report$Fone1_elastnet_child = disc1_elastnet_child$f1 
  report$Fone2_elastnet_child = disc2_elastnet_child$f1
  report$Fone3_elastnet_child = disc3_elastnet_child$f1
  report$Fone4_elastnet_child = disc4_elastnet_child$f1
  report$Fone5_elastnet_child = disc5_elastnet_child$f1
  
  report$balacc0_elastnet_exo = disc0_elastnet_exo$balacc #internal val.
  report$balacc1_elastnet_exo = disc1_elastnet_exo$balacc 
  report$balacc2_elastnet_exo = disc2_elastnet_exo$balacc
  report$balacc3_elastnet_exo = disc3_elastnet_exo$balacc
  report$balacc4_elastnet_exo = disc4_elastnet_exo$balacc
  report$balacc5_elastnet_exo = disc5_elastnet_exo$balacc
  report$Fone0_elastnet_exo = disc0_elastnet_exo$f1 #internal val.
  report$Fone1_elastnet_exo = disc1_elastnet_exo$f1 
  report$Fone2_elastnet_exo = disc2_elastnet_exo$f1
  report$Fone3_elastnet_exo = disc3_elastnet_exo$f1
  report$Fone4_elastnet_exo = disc4_elastnet_exo$f1
  report$Fone5_elastnet_exo = disc5_elastnet_exo$f1
  
  report$tn0_elastnet_all = disc0_elastnet_all$tn #internal val.
  report$tn1_elastnet_all = disc1_elastnet_all$tn 
  report$tn2_elastnet_all = disc2_elastnet_all$tn
  report$tn3_elastnet_all = disc3_elastnet_all$tn
  report$tn4_elastnet_all = disc4_elastnet_all$tn
  report$tn5_elastnet_all = disc5_elastnet_all$tn
  report$tp0_elastnet_all = disc0_elastnet_all$tp #internal val.
  report$tp1_elastnet_all = disc1_elastnet_all$tp 
  report$tp2_elastnet_all = disc2_elastnet_all$tp
  report$tp3_elastnet_all = disc3_elastnet_all$tp
  report$tp4_elastnet_all = disc4_elastnet_all$tp
  report$tp5_elastnet_all = disc5_elastnet_all$tp
  
  report$tn0_elastnet_par = disc0_elastnet_par$tn #internal val.
  report$tn1_elastnet_par = disc1_elastnet_par$tn 
  report$tn2_elastnet_par = disc2_elastnet_par$tn
  report$tn3_elastnet_par = disc3_elastnet_par$tn
  report$tn4_elastnet_par = disc4_elastnet_par$tn
  report$tn5_elastnet_par = disc5_elastnet_par$tn
  report$tp0_elastnet_par = disc0_elastnet_par$tp #internal val.
  report$tp1_elastnet_par = disc1_elastnet_par$tp 
  report$tp2_elastnet_par = disc2_elastnet_par$tp
  report$tp3_elastnet_par = disc3_elastnet_par$tp
  report$tp4_elastnet_par = disc4_elastnet_par$tp
  report$tp5_elastnet_par = disc5_elastnet_par$tp
  
  report$tn0_elastnet_child = disc0_elastnet_child$tn #internal val.
  report$tn1_elastnet_child = disc1_elastnet_child$tn 
  report$tn2_elastnet_child = disc2_elastnet_child$tn
  report$tn3_elastnet_child = disc3_elastnet_child$tn
  report$tn4_elastnet_child = disc4_elastnet_child$tn
  report$tn5_elastnet_child = disc5_elastnet_child$tn
  report$tp0_elastnet_child = disc0_elastnet_child$tp #internal val.
  report$tp1_elastnet_child = disc1_elastnet_child$tp 
  report$tp2_elastnet_child = disc2_elastnet_child$tp
  report$tp3_elastnet_child = disc3_elastnet_child$tp
  report$tp4_elastnet_child = disc4_elastnet_child$tp
  report$tp5_elastnet_child = disc5_elastnet_child$tp
  
  report$tn0_elastnet_exo = disc0_elastnet_exo$tn #internal val.
  report$tn1_elastnet_exo = disc1_elastnet_exo$tn 
  report$tn2_elastnet_exo = disc2_elastnet_exo$tn
  report$tn3_elastnet_exo = disc3_elastnet_exo$tn
  report$tn4_elastnet_exo = disc4_elastnet_exo$tn
  report$tn5_elastnet_exo = disc5_elastnet_exo$tn
  report$tp0_elastnet_exo = disc0_elastnet_exo$tp #internal val.
  report$tp1_elastnet_exo = disc1_elastnet_exo$tp 
  report$tp2_elastnet_exo = disc2_elastnet_exo$tp
  report$tp3_elastnet_exo = disc3_elastnet_exo$tp
  report$tp4_elastnet_exo = disc4_elastnet_exo$tp
  report$tp5_elastnet_exo = disc5_elastnet_exo$tp
  
  ###
  # Random Forest
  ###
  # Models
  rf_all = randomForest::randomForest(formula = formula_all, data = dev_sim_t, keep.forest=TRUE)
  rf_par = randomForest::randomForest(formula = formula_par, data = dev_sim_t, keep.forest=TRUE)
  rf_child = randomForest::randomForest(formula = formula_child, data = dev_sim_t, keep.forest=TRUE)
  rf_exo = randomForest::randomForest(formula = formula_exo, data = dev_sim_t, keep.forest=TRUE)
  
  #for train() to work we have to relevel the dx variable
  levels(dev_sim_t$dx) <- c("zero","one")
  
    
  #Predictions
  pred_rf_all = predict_rf(rf_all, set=all, test_data_predictors=test_sim_t, 
                           inta_predictors=synth_df_inta, inta2_predictors=synth_df_inta2, 
                           intap_predictors=synth_df_intap, ints_predictors=synth_df_ints,
                           inttau_predictors=synth_df_inttau)
  pred_rf_par = predict_rf(rf_par, set=parents, test_data_predictors=test_sim_t, 
                           inta_predictors=synth_df_inta, inta2_predictors=synth_df_inta2, 
                           intap_predictors=synth_df_intap, ints_predictors=synth_df_ints,
                           inttau_predictors=synth_df_inttau)
  pred_rf_child = predict_rf(rf_child, set=children, test_data_predictors=test_sim_t, 
                             inta_predictors=synth_df_inta, inta2_predictors=synth_df_inta2, 
                             intap_predictors=synth_df_intap, ints_predictors=synth_df_ints,
                             inttau_predictors=synth_df_inttau)
  pred_rf_exo = predict_rf(rf_exo, set=exo, test_data_predictors=test_sim_t, 
                             inta_predictors=synth_df_inta, inta2_predictors=synth_df_inta2, 
                             intap_predictors=synth_df_intap, ints_predictors=synth_df_ints,
                             inttau_predictors=synth_df_inttau)
  
  #Evaluation
  report$ICI0_rf_all <- ICI(outcome = test_num[, outcome], prediction = pred_rf_all$pred_0)
  report$ICI1_rf_all <- ICI(outcome = inta_num[, outcome], prediction = pred_rf_all$pred_1)
  report$ICI2_rf_all <- ICI(outcome = intap_num[, outcome], prediction = pred_rf_all$pred_2)
  report$ICI3_rf_all <- ICI(outcome = ints_num[, outcome], prediction = pred_rf_all$pred_3)
  report$ICI4_rf_all <- ICI(outcome = inttau_num[, outcome], prediction = pred_rf_all$pred_4)
  report$ICI5_rf_all <- ICI(outcome = inta2_num[, outcome], prediction = pred_rf_all$pred_5)
  report$Brier0_rf_all = brierscore(test_num[, outcome], pred=pred_rf_all$pred_0)
  report$Brier1_rf_all = brierscore(inta_num[, outcome], pred=pred_rf_all$pred_1)
  report$Brier2_rf_all = brierscore(intap_num[, outcome], pred=pred_rf_all$pred_2)
  report$Brier3_rf_all = brierscore(ints_num[, outcome], pred=pred_rf_all$pred_3)
  report$Brier4_rf_all = brierscore(inttau_num[, outcome], pred=pred_rf_all$pred_4)
  report$Brier5_rf_all = brierscore(inta2_num[, outcome], pred=pred_rf_all$pred_5)
  
  report$ICI0_rf_par <- ICI(outcome = test_num[, outcome], prediction = pred_rf_par$pred_0)
  report$ICI1_rf_par <- ICI(outcome = inta_num[, outcome], prediction = pred_rf_par$pred_1)
  report$ICI2_rf_par <- ICI(outcome = intap_num[, outcome], prediction = pred_rf_par$pred_2)
  report$ICI3_rf_par <- ICI(outcome = ints_num[, outcome], prediction = pred_rf_par$pred_3)
  report$ICI4_rf_par <- ICI(outcome = inttau_num[, outcome], prediction = pred_rf_par$pred_4)
  report$ICI5_rf_par <- ICI(outcome = inta2_num[, outcome], prediction = pred_rf_par$pred_5)
  
  report$Brier0_rf_par = brierscore(test_num[, outcome], pred=pred_rf_par$pred_0)
  report$Brier1_rf_par = brierscore(inta_num[, outcome], pred=pred_rf_par$pred_1)
  report$Brier2_rf_par = brierscore(intap_num[, outcome], pred=pred_rf_par$pred_2)
  report$Brier3_rf_par = brierscore(ints_num[, outcome], pred=pred_rf_par$pred_3)
  report$Brier4_rf_par = brierscore(inttau_num[, outcome], pred=pred_rf_par$pred_4)
  report$Brier5_rf_par = brierscore(inta2_num[, outcome], pred=pred_rf_par$pred_5)
  
  report$ICI0_rf_child <- ICI(outcome = test_num[, outcome], prediction = pred_rf_child$pred_0)
  report$ICI1_rf_child <- ICI(outcome = inta_num[, outcome], prediction = pred_rf_child$pred_1)
  report$ICI2_rf_child <- ICI(outcome = intap_num[, outcome], prediction = pred_rf_child$pred_2)
  report$ICI3_rf_child <- ICI(outcome = ints_num[, outcome], prediction = pred_rf_child$pred_3)
  report$ICI4_rf_child <- ICI(outcome = inttau_num[, outcome], prediction = pred_rf_child$pred_4)
  report$ICI5_rf_child <- ICI(outcome = inta2_num[, outcome], prediction = pred_rf_child$pred_5)
  
  report$Brier0_rf_child = brierscore(test_num[, outcome], pred=pred_rf_child$pred_0)
  report$Brier1_rf_child = brierscore(inta_num[, outcome], pred=pred_rf_child$pred_1)
  report$Brier2_rf_child = brierscore(intap_num[, outcome], pred=pred_rf_child$pred_2)
  report$Brier3_rf_child = brierscore(ints_num[, outcome], pred=pred_rf_child$pred_3)
  report$Brier4_rf_child = brierscore(inttau_num[, outcome], pred=pred_rf_child$pred_4)
  report$Brier5_rf_child = brierscore(inta2_num[, outcome], pred=pred_rf_child$pred_5)
  
  report$ICI0_rf_exo <- ICI(outcome = test_num[, outcome], prediction = pred_rf_exo$pred_0)
  report$ICI1_rf_exo <- ICI(outcome = inta_num[, outcome], prediction = pred_rf_exo$pred_1)
  report$ICI2_rf_exo <- ICI(outcome = intap_num[, outcome], prediction = pred_rf_exo$pred_2)
  report$ICI3_rf_exo <- ICI(outcome = ints_num[, outcome], prediction = pred_rf_exo$pred_3)
  report$ICI4_rf_exo <- ICI(outcome = inttau_num[, outcome], prediction = pred_rf_exo$pred_4)
  report$ICI5_rf_exo <- ICI(outcome = inta2_num[, outcome], prediction = pred_rf_exo$pred_5)
  
  report$Brier0_rf_exo = brierscore(test_num[, outcome], pred=pred_rf_exo$pred_0)
  report$Brier1_rf_exo = brierscore(inta_num[, outcome], pred=pred_rf_exo$pred_1)
  report$Brier2_rf_exo = brierscore(intap_num[, outcome], pred=pred_rf_exo$pred_2)
  report$Brier3_rf_exo = brierscore(ints_num[, outcome], pred=pred_rf_exo$pred_3)
  report$Brier4_rf_exo = brierscore(inttau_num[, outcome], pred=pred_rf_exo$pred_4)
  report$Brier5_rf_exo = brierscore(inta2_num[, outcome], pred=pred_rf_exo$pred_5)
  
  
  #Discrimination metrics
  disc0_rf_all = get_disc(pred_cl = pred_classes(pred_rf_all$pred_0), true_classes=test_num[, outcome])
  disc1_rf_all = get_disc(pred_cl = pred_classes(pred_rf_all$pred_1), true_classes=inta_num[, outcome])
  disc2_rf_all = get_disc(pred_cl = pred_classes(pred_rf_all$pred_2), true_classes=intap_num[, outcome])
  disc3_rf_all = get_disc(pred_cl = pred_classes(pred_rf_all$pred_3), true_classes=ints_num[, outcome])
  disc4_rf_all = get_disc(pred_cl = pred_classes(pred_rf_all$pred_4), true_classes=inttau_num[, outcome])
  disc5_rf_all = get_disc(pred_cl = pred_classes(pred_rf_all$pred_5), true_classes=inta2_num[, outcome])
  disc0_rf_par = get_disc(pred_cl = pred_classes(pred_rf_par$pred_0), true_classes=test_num[, outcome])
  disc1_rf_par = get_disc(pred_cl = pred_classes(pred_rf_par$pred_1), true_classes=inta_num[, outcome])
  disc2_rf_par = get_disc(pred_cl = pred_classes(pred_rf_par$pred_2), true_classes=intap_num[, outcome])
  disc3_rf_par = get_disc(pred_cl = pred_classes(pred_rf_par$pred_3), true_classes=ints_num[, outcome])
  disc4_rf_par = get_disc(pred_cl = pred_classes(pred_rf_par$pred_4), true_classes=inttau_num[, outcome])
  disc5_rf_par = get_disc(pred_cl = pred_classes(pred_rf_par$pred_5), true_classes=inta2_num[, outcome])
  disc0_rf_child = get_disc(pred_cl = pred_classes(pred_rf_child$pred_0), true_classes=test_num[, outcome])
  disc1_rf_child = get_disc(pred_cl = pred_classes(pred_rf_child$pred_1), true_classes=inta_num[, outcome])
  disc2_rf_child = get_disc(pred_cl = pred_classes(pred_rf_child$pred_2), true_classes=intap_num[, outcome])
  disc3_rf_child = get_disc(pred_cl = pred_classes(pred_rf_child$pred_3), true_classes=ints_num[, outcome])
  disc4_rf_child = get_disc(pred_cl = pred_classes(pred_rf_child$pred_4), true_classes=inttau_num[, outcome])
  disc5_rf_child = get_disc(pred_cl = pred_classes(pred_rf_child$pred_5), true_classes=inta2_num[, outcome])
  disc0_rf_exo = get_disc(pred_cl = pred_classes(pred_rf_exo$pred_0), true_classes=test_num[, outcome])
  disc1_rf_exo = get_disc(pred_cl = pred_classes(pred_rf_exo$pred_1), true_classes=inta_num[, outcome])
  disc2_rf_exo = get_disc(pred_cl = pred_classes(pred_rf_exo$pred_2), true_classes=intap_num[, outcome])
  disc3_rf_exo = get_disc(pred_cl = pred_classes(pred_rf_exo$pred_3), true_classes=ints_num[, outcome])
  disc4_rf_exo = get_disc(pred_cl = pred_classes(pred_rf_exo$pred_4), true_classes=inttau_num[, outcome])
  disc5_rf_exo = get_disc(pred_cl = pred_classes(pred_rf_exo$pred_5), true_classes=inta2_num[, outcome])
  
  
  #AUC
  report$auc0_rf_all = disc0_rf_all$auc #internal val.
  report$auc1_rf_all = disc1_rf_all$auc 
  report$auc2_rf_all = disc2_rf_all$auc
  report$auc3_rf_all = disc3_rf_all$auc
  report$auc4_rf_all = disc4_rf_all$auc
  report$auc5_rf_all = disc5_rf_all$auc
  report$auc0_rf_par = disc0_rf_par$auc #internal val.
  report$auc1_rf_par = disc1_rf_par$auc 
  report$auc2_rf_par = disc2_rf_par$auc
  report$auc3_rf_par = disc3_rf_par$auc
  report$auc4_rf_par = disc4_rf_par$auc
  report$auc5_rf_par = disc5_rf_par$auc
  report$auc0_rf_child = disc0_rf_child$auc #internal val.
  report$auc1_rf_child = disc1_rf_child$auc 
  report$auc2_rf_child = disc2_rf_child$auc
  report$auc3_rf_child = disc3_rf_child$auc
  report$auc4_rf_child = disc4_rf_child$auc
  report$auc5_rf_child = disc5_rf_child$auc
  report$auc0_rf_exo = disc0_rf_exo$auc #internal val.
  report$auc1_rf_exo = disc1_rf_exo$auc 
  report$auc2_rf_exo = disc2_rf_exo$auc
  report$auc3_rf_exo = disc3_rf_exo$auc
  report$auc4_rf_exo = disc4_rf_exo$auc
  report$auc5_rf_exo = disc5_rf_exo$auc
  
  
  #Balacc & F1
  report$balacc0_rf_all = disc0_rf_all$balacc #internal val.
  report$balacc1_rf_all = disc1_rf_all$balacc 
  report$balacc2_rf_all = disc2_rf_all$balacc
  report$balacc3_rf_all = disc3_rf_all$balacc
  report$balacc4_rf_all = disc4_rf_all$balacc
  report$balacc5_rf_all = disc5_rf_all$balacc
  report$Fone0_rf_all = disc0_rf_all$f1 #internal val.
  report$Fone1_rf_all = disc1_rf_all$f1 
  report$Fone2_rf_all = disc2_rf_all$f1
  report$Fone3_rf_all = disc3_rf_all$f1
  report$Fone4_rf_all = disc4_rf_all$f1
  report$Fone5_rf_all = disc5_rf_all$f1
  
  report$tn0_rf_all = disc0_rf_all$tn #internal val.
  report$tn1_rf_all = disc1_rf_all$tn 
  report$tn2_rf_all = disc2_rf_all$tn
  report$tn3_rf_all = disc3_rf_all$tn
  report$tn4_rf_all = disc4_rf_all$tn
  report$tn5_rf_all = disc5_rf_all$tn
  report$tp0_rf_all = disc0_rf_all$tp #internal val.
  report$tp1_rf_all = disc1_rf_all$tp 
  report$tp2_rf_all = disc2_rf_all$tp
  report$tp3_rf_all = disc3_rf_all$tp
  report$tp4_rf_all = disc4_rf_all$tp
  report$tp5_rf_all = disc5_rf_all$tp
  
  report$balacc0_rf_par = disc0_rf_par$balacc #internal val.
  report$balacc1_rf_par = disc1_rf_par$balacc 
  report$balacc2_rf_par = disc2_rf_par$balacc
  report$balacc3_rf_par = disc3_rf_par$balacc
  report$balacc4_rf_par = disc4_rf_par$balacc
  report$balacc5_rf_par = disc5_rf_par$balacc
  report$Fone0_rf_par = disc0_rf_par$f1 #internal val.
  report$Fone1_rf_par = disc1_rf_par$f1 
  report$Fone2_rf_par = disc2_rf_par$f1
  report$Fone3_rf_par = disc3_rf_par$f1
  report$Fone4_rf_par = disc4_rf_par$f1
  report$Fone5_rf_par = disc5_rf_par$f1
  
  report$tn0_rf_par = disc0_rf_par$tn #internal val.
  report$tn1_rf_par = disc1_rf_par$tn 
  report$tn2_rf_par = disc2_rf_par$tn
  report$tn3_rf_par = disc3_rf_par$tn
  report$tn4_rf_par = disc4_rf_par$tn
  report$tn5_rf_par = disc5_rf_par$tn
  report$tp0_rf_par = disc0_rf_par$tp #internal val.
  report$tp1_rf_par = disc1_rf_par$tp 
  report$tp2_rf_par = disc2_rf_par$tp
  report$tp3_rf_par = disc3_rf_par$tp
  report$tp4_rf_par = disc4_rf_par$tp
  report$tp5_rf_par = disc5_rf_par$tp
  
  
  report$balacc0_rf_child = disc0_rf_child$balacc #internal val.
  report$balacc1_rf_child = disc1_rf_child$balacc 
  report$balacc2_rf_child = disc2_rf_child$balacc
  report$balacc3_rf_child = disc3_rf_child$balacc
  report$balacc4_rf_child = disc4_rf_child$balacc
  report$balacc5_rf_child = disc5_rf_child$balacc
  report$Fone0_rf_child = disc0_rf_child$f1 #internal val.
  report$Fone1_rf_child = disc1_rf_child$f1 
  report$Fone2_rf_child = disc2_rf_child$f1
  report$Fone3_rf_child = disc3_rf_child$f1
  report$Fone4_rf_child = disc4_rf_child$f1
  report$Fone5_rf_child = disc5_rf_child$f1
  
  report$tn0_rf_child = disc0_rf_child$tn #internal val.
  report$tn1_rf_child = disc1_rf_child$tn 
  report$tn2_rf_child = disc2_rf_child$tn
  report$tn3_rf_child = disc3_rf_child$tn
  report$tn4_rf_child = disc4_rf_child$tn
  report$tn5_rf_child = disc5_rf_child$tn
  report$tp0_rf_child = disc0_rf_child$tp #internal val.
  report$tp1_rf_child = disc1_rf_child$tp 
  report$tp2_rf_child = disc2_rf_child$tp
  report$tp3_rf_child = disc3_rf_child$tp
  report$tp4_rf_child = disc4_rf_child$tp
  report$tp5_rf_child = disc5_rf_child$tp
  
  
  report$balacc0_rf_exo = disc0_rf_exo$balacc #internal val.
  report$balacc1_rf_exo = disc1_rf_exo$balacc 
  report$balacc2_rf_exo = disc2_rf_exo$balacc
  report$balacc3_rf_exo = disc3_rf_exo$balacc
  report$balacc4_rf_exo = disc4_rf_exo$balacc
  report$balacc5_rf_exo = disc5_rf_exo$balacc
  report$Fone0_rf_exo = disc0_rf_exo$f1 #internal val.
  report$Fone1_rf_exo = disc1_rf_exo$f1 
  report$Fone2_rf_exo = disc2_rf_exo$f1
  report$Fone3_rf_exo = disc3_rf_exo$f1
  report$Fone4_rf_exo = disc4_rf_exo$f1
  report$Fone5_rf_exo = disc5_rf_exo$f1
  
  report$tn0_rf_exo = disc0_rf_exo$tn #internal val.
  report$tn1_rf_exo = disc1_rf_exo$tn 
  report$tn2_rf_exo = disc2_rf_exo$tn
  report$tn3_rf_exo = disc3_rf_exo$tn
  report$tn4_rf_exo = disc4_rf_exo$tn
  report$tn5_rf_exo = disc5_rf_exo$tn
  report$tp0_rf_exo = disc0_rf_exo$tp #internal val.
  report$tp1_rf_exo = disc1_rf_exo$tp 
  report$tp2_rf_exo = disc2_rf_exo$tp
  report$tp3_rf_exo = disc3_rf_exo$tp
  report$tp4_rf_exo = disc4_rf_exo$tp
  report$tp5_rf_exo = disc5_rf_exo$tp
  
  
  ###
  # GBM
  ###
  # Models
  gbmMod_all <- gbm(formula_all, data=dev_num, distribution = "bernoulli", bag.fraction=1)
  gbmMod_par <- gbm(formula_par, data=dev_num, distribution = "bernoulli", bag.fraction=1)
  gbmMod_child <- gbm(formula_child, data=dev_num, distribution = "bernoulli", bag.fraction=1)
  gbmMod_exo <- gbm(formula_exo, data=dev_num, distribution = "bernoulli", bag.fraction=1)
   
  #Predictions
  pred_gbm_all = predict_gbm(gbm_mod=gbmMod_all, set=all, 
                              train_data_predictors=dev_num,
                              test_data_predictors=test_num, 
                              inta_predictors=inta_num, inta2_predictors=inta2_num, 
                              intap_predictors=intap_num, ints_predictors=ints_num, 
                              inttau_predictors=inttau_num)
   
  pred_gbm_parents = predict_gbm(gbm_mod=gbmMod_par, set=parents, 
                                 train_data_predictors=dev_num,
                                 test_data_predictors=test_num, 
                                  inta_predictors=inta_num, inta2_predictors=inta2_num, 
                                  intap_predictors=intap_num, ints_predictors=ints_num, 
                                 inttau_predictors=inttau_num)
   
  pred_gbm_child = predict_gbm(gbm_mod=gbmMod_child, set=children, 
                                train_data_predictors=dev_num,
                                test_data_predictors=test_num, 
                                inta_predictors=inta_num, inta2_predictors=inta2_num, 
                                intap_predictors=intap_num, ints_predictors=ints_num, 
                               inttau_predictors=inttau_num)
  
  pred_gbm_exo = predict_gbm(gbm_mod=gbmMod_exo, set=exo, 
                               train_data_predictors=dev_num,
                               test_data_predictors=test_num, 
                               inta_predictors=inta_num, inta2_predictors=inta2_num, 
                               intap_predictors=intap_num, ints_predictors=ints_num, 
                               inttau_predictors=inttau_num)
  
  #Validation on internal and external data
  report$ICI0_gbm_all <-ICI(outcome = test_num[, outcome], prediction = pred_gbm_all$pred_0)
  report$ICI1_gbm_all <-ICI(outcome = inta_num[, outcome], prediction = pred_gbm_all$pred_1)
  report$ICI2_gbm_all <-ICI(outcome = intap_num[, outcome], prediction = pred_gbm_all$pred_2)
  report$ICI3_gbm_all <-ICI(outcome = ints_num[, outcome], prediction = pred_gbm_all$pred_3)
  report$ICI4_gbm_all <-ICI(outcome = inttau_num[, outcome], prediction = pred_gbm_all$pred_4)
  report$ICI5_gbm_all <-ICI(outcome = inta2_num[, outcome], prediction = pred_gbm_all$pred_5)
  report$Brier0_gbm_all = brierscore(test_num[, outcome], pred=pred_gbm_all$pred_0)
  report$Brier1_gbm_all = brierscore(inta_num[, outcome], pred=pred_gbm_all$pred_1)
  report$Brier2_gbm_all = brierscore(intap_num[, outcome], pred=pred_gbm_all$pred_2)
  report$Brier3_gbm_all = brierscore(ints_num[, outcome], pred=pred_gbm_all$pred_3)
  report$Brier4_gbm_all = brierscore(inttau_num[, outcome], pred=pred_gbm_all$pred_4)
  report$Brier5_gbm_all = brierscore(inta2_num[, outcome], pred=pred_gbm_all$pred_5)
  
  report$ICI0_gbm_par <-ICI(outcome = test_num[, outcome], prediction = pred_gbm_parents$pred_0)
  report$ICI1_gbm_par <-ICI(outcome = inta_num[, outcome], prediction = pred_gbm_parents$pred_1)
  report$ICI2_gbm_par <-ICI(outcome = intap_num[, outcome], prediction = pred_gbm_parents$pred_2)
  report$ICI3_gbm_par <-ICI(outcome = ints_num[, outcome], prediction = pred_gbm_parents$pred_3)
  report$ICI4_gbm_par <-ICI(outcome = inttau_num[, outcome], prediction = pred_gbm_parents$pred_4)
  report$ICI5_gbm_par <-ICI(outcome = inta2_num[, outcome], prediction = pred_gbm_parents$pred_5)
  report$Brier0_gbm_par = brierscore(test_num[, outcome], pred=pred_gbm_parents$pred_0)
  report$Brier1_gbm_par = brierscore(inta_num[, outcome], pred=pred_gbm_parents$pred_1)
  report$Brier2_gbm_par = brierscore(intap_num[, outcome], pred=pred_gbm_parents$pred_2)
  report$Brier3_gbm_par = brierscore(ints_num[, outcome], pred=pred_gbm_parents$pred_3)
  report$Brier4_gbm_par = brierscore(inttau_num[, outcome], pred=pred_gbm_parents$pred_4)
  report$Brier5_gbm_par = brierscore(inta2_num[, outcome], pred=pred_gbm_parents$pred_5)
   
  report$ICI0_gbm_child <-ICI(outcome = test_num[, outcome], prediction = pred_gbm_child$pred_0)
  report$ICI1_gbm_child <-ICI(outcome = inta_num[, outcome], prediction = pred_gbm_child$pred_1)
  report$ICI2_gbm_child <-ICI(outcome = intap_num[, outcome], prediction = pred_gbm_child$pred_2)
  report$ICI3_gbm_child <-ICI(outcome = ints_num[, outcome], prediction = pred_gbm_child$pred_3)
  report$ICI4_gbm_child <-ICI(outcome = inttau_num[, outcome], prediction = pred_gbm_child$pred_4)
  report$ICI5_gbm_child <-ICI(outcome = inta2_num[, outcome], prediction = pred_gbm_child$pred_5)
  report$Brier0_gbm_child = brierscore(test_num[, outcome], pred=pred_gbm_child$pred_0)
  report$Brier1_gbm_child = brierscore(inta_num[, outcome], pred=pred_gbm_child$pred_1)
  report$Brier2_gbm_child = brierscore(intap_num[, outcome], pred=pred_gbm_child$pred_2)
  report$Brier3_gbm_child = brierscore(ints_num[, outcome], pred=pred_gbm_child$pred_3)
  report$Brier4_gbm_child = brierscore(inttau_num[, outcome], pred=pred_gbm_child$pred_4)
  report$Brier5_gbm_child = brierscore(inta2_num[, outcome], pred=pred_gbm_child$pred_5)
  
  report$ICI0_gbm_exo <-ICI(outcome = test_num[, outcome], prediction = pred_gbm_exo$pred_0)
  report$ICI1_gbm_exo <-ICI(outcome = inta_num[, outcome], prediction = pred_gbm_exo$pred_1)
  report$ICI2_gbm_exo <-ICI(outcome = intap_num[, outcome], prediction = pred_gbm_exo$pred_2)
  report$ICI3_gbm_exo <-ICI(outcome = ints_num[, outcome], prediction = pred_gbm_exo$pred_3)
  report$ICI4_gbm_exo <-ICI(outcome = inttau_num[, outcome], prediction = pred_gbm_exo$pred_4)
  report$ICI5_gbm_exo <-ICI(outcome = inta2_num[, outcome], prediction = pred_gbm_exo$pred_5)
  report$Brier0_gbm_exo = brierscore(test_num[, outcome], pred=pred_gbm_exo$pred_0)
  report$Brier1_gbm_exo = brierscore(inta_num[, outcome], pred=pred_gbm_exo$pred_1)
  report$Brier2_gbm_exo = brierscore(intap_num[, outcome], pred=pred_gbm_exo$pred_2)
  report$Brier3_gbm_exo = brierscore(ints_num[, outcome], pred=pred_gbm_exo$pred_3)
  report$Brier4_gbm_exo = brierscore(inttau_num[, outcome], pred=pred_gbm_exo$pred_4)
  report$Brier5_gbm_exo = brierscore(inta2_num[, outcome], pred=pred_gbm_exo$pred_5)
  
  
  #Discrimination metrics
  disc0_gbm_all = get_disc(pred_cl = pred_classes(pred_gbm_all$pred_0), true_classes=test_num[, outcome])
  disc1_gbm_all = get_disc(pred_cl = pred_classes(pred_gbm_all$pred_1), true_classes=inta_num[, outcome])
  disc2_gbm_all = get_disc(pred_cl = pred_classes(pred_gbm_all$pred_2), true_classes=intap_num[, outcome])
  disc3_gbm_all = get_disc(pred_cl = pred_classes(pred_gbm_all$pred_3), true_classes=ints_num[, outcome])
  disc4_gbm_all = get_disc(pred_cl = pred_classes(pred_gbm_all$pred_4), true_classes=inttau_num[, outcome])
  disc5_gbm_all = get_disc(pred_cl = pred_classes(pred_gbm_all$pred_5), true_classes=inta2_num[, outcome])
  disc0_gbm_par = get_disc(pred_cl = pred_classes(pred_gbm_parents$pred_0), true_classes=test_num[, outcome])
  disc1_gbm_par = get_disc(pred_cl = pred_classes(pred_gbm_parents$pred_1), true_classes=inta_num[, outcome])
  disc2_gbm_par = get_disc(pred_cl = pred_classes(pred_gbm_parents$pred_2), true_classes=intap_num[, outcome])
  disc3_gbm_par = get_disc(pred_cl = pred_classes(pred_gbm_parents$pred_3), true_classes=ints_num[, outcome])
  disc4_gbm_par = get_disc(pred_cl = pred_classes(pred_gbm_parents$pred_4), true_classes=inttau_num[, outcome])
  disc5_gbm_par = get_disc(pred_cl = pred_classes(pred_gbm_parents$pred_5), true_classes=inta2_num[, outcome])
  disc0_gbm_child = get_disc(pred_cl = pred_classes(pred_gbm_child$pred_0), true_classes=test_num[, outcome])
  disc1_gbm_child = get_disc(pred_cl = pred_classes(pred_gbm_child$pred_1), true_classes=inta_num[, outcome])
  disc2_gbm_child = get_disc(pred_cl = pred_classes(pred_gbm_child$pred_2), true_classes=intap_num[, outcome])
  disc3_gbm_child = get_disc(pred_cl = pred_classes(pred_gbm_child$pred_3), true_classes=ints_num[, outcome])
  disc4_gbm_child = get_disc(pred_cl = pred_classes(pred_gbm_child$pred_4), true_classes=inttau_num[, outcome])
  disc5_gbm_child = get_disc(pred_cl = pred_classes(pred_gbm_child$pred_5), true_classes=inta2_num[, outcome])
  disc0_gbm_exo = get_disc(pred_cl = pred_classes(pred_gbm_exo$pred_0), true_classes=test_num[, outcome])
  disc1_gbm_exo = get_disc(pred_cl = pred_classes(pred_gbm_exo$pred_1), true_classes=inta_num[, outcome])
  disc2_gbm_exo = get_disc(pred_cl = pred_classes(pred_gbm_exo$pred_2), true_classes=intap_num[, outcome])
  disc3_gbm_exo = get_disc(pred_cl = pred_classes(pred_gbm_exo$pred_3), true_classes=ints_num[, outcome])
  disc4_gbm_exo = get_disc(pred_cl = pred_classes(pred_gbm_exo$pred_4), true_classes=inttau_num[, outcome])
  disc5_gbm_exo = get_disc(pred_cl = pred_classes(pred_gbm_exo$pred_5), true_classes=inta2_num[, outcome])
  
  #AUC
  report$auc0_gbm_all = disc0_gbm_all$auc #internal val.
  report$auc1_gbm_all = disc1_gbm_all$auc 
  report$auc2_gbm_all = disc2_gbm_all$auc
  report$auc3_gbm_all = disc3_gbm_all$auc
  report$auc4_gbm_all = disc4_gbm_all$auc
  report$auc5_gbm_all = disc5_gbm_all$auc
  report$auc0_gbm_par = disc0_gbm_par$auc #internal val.
  report$auc1_gbm_par = disc1_gbm_par$auc 
  report$auc2_gbm_par = disc2_gbm_par$auc
  report$auc3_gbm_par = disc3_gbm_par$auc
  report$auc4_gbm_par = disc4_gbm_par$auc
  report$auc5_gbm_par = disc5_gbm_par$auc
  report$auc0_gbm_child = disc0_gbm_child$auc #internal val.
  report$auc1_gbm_child = disc1_gbm_child$auc 
  report$auc2_gbm_child = disc2_gbm_child$auc
  report$auc3_gbm_child = disc3_gbm_child$auc
  report$auc4_gbm_child = disc4_gbm_child$auc
  report$auc5_gbm_child = disc5_gbm_child$auc
  report$auc0_gbm_exo = disc0_gbm_exo$auc #internal val.
  report$auc1_gbm_exo = disc1_gbm_exo$auc 
  report$auc2_gbm_exo = disc2_gbm_exo$auc
  report$auc3_gbm_exo = disc3_gbm_exo$auc
  report$auc4_gbm_exo = disc4_gbm_exo$auc
  report$auc5_gbm_exo = disc5_gbm_exo$auc
   
  #Balacc & F1
  report$balacc0_gbm_all = disc0_gbm_all$balacc #internal val.
  report$balacc1_gbm_all = disc1_gbm_all$balacc 
  report$balacc2_gbm_all = disc2_gbm_all$balacc
  report$balacc3_gbm_all = disc3_gbm_all$balacc
  report$balacc4_gbm_all = disc4_gbm_all$balacc
  report$balacc5_gbm_all = disc5_gbm_all$balacc
  report$Fone0_gbm_all = disc0_gbm_all$f1 #internal val.
  report$Fone1_gbm_all = disc1_gbm_all$f1 
  report$Fone2_gbm_all = disc2_gbm_all$f1
  report$Fone3_gbm_all = disc3_gbm_all$f1
  report$Fone4_gbm_all = disc4_gbm_all$f1
  report$Fone5_gbm_all = disc5_gbm_all$f1
  
  report$balacc0_gbm_par = disc0_gbm_par$balacc #internal val.
  report$balacc1_gbm_par = disc1_gbm_par$balacc 
  report$balacc2_gbm_par = disc2_gbm_par$balacc
  report$balacc3_gbm_par = disc3_gbm_par$balacc
  report$balacc4_gbm_par = disc4_gbm_par$balacc
  report$balacc5_gbm_par = disc5_gbm_par$balacc
  report$Fone0_gbm_par = disc0_gbm_par$f1 #internal val.
  report$Fone1_gbm_par = disc1_gbm_par$f1 
  report$Fone2_gbm_par = disc2_gbm_par$f1
  report$Fone3_gbm_par = disc3_gbm_par$f1
  report$Fone4_gbm_par = disc4_gbm_par$f1
  report$Fone5_gbm_par = disc5_gbm_par$f1
  
  report$balacc0_gbm_child = disc0_gbm_child$balacc #internal val.
  report$balacc1_gbm_child = disc1_gbm_child$balacc 
  report$balacc2_gbm_child = disc2_gbm_child$balacc
  report$balacc3_gbm_child = disc3_gbm_child$balacc
  report$balacc4_gbm_child = disc4_gbm_child$balacc
  report$balacc5_gbm_child = disc5_gbm_child$balacc
  report$Fone0_gbm_child = disc0_gbm_child$f1 #internal val.
  report$Fone1_gbm_child = disc1_gbm_child$f1 
  report$Fone2_gbm_child = disc2_gbm_child$f1
  report$Fone3_gbm_child = disc3_gbm_child$f1
  report$Fone4_gbm_child = disc4_gbm_child$f1
  report$Fone5_gbm_child = disc5_gbm_child$f1
  
  report$balacc0_gbm_exo = disc0_gbm_exo$balacc #internal val.
  report$balacc1_gbm_exo = disc1_gbm_exo$balacc 
  report$balacc2_gbm_exo = disc2_gbm_exo$balacc
  report$balacc3_gbm_exo = disc3_gbm_exo$balacc
  report$balacc4_gbm_exo = disc4_gbm_exo$balacc
  report$balacc5_gbm_exo = disc5_gbm_exo$balacc
  report$Fone0_gbm_exo = disc0_gbm_exo$f1 #internal val.
  report$Fone1_gbm_exo = disc1_gbm_exo$f1 
  report$Fone2_gbm_exo = disc2_gbm_exo$f1
  report$Fone3_gbm_exo = disc3_gbm_exo$f1
  report$Fone4_gbm_exo = disc4_gbm_exo$f1
  report$Fone5_gbm_exo = disc5_gbm_exo$f1
  
  
  report$tn0_gbm_all = disc0_gbm_all$tn #internal val.
  report$tn1_gbm_all = disc1_gbm_all$tn 
  report$tn2_gbm_all = disc2_gbm_all$tn
  report$tn3_gbm_all = disc3_gbm_all$tn
  report$tn4_gbm_all = disc4_gbm_all$tn
  report$tn5_gbm_all = disc5_gbm_all$tn
  report$tp0_gbm_all = disc0_gbm_all$tp #internal val.
  report$tp1_gbm_all = disc1_gbm_all$tp 
  report$tp2_gbm_all = disc2_gbm_all$tp
  report$tp3_gbm_all = disc3_gbm_all$tp
  report$tp4_gbm_all = disc4_gbm_all$tp
  report$tp5_gbm_all = disc5_gbm_all$tp

  report$tn0_gbm_par = disc0_gbm_par$tn #internal val.
  report$tn1_gbm_par = disc1_gbm_par$tn 
  report$tn2_gbm_par = disc2_gbm_par$tn
  report$tn3_gbm_par = disc3_gbm_par$tn
  report$tn4_gbm_par = disc4_gbm_par$tn
  report$tn5_gbm_par = disc5_gbm_par$tn
  report$tp0_gbm_par = disc0_gbm_par$tp #internal val.
  report$tp1_gbm_par = disc1_gbm_par$tp 
  report$tp2_gbm_par = disc2_gbm_par$tp
  report$tp3_gbm_par = disc3_gbm_par$tp
  report$tp4_gbm_par = disc4_gbm_par$tp
  report$tp5_gbm_par = disc5_gbm_par$tp
  
  report$tn0_gbm_child = disc0_gbm_child$tn #internal val.
  report$tn1_gbm_child = disc1_gbm_child$tn 
  report$tn2_gbm_child = disc2_gbm_child$tn
  report$tn3_gbm_child = disc3_gbm_child$tn
  report$tn4_gbm_child = disc4_gbm_child$tn
  report$tn5_gbm_child = disc5_gbm_child$tn
  report$tp0_gbm_child = disc0_gbm_child$tp #internal val.
  report$tp1_gbm_child = disc1_gbm_child$tp 
  report$tp2_gbm_child = disc2_gbm_child$tp
  report$tp3_gbm_child = disc3_gbm_child$tp
  report$tp4_gbm_child = disc4_gbm_child$tp
  report$tp5_gbm_child = disc5_gbm_child$tp
  
  report$tn0_gbm_exo = disc0_gbm_exo$tn #internal val.
  report$tn1_gbm_exo = disc1_gbm_exo$tn 
  report$tn2_gbm_exo = disc2_gbm_exo$tn
  report$tn3_gbm_exo = disc3_gbm_exo$tn
  report$tn4_gbm_exo = disc4_gbm_exo$tn
  report$tn5_gbm_exo = disc5_gbm_exo$tn
  report$tp0_gbm_exo = disc0_gbm_exo$tp #internal val.
  report$tp1_gbm_exo = disc1_gbm_exo$tp 
  report$tp2_gbm_exo = disc2_gbm_exo$tp
  report$tp3_gbm_exo = disc3_gbm_exo$tp
  report$tp4_gbm_exo = disc4_gbm_exo$tp
  report$tp5_gbm_exo = disc5_gbm_exo$tp
  
  return(report)
  
} 
stopCluster(cl)


write.csv(report, file = paste0(output_path, 'report_', Sys.Date(), Sys.time(), '.csv'))

print('Done...\n')

