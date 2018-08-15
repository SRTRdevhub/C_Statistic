


sim_c_effect_iter_fun <- function(curr_iter,
                                  lp_var_inflate, 
                                  ctr_effect_var_inflate,
                                  ctr_effect_confounder_prop = NULL,
                                  program_var_proportion = NULL){
  
  
  if(! is.null(ctr_effect_confounder_prop)){
    if({ctr_effect_confounder_prop < 0} | {ctr_effect_confounder_prop > 1}){ stop("A non-null ctr_effect_confounder_prop should be 0-1") }
  }
  
  
  ### The number of simulated transplants is equal to the sum of transplant program sizes
  tx_n <- sum(c_effect_est_params$program_size)
  
  
  ### If the program_var_proportion is not NULL, then redistribute the overall variance
  if(is.null(program_var_proportion)){
    ctr_var <- c_effect_est_params$var_components[1]
    px_var  <- c_effect_est_params$var_components[2]
  }
  if({! is.null(program_var_proportion)}){
    if({{program_var_proportion < 0} | {program_var_proportion > 1}}){
      stop("Invalid program_var_proportion")
    }
    if({program_var_proportion >= 0} & {program_var_proportion <= 1}){
      ctr_var <- program_var_proportion * sum(c_effect_est_params$var_components)
      px_var  <- (1 - program_var_proportion) * sum(c_effect_est_params$var_components)
    }
  }

  
  
  ### Setting the seed for the current iteration
  set.seed(3 * curr_iter + 2)
  
  
  
  ### Step 1: Simulate the components of the linear predictors
  curr_lp_program_component <- rnorm(length(c_effect_est_params$program_size),
                                     mean = 0,
                                     sd = lp_var_inflate * sqrt(ctr_var))
  
  curr_lp_patient_component <- rnorm(tx_n,
                                     mean = 0,
                                     sd = lp_var_inflate * sqrt(px_var))
  
  
  # If there is unmeasured confounding, proportion the program variance
  if(is.null(ctr_effect_confounder_prop)){
    curr_program_effect <- rnorm(length(c_effect_est_params$program_size),
                                 mean = 0,
                                 sd = sqrt(ctr_effect_var_inflate * c_effect_est_params$program_var))
  }
  if(! is.null(ctr_effect_confounder_prop)){
    curr_program_effect <- rnorm(length(c_effect_est_params$program_size),
                                 mean = 0,
                                 sd = sqrt((1 - ctr_effect_confounder_prop) * 
                                             ctr_effect_var_inflate * 
                                             c_effect_est_params$program_var))
    
    curr_confound_effect <- rnorm(length(c_effect_est_params$program_size),
                                  mean = 0,
                                  sd = sqrt(ctr_effect_confounder_prop * 
                                              ctr_effect_var_inflate * 
                                              c_effect_est_params$program_var))
  }

  
  
  ### Step 2: Expand the program indicators to the length of the patient components
  curr_programs <- rep(1:length(c_effect_est_params$program_size),
                       times = c_effect_est_params$program_size)
  
  long_lp_program_component <- rep(curr_lp_program_component,
                                   times = c_effect_est_params$program_size)
  
  long_program_effect <- rep(curr_program_effect,
                             times = c_effect_est_params$program_size)
  
  # If it exists, need to track the confounder program effects
  if(! is.null(ctr_effect_confounder_prop)){
    long_confound_effect <- rep(curr_confound_effect,
                                times = c_effect_est_params$program_size)
  }
  
  
  ### Step 3: Determine the scaler to add to the linear predictor to ensure a similar number of expected events
  bl_scaler <- uniroot(function(mu){
    integrate(function(x){
      dnorm(x, 
            0,
            sd = sqrt((lp_var_inflate ^ 2) * ctr_var + (lp_var_inflate ^ 2) * px_var + ctr_effect_var_inflate * c_effect_est_params$program_var)) *
        exp(-c_effect_est_params$baseline_hazard$cumulative_hazard[366] * exp(x + mu))
    },
    lower = -Inf,
    upper = Inf)$value - exp(-c_effect_est_params$baseline_hazard$cumulative_hazard[366])
  },
  interval = c(-100, 100))
  
  
  
  ### Step 4: Combine the program-level and patient-level linear predictors with the program-specific HR. Depends on existance of confounder
  if(is.null(ctr_effect_confounder_prop)){
    curr_lp <- long_lp_program_component + curr_lp_patient_component + long_program_effect + bl_scaler$root
  }
  if(! is.null(ctr_effect_confounder_prop)){
    curr_lp <- long_lp_program_component + curr_lp_patient_component + long_program_effect + long_confound_effect + bl_scaler$root
  }
  
  
  ### Step 5: Generate survival times and events
  random_uniform_draws <- runif(tx_n,
                                0,
                                1)
  
  curr_surv <- sapply(1:tx_n,
                      function(x){
                        gen_surv_time(random_uniform_draws[x],
                                      c_effect_est_params$baseline_hazard$cumulative_hazard * exp(curr_lp[x]),
                                      c_effect_est_params$baseline_hazard$time)
                      })
  
  
  ### Step 6: Creating the data set but randomly censoring some observations at six months
  six_month_cens <- rbinom(tx_n,
                           1, 
                           0.2)
  
  
  
  curr_data <- data.frame(observed_time = ifelse(six_month_cens == 1,
                                                 ifelse(curr_surv[1,] <= 183,
                                                        curr_surv[1,],
                                                        183),
                                                 curr_surv[1,]),
                          observed_event = ifelse(six_month_cens == 1,
                                                  ifelse(curr_surv[1,] <= 183,
                                                         curr_surv[2,],
                                                         0),
                                                  curr_surv[2,]),
                          six_month_censor = six_month_cens,
                          observed_lp = long_lp_program_component + curr_lp_patient_component,
                          observed_ctr = curr_programs,
                          true_program_effect = long_program_effect)
  
  
  #########################################################################################################
  #########################################################################################################
  #########################################################################################################
  
  ### Calculating the performance metrics
  
  curr_model <- coxph(Surv(observed_time, observed_event) ~ observed_lp,
                      data = curr_data)
  
  
  curr_km <- survfit(Surv(observed_time, observed_event) ~ 1,
                     data = curr_data)
  
  
  pred_data <- data.frame(observed = curr_data$observed_event,
                          expected = predict(curr_model, type = "expected"),
                          ctr = curr_data$observed_ctr,
                          true_HR = exp(curr_data$true_program_effect))
  

  ctr_summary <- pred_data %>%
    group_by(ctr) %>%
    summarise(ctr_obs = sum(observed),
              ctr_exp = sum(expected),
              true_hr = first(true_HR)) %>%
    mutate(HR_unshr = ctr_obs / ctr_exp,
           HR_shr = (ctr_obs + 2) / (ctr_exp + 2),
           obs_mpsc_flag_c1 = pgamma(1.2, ctr_obs + 2, ctr_exp + 2, lower.tail = FALSE) > 0.75,
           obs_mpsc_flag_c2 = pgamma(2.5, ctr_obs + 2, ctr_exp + 2, lower.tail = FALSE) > 0.10,
           obs_mpsc_flag = 1 * {obs_mpsc_flag_c1 | obs_mpsc_flag_c2},
           obs_cms_flag = 1 * {{(ctr_obs - ctr_exp) > 3} & 
               {ctr_obs / ctr_exp} > 1.85 & 
               {(1 - ppois(ctr_obs - 1, ctr_exp)) < 0.05}})
  

  ctr_summary$obs_5tier_score <- sapply(1:nrow(ctr_summary),
                                        function(x){ 
                                          integrate(
                                            function(y){ 
                                              (1 / (1 + y ^ 10)) * dgamma(y, ctr_summary$ctr_obs[x] + 2, ctr_summary$ctr_exp[x] + 2)
                                            },
                                            lower = 0,
                                            upper = Inf)$value
                                          })
  
  
  ### Calculating the observed 5 tier
  ctr_summary <- ctr_summary %>%
    mutate(obs_5tier = ifelse(obs_5tier_score < 0.125, 1,
                              ifelse({obs_5tier_score >= 0.125} & {obs_5tier_score < 0.375}, 2,
                                     ifelse({obs_5tier_score >= 0.375} & {obs_5tier_score < 0.625}, 3,
                                            ifelse({obs_5tier_score >= 0.625} & {obs_5tier_score < 0.875}, 4,
                                                   ifelse(obs_5tier_score >= 0.875, 5,
                                                          NA))))))

  
  
  ctr_eval_0_3  <- ctr_summary$ctr[{ctr_summary$ctr_exp >= 0} & {ctr_summary$ctr_exp < 3}]
  ctr_eval_3_10 <- ctr_summary$ctr[{ctr_summary$ctr_exp >= 3} & {ctr_summary$ctr_exp < 10}]
  ctr_eval_10_  <- ctr_summary$ctr[{ctr_summary$ctr_exp >= 10}]
  
  
  
  op_char_mar <- c(mean(ctr_summary$ctr_exp),
                   mean(ctr_summary$true_hr),
                   mean((ctr_summary$HR_unshr - ctr_summary$true_hr)),
                   mean((ctr_summary$HR_shr - ctr_summary$true_hr)),
                   mean((ctr_summary$HR_unshr - ctr_summary$true_hr) ^ 2),
                   mean((ctr_summary$HR_shr - ctr_summary$true_hr) ^ 2),
                   cor(ctr_summary$HR_unshr, ctr_summary$true_hr, method = "kendall"),
                   cor(ctr_summary$HR_shr, ctr_summary$true_hr, method = "kendall"),
                   cor(ctr_summary$obs_5tier_score, ctr_summary$true_hr, method = "kendall"),
                   cor(ctr_summary$obs_5tier, ctr_summary$true_hr, method = "kendall"),
                   cor(ctr_summary$obs_cms_flag, ctr_summary$true_hr, method = "kendall"),
                   cor(ctr_summary$obs_mpsc_flag, ctr_summary$true_hr, method = "kendall"),
                   
                   mean(ctr_summary$true_hr[ctr_summary$obs_5tier == 1]),
                   mean(ctr_summary$true_hr[ctr_summary$obs_5tier == 2]),
                   mean(ctr_summary$true_hr[ctr_summary$obs_5tier == 3]),
                   mean(ctr_summary$true_hr[ctr_summary$obs_5tier == 4]),
                   mean(ctr_summary$true_hr[ctr_summary$obs_5tier == 5]),
                   
                   mean(ctr_summary$true_hr[ctr_summary$obs_cms_flag == 1] > 1.25),
                   mean(ctr_summary$true_hr[ctr_summary$obs_mpsc_flag == 1] > 1.25))
  
  
  op_char_0_3 <- c(mean(ctr_summary$ctr_exp[ctr_summary$ctr %in% ctr_eval_0_3]),
                   mean(ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_0_3]),
                   mean((ctr_summary$HR_unshr[ctr_summary$ctr %in% ctr_eval_0_3] - 
                              ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_0_3])),
                   mean((ctr_summary$HR_shr[ctr_summary$ctr %in% ctr_eval_0_3] - 
                              ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_0_3])),
                   mean((ctr_summary$HR_unshr[ctr_summary$ctr %in% ctr_eval_0_3] - 
                           ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_0_3]) ^ 2),
                   mean((ctr_summary$HR_shr[ctr_summary$ctr %in% ctr_eval_0_3] - 
                           ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_0_3]) ^ 2),
                   cor(ctr_summary$HR_unshr[ctr_summary$ctr %in% ctr_eval_0_3], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_0_3], method = "kendall"),
                   cor(ctr_summary$HR_shr[ctr_summary$ctr %in% ctr_eval_0_3], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_0_3], method = "kendall"),
                   cor(ctr_summary$obs_5tier_score[ctr_summary$ctr %in% ctr_eval_0_3], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_0_3], method = "kendall"),
                   cor(ctr_summary$obs_5tier[ctr_summary$ctr %in% ctr_eval_0_3], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_0_3], method = "kendall"),
                   cor(ctr_summary$obs_cms_flag[ctr_summary$ctr %in% ctr_eval_0_3], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_0_3], method = "kendall"),
                   cor(ctr_summary$obs_mpsc_flag[ctr_summary$ctr %in% ctr_eval_0_3], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_0_3], method = "kendall"),
                   
                   mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 1} & {ctr_summary$ctr %in% ctr_eval_0_3}]),
                   mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 2} & {ctr_summary$ctr %in% ctr_eval_0_3}]),
                   mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 3} & {ctr_summary$ctr %in% ctr_eval_0_3}]),
                   mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 4} & {ctr_summary$ctr %in% ctr_eval_0_3}]),
                   mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 5} & {ctr_summary$ctr %in% ctr_eval_0_3}]),
                   
                   mean(ctr_summary$true_hr[{ctr_summary$obs_cms_flag == 1} & {ctr_summary$ctr %in% ctr_eval_0_3}] > 1.25),
                   mean(ctr_summary$true_hr[{ctr_summary$obs_mpsc_flag == 1} & {ctr_summary$ctr %in% ctr_eval_0_3}] > 1.25))
  
  
  op_char_3_10 <- c(mean(ctr_summary$ctr_exp[ctr_summary$ctr %in% ctr_eval_3_10]),
                    mean(ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_3_10]),
                    mean((ctr_summary$HR_unshr[ctr_summary$ctr %in% ctr_eval_3_10] - 
                            ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_3_10])),
                    mean((ctr_summary$HR_shr[ctr_summary$ctr %in% ctr_eval_3_10] - 
                            ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_3_10])),
                    mean((ctr_summary$HR_unshr[ctr_summary$ctr %in% ctr_eval_3_10] - 
                            ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_3_10]) ^ 2),
                    mean((ctr_summary$HR_shr[ctr_summary$ctr %in% ctr_eval_3_10] - 
                            ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_3_10]) ^ 2),
                    cor(ctr_summary$HR_unshr[ctr_summary$ctr %in% ctr_eval_3_10], 
                        ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_3_10], method = "kendall"),
                    cor(ctr_summary$HR_shr[ctr_summary$ctr %in% ctr_eval_3_10], 
                        ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_3_10], method = "kendall"),
                    cor(ctr_summary$obs_5tier_score[ctr_summary$ctr %in% ctr_eval_3_10], 
                        ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_3_10], method = "kendall"),
                    cor(ctr_summary$obs_5tier[ctr_summary$ctr %in% ctr_eval_3_10], 
                        ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_3_10], method = "kendall"),
                    cor(ctr_summary$obs_cms_flag[ctr_summary$ctr %in% ctr_eval_3_10], 
                        ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_3_10], method = "kendall"),
                    cor(ctr_summary$obs_mpsc_flag[ctr_summary$ctr %in% ctr_eval_3_10], 
                        ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_3_10], method = "kendall"),
                    
                    mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 1} & {ctr_summary$ctr %in% ctr_eval_3_10}]),
                    mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 2} & {ctr_summary$ctr %in% ctr_eval_3_10}]),
                    mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 3} & {ctr_summary$ctr %in% ctr_eval_3_10}]),
                    mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 4} & {ctr_summary$ctr %in% ctr_eval_3_10}]),
                    mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 5} & {ctr_summary$ctr %in% ctr_eval_3_10}]),
                    
                    mean(ctr_summary$true_hr[{ctr_summary$obs_cms_flag == 1} & {ctr_summary$ctr %in% ctr_eval_3_10}] > 1.25),
                    mean(ctr_summary$true_hr[{ctr_summary$obs_mpsc_flag == 1} & {ctr_summary$ctr %in% ctr_eval_3_10}] > 1.25))
  
  
  op_char_10_ <- c(mean(ctr_summary$ctr_exp[ctr_summary$ctr %in% ctr_eval_10_]),
                   mean(ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_10_]),
                   mean((ctr_summary$HR_unshr[ctr_summary$ctr %in% ctr_eval_10_] - 
                           ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_10_])),
                   mean((ctr_summary$HR_shr[ctr_summary$ctr %in% ctr_eval_10_] - 
                           ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_10_])),
                   mean((ctr_summary$HR_unshr[ctr_summary$ctr %in% ctr_eval_10_] - 
                           ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_10_]) ^ 2),
                   mean((ctr_summary$HR_shr[ctr_summary$ctr %in% ctr_eval_10_] - 
                           ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_10_]) ^ 2),
                   cor(ctr_summary$HR_unshr[ctr_summary$ctr %in% ctr_eval_10_], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_10_], method = "kendall"),
                   cor(ctr_summary$HR_shr[ctr_summary$ctr %in% ctr_eval_10_], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_10_], method = "kendall"),
                   cor(ctr_summary$obs_5tier_score[ctr_summary$ctr %in% ctr_eval_10_], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_10_], method = "kendall"),
                   cor(ctr_summary$obs_5tier[ctr_summary$ctr %in% ctr_eval_10_], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_10_], method = "kendall"),
                   cor(ctr_summary$obs_cms_flag[ctr_summary$ctr %in% ctr_eval_10_], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_10_], method = "kendall"),
                   cor(ctr_summary$obs_mpsc_flag[ctr_summary$ctr %in% ctr_eval_10_], 
                       ctr_summary$true_hr[ctr_summary$ctr %in% ctr_eval_10_], method = "kendall"),
                   
                   mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 1} & {ctr_summary$ctr %in% ctr_eval_10_}]),
                   mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 2} & {ctr_summary$ctr %in% ctr_eval_10_}]),
                   mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 3} & {ctr_summary$ctr %in% ctr_eval_10_}]),
                   mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 4} & {ctr_summary$ctr %in% ctr_eval_10_}]),
                   mean(ctr_summary$true_hr[{ctr_summary$obs_5tier == 5} & {ctr_summary$ctr %in% ctr_eval_10_}]),
                   
                   mean(ctr_summary$true_hr[{ctr_summary$obs_cms_flag == 1} & {ctr_summary$ctr %in% ctr_eval_10_}] > 1.25),
                   mean(ctr_summary$true_hr[{ctr_summary$obs_mpsc_flag == 1} & {ctr_summary$ctr %in% ctr_eval_10_}] > 1.25))
  
  
  
  fnl <- c(survConcordance(Surv(observed_time, observed_event) ~ predict(curr_model),
                           data = curr_data)$concordance,
           sum(ctr_summary$ctr_exp),
           curr_km$surv[length(curr_km$surv)],
           exp(-c_effect_est_params$baseline_hazard$cumulative_hazard[366] * 
                 exp(quantile(curr_data$observed_lp, probs = c(0.75, 0.90, 0.95, 0.99)) + bl_scaler$root)),
           op_char_mar, 
           op_char_0_3, 
           op_char_3_10, 
           op_char_10_)
  
  
  return(fnl)
}







