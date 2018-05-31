
gen_surv_time <- function(r_unif,
                          bl_hazard,
                          bl_times){
  
  
  if(any(diff(bl_hazard) < 0)){ stop("Improper cumulative hazard function") }
  
  
  # If r_unif is greater than the cumulative density function at 365, then the observed time is 365 with an event indicator equal to 0
  if(r_unif >= (1 - exp(-bl_hazard[length(bl_hazard)]))){
    fnl <- c(365, 0)
  }
  
  # If r_unif is less than the cumulative density function at 365, then we have to find the appropriate day and set the event indicator to 1
  if(r_unif < (1 - exp(-bl_hazard[length(bl_hazard)]))){
    fnl <- c(bl_times[(1 - exp(-bl_hazard)) > r_unif][1],
             1)
  }
  
  return(fnl)
}

