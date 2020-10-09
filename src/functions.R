library(dplyr)


par_time_series_catch <- function(original, freq = "month", steps = 3, dlmPoly = 2, dlmSeas = 12,  num.cores = 2, error = "mape", xreg = NULL, a.a.args = list(NULL),
                                  ets.args = list(NULL), tbats.args = list(NULL),  n.n.args = list(NULL)){
  auto_args <- list(max.p = 5, max.q = 5, max.P = 2,
                    max.Q = 2, max.order = 5, max.d = 2, max.D = 1, start.p = 2,
                    start.q = 2, start.P = 1, start.Q = 1, stationary = FALSE,
                    seasonal = TRUE, nmodels = 94)
  ets_args <- list(model = "ZZZ", damped = NULL, alpha = NULL, beta = NULL,
                   gamma = NULL, phi = NULL, additive.only = FALSE, lambda = NULL,
                   biasadj = FALSE)
  tbats_args <- list(use.box.cox = NULL, use.trend = NULL,
                     use.damped.trend = NULL, seasonal.periods = NULL,
                     use.arma.errors = TRUE)
  n_n_args <- list(P = 1, repeats = 20, xreg = NULL, lambda = NULL,
                   model = NULL, subset = NULL, scale.inputs = TRUE)
  
  ## modifying lists to update code
  auto_args <- utils::modifyList(auto_args, a.a.args)
  ets_args <- utils::modifyList(ets_args, ets.args)
  tbats_args <- utils::modifyList(tbats_args, tbats.args)
  n_n_args <- utils::modifyList(n_n_args, n.n.args)
  startDate = zoo::as.Date(format(lubridate::date_decimal(utils::head(zoo::index(original[[1]]),1)), "%Y-%m-%d"))
  
  if (class(original)!="list"){
    print("Put the data into a list and retry dummy")
    model_list <- NULL
    stop()
  } else {
    ## forced 3 to account for MoQ movement to be recognized
    x <- lapply(original, utils::head, -steps)
    x <- lapply(x, stats::as.ts)
    non_ts <- unlist(lapply(x, stats::is.ts))
    if (FALSE %in% non_ts){
      print("Your list contains non time-series, please change them all to time-series")
      stop()
    }
    ##modelling functions
    cl <- parallel::makeCluster(getOption("cl.cores", num.cores))
    par_auto <- parallel::parLapply(cl,x,forecast::auto.arima, max.p = auto_args$max.p, max.q = auto_args$max.q, max.P = auto_args$max.P,
                                    max.Q = auto_args$max.Q, max.order = auto_args$max.order, max.d = auto_args$max.d, max.D = auto_args$max.D,
                                    start.p = auto_args$start.p, start.q = auto_args$start.q, start.P = auto_args$start.P, start.Q = auto_args$start.Q,
                                    stationary = auto_args$stationary, seasonal = auto_args$seasonal, nmodels = auto_args$nmodels)
    parallel::stopCluster(cl)
    cl <- parallel::makeCluster(getOption("cl.cores", num.cores))
    par_ets <- parallel::parLapply(cl, x, forecast::ets, model = ets_args$model, damped = ets_args$damped, alpha = ets_args$alpha, beta = ets_args$beta,
                                   gamma = ets_args$gamma, phi = ets_args$phi, additive.only = ets_args$additive.only, lambda = ets_args$lambda,
                                   biasadj = ets_args$biasadj)
    parallel::stopCluster(cl)
    cl <- parallel::makeCluster(getOption("cl.cores", num.cores))
    par_tbats <- parallel::parLapply(cl, x, forecast::tbats, use.box.cox = tbats_args$use.box.cox, use.trend = tbats_args$use.trend, use.damped.trend = tbats_args$use.damped.trend,
                                     seasonal.periods = tbats_args$seasonal.periods, use.arma.errors = tbats_args$use.arma.errors)
    parallel::stopCluster(cl)
    cl <- parallel::makeCluster(getOption("cl.cores", num.cores))
    par_hybrid <- parallel::parLapply(cl, x, forecastHybrid::hybridModel,  a.args = auto_args, e.args = ets_args, t.args = tbats_args, n.args = n_n_args)
    parallel::stopCluster(cl)
    cl <- parallel::makeCluster(getOption("cl.cores", num.cores))
    par_dlm <- parallel::parLapply(cl, x, dlmPara, dlmPoly = dlmPoly, dlmSeas = dlmSeas, steps = steps, freq = freq)
    parallel::stopCluster(cl)
    model_list <- c(par_auto, par_ets, par_tbats, par_hybrid, par_dlm)
  }
  print("Modelling Done")
  model_forecast <- extract_model_fit_forecast(model_list, steps = steps, xreg = xreg, freq = freq, num.cores = num.cores)
  print("Forecasting Done")
  best_model <- best_forecast_type(model_forecast, steps = steps, origin = original)
  print("Picking Best Model")
  if (error=="mape"){
    error_df <- best_model[[2]]
  } else if (error=="smape"){
    error_df <- best_model[[1]]
  } else if (error=="mase"){
    error_df <- best_model[[3]]
  } else {
    print("Not a valid error metric")
    stop()
  }
  
  models <- colnames(error_df)
  orig_x <- lapply(original, stats::window, start = lubridate::decimal_date(startDate))
  cl <- parallel::makeCluster(getOption("cl.cores", num.cores))
  models_final <- parallel::parLapply(cl, orig_x, forecast_function, models, freq, steps, dlmPoly, dlmSeas, num.cores, error, xreg, a.a.args,
                                      ets.args, tbats.args, n.n.args)
  parallel::stopCluster(cl)
  print("Forecasting with Best Model")
  best_model_forecast <- dplyr::bind_cols(models_final)
  colnames(best_model_forecast) <- models
  return(best_model_forecast)
}

forecast_function <- function(orig_x, models, freq = freq, steps = steps, dlmPoly = dlmPoly, dlmSeas = dlmSeas,  num.cores = num.cores, error = error, xreg = xreg, a.a.args = a.a.args,
                              ets.args = ets.args, tbats.args = tbats.args, n.n.args = n.n.args) {
  auto_args <- list(max.p = 5, max.q = 5, max.P = 2,
                    max.Q = 2, max.order = 5, max.d = 2, max.D = 1, start.p = 2,
                    start.q = 2, start.P = 1, start.Q = 1, stationary = FALSE,
                    seasonal = TRUE, nmodels = 94)
  ets_args <- list(model = "ZZZ", damped = NULL, alpha = NULL, beta = NULL,
                   gamma = NULL, phi = NULL, additive.only = FALSE, lambda = NULL,
                   biasadj = FALSE)
  tbats_args <- list(use.box.cox = NULL, use.trend = NULL,
                     use.damped.trend = NULL, seasonal.periods = NULL,
                     use.arma.errors = TRUE)
  n_n_args <- list(P = 1, repeats = 20, xreg = NULL, lambda = NULL,
                   model = NULL, subset = NULL, scale.inputs = TRUE)
  auto_args <- utils::modifyList(auto_args, a.a.args)
  ets_args <- utils::modifyList(ets_args, ets.args)
  tbats_args <- utils::modifyList(tbats_args, tbats.args)
  n_n_args <- utils::modifyList(n_n_args, n.n.args)
  i <- parent.frame()$i[]
  
  models <- models[i]
  
  
  if ( grepl("Auto Arima",models)){
    model <- forecast::auto.arima(orig_x, max.p = auto_args$max.p, max.q = auto_args$max.q, max.P = auto_args$max.P,
                                  max.Q = auto_args$max.Q, max.order = auto_args$max.order, max.d = auto_args$max.d, max.D = auto_args$max.D,
                                  start.p = auto_args$start.p, start.q = auto_args$start.q, start.P = auto_args$start.P, start.Q = auto_args$start.Q,
                                  stationary = auto_args$stationary, seasonal = auto_args$seasonal, nmodels = auto_args$nmodels)
    forecast_model <- forecast::forecast(model, h=steps, xreg=xreg)
    forecast_model <- c(forecast_model$x, forecast_model$mean)
  } else if ( grepl("ETS",models)){
    model <- forecast::ets(orig_x, model = ets_args$model, damped = ets_args$damped, alpha = ets_args$alpha, beta = ets_args$beta,
                           gamma = ets_args$gamma, phi = ets_args$phi, additive.only = ets_args$additive.only, lambda = ets_args$lambda,
                           biasadj = ets_args$biasadj)
    forecast_model <- forecast::forecast(model, h=steps, xreg=xreg)
    forecast_model <- c(forecast_model$x, forecast_model$mean)
  } else if (grepl("TBATS",models)){
    model <- forecast::tbats(orig_x, use.box.cox = tbats_args$use.box.cox, use.trend = tbats_args$use.trend, use.damped.trend = tbats_args$use.damped.trend,
                             seasonal.periods = tbats_args$seasonal.periods, use.arma.errors = tbats_args$use.arma.errors)
    forecast_model <- forecast::forecast(model, h=steps, xreg=xreg)
    forecast_model <- c(forecast_model$x, forecast_model$mean)
  } else if ( grepl("Hybrid",models)){
    model <- forecastHybrid::hybridModel(orig_x ,  a.args = auto_args, e.args = ets_args, t.args = tbats_args, n.args = n_n_args)
    forecast_model <- forecast::forecast(model, h=steps, xreg=xreg)
    forecast_model <- c(forecast_model$x, forecast_model$mean)
  } else {
    dlmPara <- function ( x, dlmPoly = dlmPoly, dlmSeas = dlmSeas, steps = steps , freq = freq){
      model.build <- function(p) {
        return(
          dlm::dlmModPoly(dlmPoly, dV=p[1], dW=p[2:3]) +
            dlm::dlmModSeas(dlmSeas, dV=p[4])
        )
      }
      model.mle <- dlm::dlmMLE(x, parm=c(0.1, 0, 1, 1), build=model.build)
      model.fit <- model.build(model.mle$par)
      if (freq == "month"){
        time <- 12
      } else if (freq == "day") {
        time <- 365
      } else if (freq == "week") {
        time <- 52
      } else if (freq == "quarter") {
        time <- 4
      } else { time <- 1}
      model.filtered <- dlm::dlmFilter(x, model.fit)
      model.smoothed <- dlm::dlmSmooth(x, model.fit)
      model.forecast <- dlm::dlmForecast(model.filtered, nAhead=steps)
      xf <- seq(max(zoo::index(x)), max(zoo::index(x))+steps/time, 1/time)
      xf <- xf[(-1)]
      aa <- model.forecast$a[,-1]*(-1)
      aa <- cbind(model.forecast$a[,1], aa)
      a <- drop(model.forecast$a%*%t(dlm::FF(model.fit)))
      a <- c(x,a)
      class(a) <- "dlm"
      return(a)
    }
    forecast_model <- dlmPara(orig_x,  dlmPoly = dlmPoly, dlmSeas = dlmSeas, steps = steps, freq = freq)
    
  }
  return(forecast_model)
}

best_forecast_type <- function(model_forecast, steps = steps, origin = NULL){
  ## use a list of models
  forecasts <- length(model_forecast)/5
  best_mape <- list()
  best_smape <- list()
  best_mase <- list()
  for (q in 1:forecasts){
    seq_list <- seq(q, length(model_forecast), by = forecasts)
    first_ts <- cbind(dplyr::bind_cols(model_forecast[seq_list]), as.numeric(origin[[q]]))
    colnames(first_ts) <- c("Auto Arima", "ETS", "TBATS", "Hybrid", "DLM", "Original")
    mape <- NA
    smape <- NA
    mase <- NA
    for (i in 1:5){
      mape[i] <- Metrics::mape(utils::tail(as.numeric(first_ts$Original),steps), utils::tail(as.numeric(first_ts[,i]),steps))
      smape[i] <- Metrics::smape(utils::tail(as.numeric(first_ts$Original),steps), utils::tail(as.numeric(first_ts[,i]),steps))
      mase[i] <- Metrics::mase(utils::tail(as.numeric(first_ts$Original),steps), utils::tail(as.numeric(first_ts[,i]),steps), step_size = steps)
    }
    
    minmape <- match(min(mape), mape)
    best_forecastmape <- as.data.frame(as.double(first_ts[,minmape]))
    colnames(best_forecastmape) <- paste(c(names(first_ts)[minmape]), q, sep =", " )
    best_mape[[q]] <- best_forecastmape
    
    minsmape <- match(min(smape), smape)
    best_forecastsmape <- as.data.frame(as.double(first_ts[,minsmape]))
    colnames(best_forecastsmape) <-  paste(c(names(first_ts)[minsmape]), q, sep =", " )
    best_smape[[q]] <- best_forecastsmape
    
    minmase <- match(min(mase), mase)
    best_forecastmase <- as.data.frame(as.double(first_ts[,minmase]))
    colnames(best_forecastmase) <-  paste(c(names(first_ts)[minmase]), q, sep =", " )
    best_mase[[q]] <- best_forecastmase
    
  }
  
  smape_df <- dplyr::bind_cols(best_smape)
  mape_df <- dplyr::bind_cols(best_mape)
  mase_df <- dplyr::bind_cols(best_mase)
  
  error_list <- list("smape" = smape_df, "mape" = mape_df, "mase" = mase_df)
  return(error_list)
}

extract_model_fit_forecast <- function(x, steps = 3, xreg = NULL, freq = "month", num.cores = 2){
  ## forecasts are put into lists based on forecasting type
  aa_ets_tbats = list()
  dlm = list()
  hybrid = list()
  for (i in 1:(length(x)*(3/5))){
    aa_ets_tbats[[i]] <- x[[i]]
  }
  for (i in (length(x)*(3/5)+1):(length(x)*(4/5))){
    hybrid[[i]] <- x[[i]]
  }
  hybrid <- hybrid[lengths(hybrid) != 0]
  for (i in (length(x)*(4/5)+1):length(x)){
    dlm[[i]] <- x[[i]]
  }
  dlm <- dlm[lengths(dlm) != 0]
  cl <- parallel::makeCluster(getOption("cl.cores", num.cores))
  final_forecast <- parallel::parLapply(cl, aa_ets_tbats, forecast::forecast, h=steps, xreg=xreg)
  parallel::stopCluster(cl)
  final_hybrid <- lapply(hybrid, forecast::forecast, h=steps, xreg=xreg, PI=FALSE)
  final_dlm <- dlm
  forecast_combine <- c(final_forecast, final_hybrid)
  final_reg <- lapply(forecast_combine, function(x){c(x$x, x$mean)})
  final <- c(final_reg, final_dlm)
  return(final)
}
dlmPara <- function ( x, dlmPoly = dlmPoly, dlmSeas = dlmSeas, steps = steps , freq = freq){
  model.build <- function(p) {
    return(
      dlm::dlmModPoly(dlmPoly, dV=p[1], dW=p[2:3]) +
        dlm::dlmModSeas(dlmSeas, dV=p[4])
    )
  }
  model.mle <- dlm::dlmMLE(x, parm=c(0.1, 0, 1, 1), build=model.build)
  model.fit <- model.build(model.mle$par)
  if (freq == "month"){
    time <- 12
  } else if (freq == "day") {
    time <- 365
  } else if (freq == "week") {
    time <- 52
  } else if (freq == "quarter") {
    time <- 4
  } else { time <- 1}
  model.filtered <- dlm::dlmFilter(x, model.fit)
  model.smoothed <- dlm::dlmSmooth(x, model.fit)
  model.forecast <- dlm::dlmForecast(model.filtered, nAhead=steps)
  xf <- seq(max(zoo::index(x)), max(zoo::index(x))+steps/time, 1/time)
  xf <- xf[(-1)]
  aa <- model.forecast$a[,-1]*(-1)
  aa <- cbind(model.forecast$a[,1], aa)
  a <- drop(model.forecast$a%*%t(dlm::FF(model.fit)))
  a <- c(x,a)
  class(a) <- "dlm"
  return(a)
}



fix_negatives <- function(final, colm_means){
  ## imputing missing data based on sampeling of the data
  for (i in 1:ncol(final)){
    for (q in 1:length(final[,i])){
      if (is.na(final[q,i])){
        final[q,i] <- colm_means[i]
      }
    }
  }
  return(final)
}
#removing weird values and replacing some odd values with means
remove_and_replace <- function(final){
  ##maing values NA
  final <- as.data.frame(lapply(final, function(x) {gsub("-Inf", NA, x)}), stringsAsFactors = F)
  final <- as.data.frame(lapply(final, function(x) {gsub("NaN", NA, x)}), stringsAsFactors = F)
  final <- as.data.frame(lapply(final, function(x) {gsub("Inf", NA, x)}), stringsAsFactors = F)
  
  
  ## making the columns numeric and transposing
  final <- sapply( final, as.numeric )
  finaldf <- as.data.frame(final)
  return(finaldf)
}



seperate_df_gt <- function(df, na_count){
  
  
  incomplete_cols<- df %>%
    select_if(~ any(is.na(.)))
  
  count_nas <- as.data.frame(t(colSums(is.na(incomplete_cols))))
  
  ##making df with missing values
  gt_missing <- incomplete_cols[,colnames(count_nas[,count_nas >=na_count])]
  return(gt_missing)
  
}

seperate_df_lt <- function(df, na_count){
  
  incomplete_cols<- df %>%
    select_if(~ any(is.na(.)))
  
  count_nas <- as.data.frame(t(colSums(is.na(incomplete_cols))))
  
  ##making df with missing values
  lt_missing <- incomplete_cols[,colnames(count_nas[,count_nas <na_count])]
  
  return(lt_missing)
}

seperate_df_gt_lt <- function(df, na_count_gt , na_count_lt){
  
  
  incomplete_cols<- df %>%
    select_if(~ any(is.na(.)))
  
  count_nas <- as.data.frame(t(colSums(is.na(incomplete_cols))))
  
  ##making df with missing values
  lt_gt_missing <- incomplete_cols[,colnames(count_nas[,count_nas >=na_count_gt & count_nas <na_count_lt])]
  return(lt_gt_missing)
  
}

replace_0_with_na <- function(df){
  # replacing 0 with NA
  df[df==0] <- NA
  return(df)
}


replace_gt1_lt0 <- function(df){
  # replacing negatives with NA
  df[df<=0] <- NA
  df[df>=1] <- NA
  return(df)
}

values_replaced <- function(new, old){
  if (dim(new)[1]==dim(old)[1] && dim(new)[2]==dim(old)[2]){
  } else {
    stop("Dimensions are different for the dataframes")
  }
  for (j in 1:NCOL(new)){
    for (i in 1:NROW(new)){
      if (is.na(new[i,j])|is.nan(new[i,j])|is.infinite(new[i,j])|is.null(new[i,j])){
        print(paste(colnames(new)[j], "at row",i, "is missing!" ))
      } else if (is.na(old[i,j])){
        print(paste(colnames(new)[j], "at row",i, "change by",(new[i,j]-0)))
      } else if (new[i,j]==old[i,j]){
      } else {
        print(paste(colnames(new)[j], "at row",i, "change by",(new[i,j]-old[i,j])))
      }
    }
  }
}


all_in_one_time_series <- function(original, freq = "month", steps = 3, dlmPoly = 2, dlmSeas = 12,  num.cores = 2, error = "mape", xreg = NULL, a.a.args = list(NULL),
                                   ets.args = list(NULL), tbats.args = list(NULL),  n.n.args = list(NULL)){
  
  if (class(original)!="data.frame"){
    stop("data isnt in data frame, please change it!")
  }
  
  original_imputed <- original %>% 
    mutate_all( .funs= ~ifelse(is.na(.), mean(., na.rm = TRUE), .))
  
  ## then looking for anomolies
  original_imputed_anom <- original_imputed %>% 
    mutate_all( .funs= ~ifelse(abs(.)>median(.)+3*sd(.), median(.), .))
  
  ## tell you what values where replaced
  values_replaced(original_imputed_anom, original)
  
  good_format <- lapply(original_imputed_anom, as.ts)
  
  output <- par_time_series_catch(good_format, freq = freq, steps = steps, dlmPoly = dlmPoly, dlmSeas = dlmSeas,  num.cores = num.cores, error = error, xreg = xreg, a.a.args = a.a.args,
                                  ets.args = ets.args, tbats.args = tbats.args,  n.n.args = n.n.args)
  ##collecting model info and pasteing it to the bottom
  models <- colnames(output)
  models <- gsub(",.*", "", models)
  colnames(output) <- colnames(original)
  output <- rbind.data.frame(output, models)
  
  #changing row name to models
  rownames(output)[nrow(output)] <- "Models"
  
  return(output)
  
}


mean_prediction <- function(df, steps){
  col_means <- colMeans(df, na.rm = TRUE)
  col_means <- ifelse(col_means=="NaN", 0, col_means)
  
  client_level_forecast_mean <- rbind.data.frame(df, t(replicate(steps,col_means)))
  client_level_forecast_mean <- rbind.data.frame(client_level_forecast_mean, replicate(ncol(client_level_forecast_mean), "Mean Predicition"))
  rownames(client_level_forecast_mean) <-  seq(1, nrow(client_level_forecast_mean))
  
  rownames(client_level_forecast_mean)[nrow(client_level_forecast_mean)] <- "Models"
  
  
  return(client_level_forecast_mean)
}


forecast_slack_trigger <- function(project_name){
  clue <- paste("The project", project_name, " has had its data appended to the database.")
  slackr::slackr_bot(clue,
                     channel = "ds-forecast",
                     username = "Forecasting Bot", 
                     icon_emoji = "",
                     incoming_webhook_url = Sys.getenv("SLACK_WEBHOOK"))
}
