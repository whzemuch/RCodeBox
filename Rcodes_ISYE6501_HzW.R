# calculate the r squared and/or adjusted r squared
cal_r_square_adjust <- function(y, y_pred, num_vars=NULL){
  # y: vector, actual value of response variable
  # y_pred: vector, predicted value of response variable
  # numbvars: number, the number of predictors for regression
  # Return: R-squared and adjusted-R_squared if num_vars provided
  n = length(y)
  RSS = sum((y - y_pred)^2) # Residual Sum of Squares
  TSS = sum((y - mean(y))^2) # Total Sum of Squares
  
  r_squared = 1 - RSS/TSS
  if (is.null(num_vars)) {
    return(r_squared)
  } else{
    k = num_vars
    return(list(r_squared=r_squared,
                adjust_r_squared = 1 - (1 - r_squared)*(n - 1)/(n - k - 1)))
  }
}

# scale all the columns except columns in col_skip
scale_train_df <- function(data, cols_skip) {
  # scale all the columns except the response variable Y
  # data: data frame include the response variable Y
  # Y: the column name of the response variable
  # Return: all scaled column except the response variable
  if (!all(cols_skip %in% names(data))) {
    stop("All the cols_skip must be  the columns of data!")
  }
  
  df_scaled_cols = dplyr::select(data, -any_of(cols_skip))
  df_no_scaled_cols = data[, cols_skip]
  
  return(cbind(scale(df_scaled_cols), df_no_scaled_cols))
}

# scale test data using the mean and sd from train data

# helper function for performing leave one out cross validation on a lm model
loo_cv_lm <- function(data, Y,  lm_formula) {
  # data: scaled data frame
  # Y: the column name of the responsive variable
  # lm_formula: the formula object for lm
  # Return: R-squared
  rows <- nrow(data)
  
  SSE_lst = sapply(1:rows, function(row){
    # Split data
    train <- data[-row, ]
    test <- data[row, ]
    
    lm_fit <- lm(formula=lm_formula, data = train)
    
    y_test <- test$Crime
    y_pred <- predict(lm_fit, newdata = test)
    
    SSE = (y_test - y_pred)^2
  })
  
  # calculated the R^2
  R_squared = 1 - sum(SSE_lst) / sum( (data[[Y]]-mean(data[[Y]]))^2)
  
  return(R_squared)
  
}

# scale test data using the mean and sd from train data
scale_test_df <- function(test_data, train_df_wo_y) {
  # scale test data using the mean and sd from train data
  # notes: do not include the response variable Y
  # test_data: data frame of test data
  # train_df_wo_y:  dataframe
  # Return the scaled data frame
  
  scaled_obj = scale(train_df_wo_y)
  mean_train = attr(scaled_obj, "scaled:center")
  sd_train = attr(scaled_obj, "scaled:scale")
  
  data_1 <- sweep(test_data, 2, STATS = mean_train, FUN = "-")
  data_1 <- sweep(data_1, 2, STATS = sd_train, FUN = "/")
  
  return (data_1)
  
}

# scale all the columns except the response variable Y
scale_train_df <- function(data, Y) {
  # scale all the columns except the response variable Y
  # data: data frame include the response variable Y
  # Y: the column name of the response variable
  # Return: all scaled column except the response variable
  if (!Y %in% names(data)) {
    stop("The Y must is not a column of data!")
  }
  
  df_wo_y = data[, -which(names(data) == "Crime")]
  df_y = data[Y]
  
  return(cbind(scale(df_wo_y), df_y))
}

# construct a new formula based on a lm object and p value
make_formula_obj <- function(lm_obj, pvalue_threshold=NULL){
  # construct a new formula by filtering some response varialbes based on p value
  # lm_obj: lm object
  # p_value_threshold: number
  # Return: a formula that can be used to train a new lm model
  if(is.null(pvalue_threshold)){
    return(formula(lm_obj))
  }else{
    terms_filter <- broom::tidy(lm_obj) %>% 
      filter(p.value < pvalue_threshold & term != "(Intercept)") %>%
      pull(term)
    
    formula_obj <- as.formula(paste("Crime ~ ", 
                                    paste0(terms_filter, collapse = " + ")))
    return(formula_obj)
  }
  
  
}


# split data into three parts:train, validation and test dataset
split.train_val_test <- function(data){
  n_sample = nrow(data)
  idx_train = sample(n_sample, size = floor(n_sample * 0.6))
  train = data[idx_train, ] # training data set 
  val_test = data[-idx_train, ]
  n_val = nrow(val_test)
  idx_val = sample(n_val, size = floor(n_val * 0.5))
  val = val_test[idx_val, ]
  test = val_test[-idx_val, ]
  
  return(list(train=train, val=val, test=test))
}

# Calculate the accuracy of a KNN model 
cal_kknn_accuracy <- function(model, data){
  y = data[, 11]
  pred = as.integer(fitted(model) + 0.5)
  return(sum(pred == y) / nrow(data))
}

# Calculate the accuracy of a SVM model
cal_ksvm_accuracy <- function(model, data) {
  X = data[, -11]
  y = data[, 11]
  pred <- predict(model, X)
  
  return(sum(pred == y ) / nrow(X))
}


# calculate the accuracy of a ksvm model
cal_model_accuracy <- function(fit, X, y) {
  pred <- predict(fit, X)
  accuracy <- sum(pred == y ) / nrow(X)
  return(accuracy)
}

# show the equation of a ksvm model
show_model_equation <- function(fit){
  
  a <- colSums(fit@xmatrix[[1]] * fit@coef[[1]])
  a0 <- -fit@b
  
  col_names <- colnames(data[, 1:10])
  equation = paste0(paste(a, col_names, collapse = " + ", sep = ""), " + ", a0, " = 0")
  
  return(equation)
}
