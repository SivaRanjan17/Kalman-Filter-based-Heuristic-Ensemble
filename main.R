require (caret)
require (rpart)
require (PerfMeas)
library(adabag)

####             Training             #####

kfhe_train <- function (X, Y, max_iter, rpart_control = NA, blend_with_class = TRUE, early_stop = TRUE, wt_mode = "exp", reset_dist = TRUE)
{

  # Kalman filter variables for the model Kalman Filter, kf-m
  kalmanFilter1 <- list ()
  kalmanFilter1$pVariance <- NULL     # Posterior variance/error
  kalmanFilter1$mVariance <- NULL     # Measurement variance/error
  kalmanFilter1$gain <- NULL     # Kalman Gain
  kalmanFilter1$learntList <- list ()  # List of learnt models
  kalmanFilter1$intialModel <- NULL     # The initial model.
  kalmanFilter1$internalState <- NULL     # Internal state

  # Kalman filter variables for the kf-w
  kalmanFilter2 <- list ()
  kalmanFilter2$pDistVariance <- NULL     # Posterior variance/error
  kalmanFilter2$mDistVariance <- NULL     # Measurement variance/error
  kalmanFilter2$distGain <- list ()  # Kalman Gain
  kalmanFilter2$distribution <- NULL     # Distribution

  # Debugging information
  info <- list ()
  info$train_accuracy <- NULL     # Accuracy on the train blend error
  info$fScore <- NULL     # F-Score  on the train blend error
  info$state_err <- list ()  # Error function evaluation of the train blend with respect to ground truth
  info$blendError <- NULL     # Uniformly weighted training blend error
  info$distribution <- list ()  # Store distribution of all iterations
  info$reset <- rep (FALSE, max_iter)  # Trace in which iteration a distribution reset occurred
  info$blend_err_delta <- NULL

  # Final model to pack, stamp and return
  retmodel <- list ()

  if (is.na (rpart_control))
  {
    rpart_control <- rpart.control(minsplit = 20, cp = -1, maxdepth = 30)  # Default
  }

  # Initialise the kf-m statetrain_blend, this_err, this_blend_err
  kalmanFilter1$intialModel      <- rpart (formula = Y ~ ., data = cbind (X, Y), control = rpart_control)
  kalmanFilter1$internalState <- predict (kalmanFilter1$intialModel, X, type = "vector")

  unwt_comp_err    <- err_fun (kalmanFilter1$internalState, Y, NULL)  # Find the per datapoint error. This is a vector.

  # Initialise the state variance for kf-m
  kalmanFilter1$pVariance[1] <- 1.0  # No confidence

  # Initialise the state variance for kf-w
  kalmanFilter2$distribution <- rep (1 / nrow (X), nrow (X))
  kalmanFilter2$pDistVariance <- 1.0

  info$distribution[[1]] <- kalmanFilter2$distribution  #initial training distribution

  this_blend_err  <- sum (err_fun(kalmanFilter1$internalState, Y, NULL) * (1/nrow (X)))

  #### Start ensemble iterations ####
  for (t in 1:max_iter)
  {
    #### Model Kalman filter (kf-m) section ####
    proc_noise <- 0
    kalmanFilter1$internalState <- kalmanFilter1$internalState
    kalmanFilter1$pVariance[t]        <- kalmanFilter1$pVariance[t] + proc_noise

    #### Measurement update kf-m ####
    # Measurement process
    # Retry loop, if the error is more than a threshold then recompute
    dist_reset_flag <- FALSE
    while (1)
    {
      resample_flag <- TRUE
      resample_count <- 0
      while (resample_flag == TRUE)
      {
        bsamp <- sample (seq_len(nrow(X)), nrow (X), replace = TRUE, prob = as.vector(kalmanFilter2$distribution / sum (kalmanFilter2$distribution)))
        resample_flag <- sum (table(Y[bsamp]) == 0) != 0
        resample_count <- resample_count + 1
        if (resample_count >= 100)
        {
          resample_flag <- FALSE
        }
      }
      kalmanFilter1$learntList[[t]] <- rpart (formula = Y ~ ., data = cbind (X, Y)[bsamp,], control = rpart_control)

      this_pred   <- predict (kalmanFilter1$learntList[[t]], X, type = "prob")
      this_pred   <- get_fixed_matrix (this_pred, levels (Y))

      if (blend_with_class == TRUE)
      {
        this_cls    <- factor (levels (Y)[apply (this_pred, 1, which.max)], levels = levels (Y))
        dummy_pred  <- matrix (0, nrow = nrow (this_pred), ncol = ncol (this_pred))
        colnames (dummy_pred) <- colnames (Y)
        for (i in seq_len(nrow(dummy_pred)))
        {
          dummy_pred[i,this_cls[i]] <- 1
        }
        this_pred <- dummy_pred
      }

      this_measurement <- (this_pred + kalmanFilter1$internalState)/2
      # Get the error for "this_measurement"
      unwt_comp_err  <- err_fun (this_measurement, Y, NULL)  # Find the per datapoint error. This is a vector.
      # (paste ("unwt2:", length(unwt_comp_err)))
      uniwt_comp_err <- unwt_comp_err * (1/nrow (X))         # Get a uniformly weighted version. This is a weighted vector.
      this_pred_err       <- err_fun (this_pred, Y, NULL)
      uniwt_this_pred_err <- this_pred_err * (1/nrow (X))

      this_m_err     <- uniwt_comp_err
      this_d_err     <- uniwt_comp_err

      # Measurement error for the model kf-m. This is a heuristic and can be computed in different ways.
      kalmanFilter1$mVariance[t]          <- sum (this_m_err)

      if (reset_dist == TRUE)
      {
        # Retry with uniform distribution if the error is more than (1 - 1/nclass)
        # Reset only if it was not reset in the immediately previous attempt. Saves from infinite loops.
        if ((sum (uniwt_this_pred_err) >= (1 - 1/length (levels(Y)))) && (dist_reset_flag == FALSE))
        {
          kalmanFilter2$distribution <- rep (1/nrow (X), nrow (X))
          kalmanFilter2$pDistVariance <- 1.0

          dist_reset_flag <- TRUE
          info$reset[t] <- dist_reset_flag
          next
        }
        else
        {
          break
        }
      }
      else
      {
        break
      }
    }

    # Compute the Kalman gain for the kf-m
    kalmanFilter1$gain[t]        <- kalmanFilter1$pVariance[t] / (kalmanFilter1$pVariance[t] + kalmanFilter1$mVariance[t] + .Machine$double.xmin)
    # Update internal state for kf-m
    kalmanFilter1$internalState <- kalmanFilter1$internalState + kalmanFilter1$gain[t] * (this_measurement - kalmanFilter1$internalState)
    prev_blend_err   <- this_blend_err
    # Update internal error for the kf-m
    P_t_pred <- kalmanFilter1$pVariance[t] - kalmanFilter1$pVariance[t] * kalmanFilter1$gain[t]

    kalmanFilter1$pVariance[t+1] <- P_t_pred

    # Compute the actual error for the internal state
    info$state_err[[t]]   <- err_fun(kalmanFilter1$internalState, Y, NULL)
    this_blend_err              <- sum (err_fun(kalmanFilter1$internalState, Y, NULL) * (1/nrow (X)))
    info$blendError[t] <- this_blend_err

    train_pred_cls               <- factor(apply (kalmanFilter1$internalState, 1, which.max), levels = levels (Y))

    info$train_accuracy[t] <- confusion_metric (train_pred_cls, Y, "A")
    info$fScore[t]        <- confusion_metric (train_pred_cls, Y, "F")

    #############################################
    #### Weight Kalman filter (kf-w) section ####
    #############################################
    proc_noise <- 0
    kalmanFilter2$pDistVariance      <- kalmanFilter2$pDistVariance + proc_noise

    # Measurement of state vector of the distribution kf-w
    if (wt_mode == "exp")
    {
      dtemp           <- unwt_comp_err
      dtemp           <- dtemp + 1/nrow(X)
      kalmanFilter2$D_t_obs    <- kalmanFilter2$distribution * exp (dtemp)
    }

    if (wt_mode == "linear") {
      dtemp           <- unwt_comp_err
      kalmanFilter2$D_t_obs    <- rep (1/nrow(X), nrow(X))
      kalmanFilter2$D_t_obs    <- kalmanFilter2$D_t_obs + kalmanFilter2$distribution * (dtemp)
      kalmanFilter2$D_t_obs    <- kalmanFilter2$D_t_obs / sum (kalmanFilter2$D_t_obs)
    }

    # Measurement error mDistVariance for the distribution kf-w
    kalmanFilter2$mDistVariance        <- sum (this_d_err)
    # Compute the Kalman gain for the distribution kf-w
    kalmanFilter2$distGain[[t]]   <- kalmanFilter2$pDistVariance / (kalmanFilter2$pDistVariance + kalmanFilter2$mDistVariance + .Machine$double.xmin)
    # Update iternal state for the distribition kf-w
    kalmanFilter2$distribution             <- kalmanFilter2$distribution + kalmanFilter2$distGain[[t]] #* (kalmanFilter2$D_t_obs - kalmanFilter2$distribution)
    # Update internal error for the distribution KF
    kalmanFilter2$pDistVariance        <- kalmanFilter2$pDistVariance - kalmanFilter2$pDistVariance * kalmanFilter2$distGain[[t]]
    # Compute the change in the internal state actual error
    info$blend_err_delta[t] <- prev_blend_err - this_blend_err
    # Normalise distribution
    kalmanFilter2$distribution  <- kalmanFilter2$distribution / sum (kalmanFilter2$distribution)

    info$distribution[[t]]        <- kalmanFilter2$distribution

    if ((early_stop == TRUE) && (kalmanFilter1$gain[t] == 0))
    {
      max_iter <- t
      break
    }
  }

  # Pack
  retmodel$kalmanFilter1 <- kalmanFilter1
  retmodel$kalmanFilter2 <- kalmanFilter2
  retmodel$info <- info
  retmodel$max_iter <- max_iter
  retmodel$cls_lvls <- levels (Y)
  retmodel$blend_with_class <- blend_with_class

  print (paste("Accuracy: ", retmodel$info$train_accuracy))
  print (paste("Fscore: ",retmodel$info$fScore))

  class (retmodel)  <- "kfhe_m"
  return (retmodel)
}

####         Utility functions         ####
get_fixed_matrix <- function (mat_to_fix, target_levels)
{
  template              <- as.data.frame (matrix (0, nrow = nrow (mat_to_fix), ncol = length (target_levels)))
  colnames (template)   <- target_levels
  temp_names            <- colnames (mat_to_fix)
  template[,temp_names] <- mat_to_fix
  return (template)
}
err_fun <- function (pred, target, wt)
{
  if (!is.data.frame(target))
  {
    target         <- as.data.frame (target)
    one_hot_target <- model.matrix (~ ., data = target, contrasts.arg = lapply (target, contrasts, contrasts = FALSE))
    one_hot_target <- one_hot_target[,-1]
    target         <- one_hot_target
  }
  else if (ncol (target) > 1)
  {

  }
  if (is.null (wt) == TRUE)
  {
    wt <- 1
  }
  final_err <- wt * (sapply (pred, which.max) != sapply (target, which.max))

  return (final_err)
}
confusion_metric <- function (pred_vec, target_vec, metric)
{
  pred_one_hot_mat   <- to_one_hot (pred_vec)
  colnames (pred_one_hot_mat)   <- LETTERS[seq_along(levels(pred_vec))]
  target_one_hot_mat <- to_one_hot (target_vec)
  colnames (target_one_hot_mat) <- LETTERS[seq_along(levels(target_vec))]
  fobj <- F.measure.single.over.classes (target_one_hot_mat, pred_one_hot_mat)
  return (fobj$average[metric])
}
to_one_hot <- function (one_cool_vec)
{
  total_lvls <- length (levels(one_cool_vec))
  one_cool_vec <- as.numeric (one_cool_vec)
  one_hot_matrix <- matrix (0, nrow = length (one_cool_vec), ncol = total_lvls, byrow = TRUE)
  for (i in seq_along(one_cool_vec))
  {
    one_hot_matrix[i,one_cool_vec[i]] <- 1
  }
  return (one_hot_matrix)
}


filename <- "cmc.data"
dataset <- read.csv(filename, header=FALSE)
dataset$V10 <-as.factor(dataset$V10)

print ("CMC Dataset")
model <- kfhe_train(dataset[1:9],  dataset$V10, 10)


filename2 <- "iris.data"
dataset2 <- read.csv(filename2, header = FALSE)
dataset2$V5 <- as.factor(dataset2$V5)

print ("Iris Dataset")
model <- kfhe_train(dataset2[1:4],  dataset2$V5, 10)


indexes <- createDataPartition(iris$Species, p=.90, list = F)
train <- iris[indexes, ]
test <- iris[-indexes, ]

model <- boosting(Species~., data=train, boos=TRUE, mfinal=50)
pred <- predict(model, test)
print ("Adaboost Iris Dataset")
#print(pred$confusion)
print(pred$error)

cvmodel <- boosting.cv(V10~., data=dataset, boos=TRUE, mfinal=10, v=10)
#print (cvmodel$confusion)
print ("Adaboost cmc Dataset")
print (cvmodel$error)




