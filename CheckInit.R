CheckInit <- function(MyFoodWeb){
  #' @title Check of food web interaction matrix
  #' @param MyFoodWeb is the interaction square matrix of the food web of interest
  #' @description
    #' Check if the interaction matrix representing the food web is correct and ready for further analysis
  #' @returns the checked interaction matrix, which can be different than in input because unconnected species are removed 
    
  if (!is.matrix(MyFoodWeb)){
    stop("Error: the variable is not a matrix")
  }
  if (!all(is.numeric(MyFoodWeb))){
    stop("Error: the matrix must contains numerics only")
  }
  if (nrow(MyFoodWeb) != ncol(MyFoodWeb)){
    stop("Error: the matrix is not square")
  }
  # No diagonal element can be empty, zero or positive
  for (element in diag(MyFoodWeb)) {
    if (is.na(element) || element == 0 || element > 0) {
      stop("Error: diagonal elements must be numeric < 0. Check the diagonal of the interaction matrix. If elements of the diagonal have not been inferred, set the value to -1.")
    }
  }
  # Verification of positive-negative symmetry (trophic links only)
  for (i in 1:nrow(MyFoodWeb)){
    for (j in 1:ncol(MyFoodWeb)){
      if (i != j){
        if ((MyFoodWeb[i, j] < 0 && MyFoodWeb[j, i] <= 0) || (MyFoodWeb[i, j] > 0 && MyFoodWeb[j, i] >= 0)){
          stop(paste("Error: Off-diagonal negative (positive) values must have a positive (negative) value at their symmetry. Problem with positions i =", i, "et j =", j))
        }
      }
    }
  }
  # Check for unconnected species: if they exist, remove them from the network
  FoodWebBis <- MyFoodWeb
  diag(FoodWebBis) <- 0 # remove self-reg to test if each species is connected with at least another one
  IdxUnconnected <- which(colSums(FoodWebBis) == 0)
  if (length(IdxUnconnected) > 0){
    MyFoodWeb <- MyFoodWeb[-IdxUnconnected,-IdxUnconnected]
    print(paste0("Species ", paste(IdxUnconnected, collapse = ", "), " are removed because unconnected"))
  }
  print("The food web is validated")
  return(MyFoodWeb)
}