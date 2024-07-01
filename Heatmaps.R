Heatmaps <- function(MyFoodWeb, NamesAndIndices, Abbrev=FALSE){
  #' @title Effects heatmaps
  #' @description
    #' Heatmaps of direct and net effects
  #' @param MyFoodWeb is the interaction square matrix of the food web of interest
  #' @param NamesAndIndices is the dataframe with the names, indices, and abbreviation of all species
  #' @param Abbrev is a boolean, TRUE if you want abbreviated names, else FALSE
  #' @returns is the two heatmaps: direct effects and net effects

  if (Abbrev==TRUE){
    rownames(MyFoodWeb) <- NamesAndIndices$Abbreviation
    colnames(MyFoodWeb) <- NamesAndIndices$Abbreviation
  }
  
  I <- diag(1, nrow = nrow(MyFoodWeb))
  D <- diag(MyFoodWeb) # extract the self-regulation (diagonal)
  MyFoodWeb <- MyFoodWeb / -D # normalize by self-regulation --> non-dimensional interaction matrix
  DirectEffects <- MyFoodWeb + I # Extract the normalized self-regulation
  NetEffects <- solve(I - DirectEffects)
  
  x <- rownames(MyFoodWeb)
  y <- colnames(MyFoodWeb)
  DirectEffects <- flipdim(DirectEffects, 1)
  x <- rev(x)
  DataForHeatmap <- cbind(expand.grid(X = x, Y = y), melt(unlist(DirectEffects)))
  DirectHeatmap <- ggplot(DataForHeatmap, aes(Y, X, fill = value)) + geom_tile() + scale_fill_gradient2(na.value = "black") +
    theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size = 7)) +
    labs(x = "Species", y = "Species", title = "Direct Effects")
  
  x <- rownames(MyFoodWeb)
  y <- colnames(MyFoodWeb)
  NetEffects <- flipdim(NetEffects, 1)
  x <- rev(x)
  DataForHeatmap <- cbind(expand.grid(X = x, Y = y), melt(unlist(NetEffects)))
  NetHeatmap <- ggplot(DataForHeatmap, aes(Y, X, fill = value)) + geom_tile() + scale_fill_gradient2(na.value = "black") +
    theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size = 7)) +
    labs(x = "Species", y = "Species", title = "Net Effects")
  
  out <- list(DirectHeatmap, NetHeatmap)
  return(out)
}













