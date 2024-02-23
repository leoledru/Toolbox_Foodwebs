FoodWebGraph <- function(MyFoodWeb, IdxFocusSpecies = NULL){
  
  # Trophic Levels for nodes' y position
  I <- diag(1, nrow = nrow(MyFoodWeb))
  DirectEffects <- MyFoodWeb + I # Extract the normalized self-regulation
  SecondOrderEffects <- DirectEffects%^%2
  Flow <- -1 * DirectEffects * (DirectEffects < 0)
  Flow[is.na(Flow)] <- 0
  TrophAndOmni <- TrophInd(Flow)
  TrophLevels <- round(TrophAndOmni$TL)
  
  # Create the graphic, emphasizes the selected chain, or species
  MatrixForGraph <- DirectEffects
  edges <- which(MatrixForGraph < 0, arr.ind = TRUE)
  MatrixForGraph <- abs(MatrixForGraph)
  MatrixForGraph <- log10(MatrixForGraph) # log10 to better visualization of different edges' weights
  EdgesAll <- as.data.frame(cbind(edges, abs(MatrixForGraph[edges]), rgb(0, 0, 0)))
  colnames(EdgesAll) <- c("prey","pred","value","color")
  
  if (!is.null(IdxFocusSpecies)){
    DirectLinks <- which(DirectEffects[, IdxFocusSpecies] != 0)
    SecondOrderLinks <- which(SecondOrderEffects[, IdxFocusSpecies] != 0)
    SecondOrderLinks <- SecondOrderLinks[!SecondOrderLinks == IdxFocusSpecies] # remove the focus species 
    # direct links in red
    EdgesAll$color <- ifelse(EdgesAll$prey == IdxFocusSpecies, rgb(1,0,0), EdgesAll$color)
    EdgesAll$color <- ifelse(EdgesAll$pred == IdxFocusSpecies, rgb(1,0,0), EdgesAll$color)
    # second-order links in blue
    EdgesAll$color <- ifelse(EdgesAll$prey %in% DirectLinks & EdgesAll$pred %in% SecondOrderLinks, rgb(0,0,1), EdgesAll$color)
    EdgesAll$color <- ifelse(EdgesAll$pred %in% DirectLinks & EdgesAll$prey %in% SecondOrderLinks, rgb(0,0,1), EdgesAll$color)
  }

  colnames(EdgesAll) <- c("from", "to","weight", "col")
  # x position des noeuds
  Nodes <- c(1:nrow(MyFoodWeb))
  x <- runif(nrow(MyFoodWeb))
  # y position des noeuds (defini par le niveau trophique)
  y <- TrophLevels
  Node_list <- data.frame(Nodes, x, y)
  
  # Créer le graphique à partir des arêtes et des positions x et y
  Graph <- graph_from_data_frame(vertices = Node_list, d = EdgesAll, directed = F)
  
  return(Graph)
}