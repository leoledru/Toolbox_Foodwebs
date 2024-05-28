FoodWebGraph <- function(MyFoodWeb, IdxFocusSpecies = NULL){
  
  # Trophic Levels for nodes' y position
  I <- diag(1, nrow = nrow(MyFoodWeb))
  D <- diag(MyFoodWeb) # extract the self-regulation (diagonal)
  MyFoodWeb <- MyFoodWeb / -D # normalize by self-regulation --> non-dimensional interaction matrix
  DirectEffects <- MyFoodWeb + I # Extract the normalized self-regulation
  SecondOrderEffects <- DirectEffects%^%2
  Flow <- -1 * DirectEffects * (DirectEffects < 0)
  Flow[is.na(Flow)] <- 0
  TrophAndOmni <- TrophInd(Flow)
  TrophLevels <- round(TrophAndOmni$TL)
  yCoord <- TrophAndOmni$TL
  
  # Nodes coordinates
  Nodes <- c(1:nrow(MyFoodWeb))
  # x <- runif(nrow(MyFoodWeb))
  xSample <- c(1:nrow(MyFoodWeb))
  x <- sample(xSample)
  # y position des noeuds (defini par le niveau trophique)
  y <- yCoord # real troph levels
  # y <- TrophLevels # round troph levels
  NodeList <- data.frame(Nodes, x, y, rgb(1,1,1,))
  colnames(NodeList) <- c("nodes","x","y","color")
  
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
    # color the focus species
    NodeList$color <- ifelse(NodeList$nodes == IdxFocusSpecies, rgb(0,1,0,.5), NodeList$color)
    # direct links & nodes in red
    EdgesAll$color <- ifelse(EdgesAll$prey == IdxFocusSpecies, rgb(1,0,0), EdgesAll$color)
    EdgesAll$color <- ifelse(EdgesAll$pred == IdxFocusSpecies, rgb(1,0,0), EdgesAll$color)
    NodeList$color <- ifelse(NodeList$nodes %in% DirectLinks, rgb(1,0,0,.5), NodeList$color)
    # second-order links in blue
    EdgesAll$color <- ifelse(EdgesAll$prey %in% DirectLinks & EdgesAll$pred %in% SecondOrderLinks, rgb(0,0,1), EdgesAll$color)
    EdgesAll$color <- ifelse(EdgesAll$pred %in% DirectLinks & EdgesAll$prey %in% SecondOrderLinks, rgb(0,0,1), EdgesAll$color)
    NodeList$color <- ifelse(NodeList$nodes != DirectLinks & NodeList$nodes %in% SecondOrderLinks, rgb(0,0,1,.5), NodeList$color)
  }

  colnames(EdgesAll) <- c("from", "to","weight", "col")
  colnames(NodeList) <- c("nodes","x","y","color")
  
  # Créer le graphique à partir des arêtes et des positions x et y
  Graph <- graph_from_data_frame(vertices = NodeList, d = EdgesAll, directed = F)
  
  return(Graph)
}