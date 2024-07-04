ChainCount <- function(MyFoodWeb, NamesAndIndices){
  #' @title Trophic chains count and length
  #' @description
    #' Identifies all the trophic chains, i.e. number of paths from apex predator(s) (predator species without predator)
    #' to basal species, stores them and measures their length. Gives also the mean length and the standard-deviation.
  #' @param MyFoodWeb is the interaction square matrix of the food web of interest
  #' @param NamesAndIndices is the dataframe with the names, indices, and abbreviation of all species
  #' @returns a list of lists giving all the trophic chains by the indices of species as well as the mean length of chains and the standard deviation
  
  # Trophic Levels
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
  xSample <- c(1:nrow(MyFoodWeb))
  x <- sample(xSample)
  y <- yCoord # real troph levels
  NodeList <- data.frame(Nodes, x, y, rgb(1,1,1,))
  colnames(NodeList) <- c("nodes","x","y","color")
  
  # Create the graph object
  MatrixForGraph <- DirectEffects
  edges <- which(MatrixForGraph < 0, arr.ind = TRUE)
  MatrixForGraph <- abs(MatrixForGraph)
  EdgesAll <- as.data.frame(cbind(edges, abs(MatrixForGraph[edges]), rgb(0, 0, 0)))
  colnames(EdgesAll) <- c("from", "to","weight", "col")
  colnames(NodeList) <- c("nodes","x","y","color")
  # Créer le graphique à partir des arêtes et des positions x et y
  Graph <- graph_from_data_frame(vertices = NodeList, d = EdgesAll, directed = T)
  
  Roots <- which(degree(Graph, v = V(Graph), mode = "in")==0, useNames = T) # all species without prey
  Leaves <- which(degree(Graph, v = V(Graph), mode = "out")==0, useNames = T) # all species without pred
  AllChains <- list()
  for (Root in Roots){
    AllChains <- append(AllChains, all_simple_paths(Graph, from = Root, to = Leaves))
  }
  
  ChainsLength <- lengths(AllChains)
  MeanChainsLength <- mean(ChainsLength)
  SdChainsLength <- sd(ChainsLength)
  
  Out <- list("AllChains" = AllChains, "MeanChainsLength" = MeanChainsLength, "SdChainsLength" = SdChainsLength)
  return(Out)
}