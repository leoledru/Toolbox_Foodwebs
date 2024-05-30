###################################################################################################################
#                    Reconstruct data set for future lake food webs                                               #
#                                                                                                                 #
# VAGNON Chloe                                                                                                    #
# 2022                                                                                                            #
# Univ. Savoie Mont Blanc                                                                                         #
#                                                                                                                 #
# This script references all the functions used for the article :                                                 #
# ""                                                                                                              #
#                                                                                                                 #
#                                                                                                                 #
# It is composed of 9 functions:                                                                                  #
#     - Niche_Attributes                                                                                          #
#     - Binary_Mat                                                                                                #
#     - Ref_L_Diet                                                                                                #
#     - Weighting                                                                                                 #
#     - Rel_Ab_Mat                                                                                                #
#     - Make_Bern                                                                                                 #
#     - Jacobian_Mat                                                                                              #
#     - Scale_JmatBS                                                                                              #
#     - Stability                                                                                                 #
#                                                                                                                 #
# Each function is detailed along the script (use the summary for quick access).                                  #
#                                                                                                                 #
###################################################################################################################





################################################## Niche_Attributes ################################################
# Infer the niche attributes of one or more consumers from Vagnon et al. 2021
# Necessitates to load the packages "stringr".
# Capital or lowercase letters work for the species category.
#
# Input:
# species_name = the name of the species to use
# body_size = log10(vector of consumer body sizes in ?m)
# species_category = "vertebrate" or "vertebre" or "invertebrate" or "invertebre" or other
#    /!\ If "other" is mentionned, the species will automatically be considered as producer and not as consumer
#
# Returns :
# a data frame with the niche attributes of a species

Niche_Attributes <- function(species_name, body_size, species_category) {

  library(stringr)

  Niche <- data.frame(name=NA, n=NA, c=NA, low=NA, high=NA)

  for(i in 1:length(body_size)){

    Niche2 <- data.frame(name=NA, n=NA, c=NA, low=NA, high=NA)
    
    # If not vertebrate or not invertebrate or not vertebre or not invertebre or not zoop
    if(!str_detect(species_category[i], regex("vertebrate", ignore_case = TRUE)) |
       !str_detect(species_category[i], regex("invertebrate", ignore_case = TRUE)) |
       !str_detect(species_category[i], regex("vertebre", ignore_case = TRUE)) |
       !str_detect(species_category[i], regex("invertebre", ignore_case = TRUE))|
       !str_detect(species_category[i], regex("zoop", ignore_case = TRUE))){
      # Estimate the parameters for the niche model
      Niche2$name = species_name[i]
      Niche2$n = body_size[i]	# The niche n
      Niche2$low = 0	# The lower limit of the range
      Niche2$high = 0	# The higher limit of the range
      Niche2$c = 0 # The centroid c
    }
  
    # If vertebrate and not inv or vertebre and not inv
    if(str_detect(species_category[i], regex("vertebrate", ignore_case = TRUE))&
       !str_detect(species_category[i], regex("inv", ignore_case = TRUE))|
       str_detect(species_category[i], regex("vertebre", ignore_case = TRUE))&
       !str_detect(species_category[i], regex("inv", ignore_case = TRUE))){
      # Unwrap the input parameters
      qrsup = Param_regvert[[2]]
      qrinf = Param_regvert[[3]]
      # Estimate the parameters for the niche model
      Niche2$name = species_name[i]
      Niche2$n = body_size[i] # The niche n
      Niche2$low = qrinf[1] + qrinf[2]*body_size[i]	# The lower limit of the range
      Niche2$high = qrsup[1] + qrsup[2]*body_size[i] # The higher limit of the range
      Niche2$c = Niche2$low + (Niche2$high - Niche2$low)/2 # The centroid c
    }

    # If invertebrate or invertebre or zoop
    if(str_detect(species_category[i], regex("invertebrate", ignore_case = TRUE))|
       str_detect(species_category[i], regex("invertebre", ignore_case = TRUE))|
       str_detect(species_category[i], regex("zoop", ignore_case = TRUE))){
      # Unwrap the input parameters
      qrsup = Param_reginvert[[2]]
      qrinf = Param_reginvert[[3]]
      # Estimate the parameters for the niche model
      Niche2$name=species_name[i]
      Niche2$n = body_size[i] # The niche n
      Niche2$low = qrinf[1] + qrinf[2]*body_size[i]	# The lower limit of the range
      Niche2$high = qrsup[1] + qrsup[2]*body_size[i]	# The higher limit of the range
      Niche2$c = Niche2$low+(Niche2$high-Niche2$low)/2 # The centroid c
    }
    Niche <- rbind(Niche, Niche2)
  }


  return(na.omit(Niche))
}


############################################### Binary_Mat ###########################################################
# Transform the parameters into a binary interaction matrix based on Gravel et al. 2013
# Modified to also obtain the data frame referencing trophic links corresponding to the binary matrix
# Input:
# name = the name of the species to use
# n = vector of size S with the parameter n for each of the S species
# c = vector of size S with the parameter c for each of the S species
# low = vector of size S with the parameter low for each of the S species
# high = vector of size S with the parameter high for each of the S species
# table = "YES" or "NO" if the matrix is also wanted in data.frame counting one line per interaction
#
#
# Returns :
# if table = "NO" a SxS matrix with 0 indicating absence of a link and 1 indicating the presence of a link
# Predators are on columns, preys are on rows
# if table = "YES" or a list with the SxS matrix + A table with the name of the prey, the predator and the body size of the prey and the predator

Binary_Mat = function(name, n, c, low, high, table) {

  S = length(n)
  L = matrix(0, nr=S, nc=S)
  colnames(L) <- name
  rownames(L) <- name

  for(i in 1:S){
    for(j in 1:S){
      # If the size of the prey j is in the range of the pred i
      if(n[j] > low[i] & n[j] < high[i]){
        L[j,i] <- 1
      }
    }
  }

  if(table == "NO"){
    return(L)
  }

  if(table == "YES"){
     Table <- data.frame(Prey=NA, Pred=NA, Log_Size_Prey=NA, Log_Size_Pred=NA)
    for(i in 1:S){
      if(length(which(L[,i] == 1))!=0){
      Table2 <- data.frame(Prey = names(which(L[,i] == 1)),
                         Pred = rep(colnames(L)[i], length(which(L[,i] == 1))),
                         Log_Size_Prey = n[which(L[,i] == 1)],
                         Log_Size_Pred = n[i])} else{
                           Table2 <- data.frame(Prey=NA, Pred=NA, Log_Size_Prey=NA, Log_Size_Pred=NA)}
    Table <- rbind(Table, Table2)}

    return(list(Bmat=L, Table=na.omit(Table)))
  }
}


###############################################  Ref_L_Diet  ##################################################
# Refine links for impossible trophic links according to species diet trait
# Fish do not eat primary producers
# Carnivorous macroinvertebrates do not eat primary producers
# Other diet refinement must be done after applying the function when they are considered as site-specific
#
# Inputs:
# Bmat = binary matrix of trophic links inferred from Binary_Mat function with :
#        colnames and rownames = names of species in the inventory
# diet = "invP" or "inv" for Invertebrates | "zoo" for zooplankton | "p" = piscivorous or "o" = omnivorous for fish | "prod" for primary producers
# Table = "YES" or "NO" to obtain or not the matrix with refined links
# LinksTab = if Table = "YES" provide the table of links obtained with the function L_fn2
#
# Returns :
# if table = "NO", returns th Binary matrix after link refinement
# if table = "YES", returns th Binary matrix after link refinement + the table of links after links refinement

Ref_L_Diet = function(Bmat, diet, Table, LinksTab) {
  if(Table == "NO"){
    for(i in 1:ncol(Bmat)){
      for(j in 1:nrow(Bmat)){
        if(length(which(Bmat[,i] == 1)) != 0){
          
          if(diet[i] == "omnivorous" & diet[j] == "prod"){
            Bmat[j,i]<-0
          }else{Bmat[j,i] <- Bmat[j,i]}

          if(diet[i] == "invP" & diet[j] != "inv"){
            Bmat[j,i] <- 0
          }else{Bmat[j,i] <- Bmat[j,i]}

           if(diet[i] == "inv" & diet[j] != "prod"){
             Bmat[j,i] <- 0
           }else{Bmat[j,i] <- Bmat[j,i]}

          if(diet[i] == "zoo" & diet[j] != "prod"){
            Bmat[j,i] <- 0
          }else{Bmat[j,i] <- Bmat[j,i]}

        }else{Bmat[,i] <- 0}
      }
    }

    return(Bmat_ref = Bmat)
  }

  if (Table == "YES") {
    for(i in 1:ncol(Bmat)){
      for (j in 1:nrow(Bmat)){
        if(length(which(Bmat[,i] == 1)) !=0 ){

          if(diet[i] == "omnivorous" & diet[j] == "prod"){
            Bmat[j,i] <- 0
          }else{Bmat[j,i] <- Bmat[j,i]}

          if(diet[i] == "invP" & diet[j] != "inv"){
            Bmat[j,i] <- 0
          }else{Bmat[j,i] <- Bmat[j,i]}

          if(diet[i] == "inv" & diet[j] != "prod"){
            Bmat[j,i] <- 0
          }else{Bmat[j,i] <- Bmat[j,i]}

          if(diet[i] == "zoo" & diet[j] != "prod"){
            Bmat[j,i] <- 0
          }else{Bmat[j,i] <- Bmat[j,i]}

        }else{Bmat[,i] <- 0}
      }
    }
    for (i in 1:ncol(Bmat)) {
      for (j in 1:nrow(Bmat)) {
        if(Bmat[j,i] == 0){
          LinksTab$Prey[LinksTab$Prey == rownames(Bmat)[j] & LinksTab$Pred == colnames(Bmat)[i]] <- NA
        }else{

        }

      }
    }
    return(list(Bmat_ref = Bmat, Table_ref=na.omit(LinksTab)))
  }

}




######################################  Weighting  ###########################################################
# Weighting procedure:
#
# Inputs:
# Niche_attributes = data frame resulting from the function "get_niche_attribute" with:
#     names = species name also used as colnames/rownames for the binary matrix
#     n = log10(species body size (?m))
#     c = optimal center of the niche (log10(?m))
#     low = lower bound of the niche range (log10(?m))
#     high = higher bound of the niche range (log10(?m))

# Bmat = initial binary interaction matrix

# Returns : the weighted interaction matrix according to the method choosen

Weighting <- function(Niche_attributes, Bmat){

    for(i in 1:ncol(Bmat)){
      j = 0
      if(length(which(Bmat[,i] == 1)) != 0){
        prob_interact <- vector()
        for(j in 1:length(which(Bmat[,i] == 1))){
        prey_range_vector <- seq(from = Niche_attributes$low[i], to = Niche_attributes$high[i], length.out=10000)
        prob <- scales::rescale(dnorm(seq(from = 5, to = 15, length.out=10000), mean = 10, sd = 2), to = c(0.1, 0.95)) # Or 0.01 and 0.99 (Pomeranz et al. 2020)
        prob <- prob[which.min(abs(Niche_attributes$n[which(Bmat[,i] == 1)[j]] - prey_range_vector))]
        prob_interact[j] <- prob
        }
        Bmat[which(Bmat[,i] == 1),i] <- prob_interact

      }else{Bmat[which(Bmat[,i] != 1),i] <- 0}
    }

  return(Bmat)
}



######################################### Rel_Ab_Mat ########################################################
# Create relative abundance matrix as the product of each species-pair density in the food web
#
# Inputs:
# Nodes_name = name of the nodes in the food web
# Nodes_Category = taxonomic category of the node (Fish, Invertebrate, Zooplankton, Primary producer)
# Nodes_Density = density of the node considered
# Bmat = binary matrix
#
# Returns :
# the quared matrix of relative abundance as SxS, with a value for each species-pair in interaction
#

# make multiple binary matrix from a weighted matrix.
Rel_Ab_Mat <- function(Nodes_name, Nodes_Category, Nodes_Density, Bmat){

  Tab <- data.frame(Node = as.character(Nodes_name),
                  Category = as.character(Nodes_Category), Density = Nodes_Density)

  # Compute each species-pair density products
  MatAb <- Bmat
  for (j in 1:ncol(MatAb)){
    for(k in 1:nrow(MatAb)){
      if(MatAb[k,j] != 0){
        MatAb[k,j] <- Tab$Density[which(Tab$Node == colnames(MatAb)[j])] * Tab$Density[which(Tab$Node == rownames(MatAb)[k])]
        }
      }
  }

  # Rescaled it from 0.5 to 1 based on local relative abundances (i.e. min and max found considering all lakes)
  # Rescaling is different if the interaction is between one fish and the rest or one inv/zoo and phyto as units
  # of abondance are not the same

  vecFish <- as.character(Tab$Node[Tab$Category == "Fish"])
  vecInv <- as.character(Tab$Node[Tab$Category != "Fish"])

  RangeFish <- c(MatAb[rownames(MatAb) %in% vecInv, colnames(MatAb) %in% vecFish])[c(MatAb[rownames(MatAb) %in% vecInv, colnames(MatAb) %in% vecFish]) > 0]
  RangeInv <- c(MatAb[rownames(MatAb) %in% vecInv, !colnames(MatAb) %in% vecFish])[c(MatAb[rownames(MatAb) %in% vecInv, !colnames(MatAb) %in% vecFish]) > 0]

  MatAb2 <- MatAb
  for (j in 1:ncol(MatAb)){
    for(k in 1:nrow(MatAb)){

      if(MatAb2[k,j] != 0){

        if(colnames(MatAb2)[j] %in% vecFish & rownames(MatAb2)[k] %in% vecInv){
          MatAb2[k,j] <- rescale(MatAb[k,j], to = c(0.5,1),
                               from = range(RangeFish, na.rm = TRUE))

        }else{
          MatAb2[k,j] <- rescale(MatAb[k,j], to = c(0.5,1),
                               from = range(RangeInv, na.rm = TRUE))
          }
      }
    }
  }

  return(as.matrix(MatAb2))
  }



########################################## Make_Bern #######################################################
# Create multiple binary matrices from on Bernouilli trials based on a matrix of interaction probability
# here computed from the product of the weighted matrices and the relative abundance matrices
# /!\ This function necessitates to load the package Rlab
# This function was already established by Vagnon et al. 2022
#
# Inputs :
# n = number of matrices to generate
# Wmat = squared matrix of interaction probabilities
#
# Returns :
# a large array of dimension n, nrow(Wmat), ncol(Wmat) where each n = a simulated interaction matrix
# filled with 0 or 1 generated from the probabilities Wmat implemented in the Bernouilli trials

library(Rlab)

# make multiple binary matrix from a weighted matrix.
Make_Bern <- function(n, Wmat){
  mat_ber <- array(NA, c(n, ncol(Wmat), nrow(Wmat)))
  for(k in 1:n){
    mat_inter <- Wmat

    for (j in 1:ncol(Wmat)){
      for (i in 1:nrow(Wmat)){

        mat_inter[i,j] <- rbern(1, Wmat[i,j])
      }
    }
    mat_ber[k,,] <- mat_inter
    dimnames(mat_ber) <- list(seq(1, n, 1), colnames(Wmat), rownames(Wmat)) 
  }

  return(mat_ber)
}


###########################################   Jacobian_Mat  ######################################################
# Create a jacobian binary matrix using the function jacobian_binary of Sauve et al (2016) modified by Pomeranz et al. (2020)
# Allows to infer interaction strengths between species-pairs based on a half normal distribution.
# Can be totally randomly distributed or can be symetric between J_ij (pred on prey and < 0) and J_ji (prey on pred and >0)
#
# Inputs :
# m = adjacency matrix from the make_bern function
# Nodes = node names
# Method = method to use from Pomeranz et al. (2020). Can be Random for the random interaction strengths or other
#
# Returns :
# If Method = "Random":
#    returns a jacobian binary matrix with interaction strengths drawn randomly following a half normal distribution with
#    values different on the upper (prey on pred and >0) and the lower triangle (pred on prey and < 0) with J_ij different J_ji
# If Method != "Random":
#    a jacobian binary matrix with interaction strength drawn randomly following a half normal distribution, antisymetrically
#    distributed on the upper and lower triangle, with J_ij = J_ji*-1

Jacobian_Mat <- function(m, Nodes, Method){

if (dim(m)[1] == dim(m)[2]){ # Is m a square matrix?
  if(Method == "Random"){

      J <- t(m)
      # transpose matrix so negative interactions are effects
      # of predators in colums on prey in rows and
      # positive interactions are effects of prey in rows
      # on predator in column
      J[which(J < t(J))] <- -t(J)[which(J < t(J))]

      L <- sum(abs(J))/2
      strengthPreyOnPred <- rnorm(L, sd = 0.1) + 1 # L values for interaction strength magnitude drawn from a normal distribution
      upper_tri_str_index <- which((upper.tri(J) == TRUE) & (J != 0), arr.ind = TRUE) # non-zero elements in the upper triangle matrix
      J[upper_tri_str_index] <- J[upper_tri_str_index]*strengthPreyOnPred

      lower_tri_str_index <- cbind(upper_tri_str_index[, 2], upper_tri_str_index[, 1]) # so indices match between upper triangle and lower triangle elements
      strengthPredOnPrey <- rnorm(L, sd = 0.1) + 1
      J[lower_tri_str_index] <- J[lower_tri_str_index]*strengthPredOnPrey

      dimnames(J)<-list(Nodes, Nodes)
      return(J)

    } else if(Method != "Random"){
      J <- t(m)
      # Modified from original function
      # transpose matrix so negative interactions are effects
      # of predators in colums on prey in rows and
      # positive interactions are effects of prey in rows
      # on predator in column
      J[which(J < t(J))] <- -t(J)[which(J < t(J))]

      L <- sum(abs(J))/2
      strength <- rnorm(L, sd = 0.1) + 1 # L values for interaction strength magnitude drawn from a normal distribution
      #strength <- rhalfnorm(L, theta=sqrt(pi/2))
      upper_tri_str_index <- which((upper.tri(J) == TRUE) & (J != 0), arr.ind = TRUE) # non-zero elements in the upper triangle matrix
      lower_tri_str_index <- cbind(upper_tri_str_index[, 2], upper_tri_str_index[, 1]) # so indices match between upper triangle and lower triangle elements

      J[upper_tri_str_index] <- J[upper_tri_str_index]*strength
      J[lower_tri_str_index] <- J[lower_tri_str_index]*strength

      dimnames(J)<-list(Nodes, Nodes)
      return(J)
    }

  }else {
    stop("m is not a square matrix.")
  }
}




#######################################    Scale_JmatBS  ######################################################
# Use the jacobian matrix from the function "Jacobian_Mat" and scale it by predator/prey body size ratios
# e.g., smallest positive and greatest negative effects between large predators and small prey
# Values on the lower triangle are the interaction strength J_ij multiplicated by the pred/prey body size ratio
# Values on the upper triangle are the interaction strength J_ji multiplicated by the prey/pred body size ratio
#
# Inputs :
# - JmatPom = jacobian matrix from the function "Jacobian_Mat" which considered random interaction strength
# - Nodes = prey body size corresponding to links present in the binary matrix
# - Node_Category = Category of each node in the food web
# - Node_size = Body size (?m) of of each node in the food web
# - RangeVert = range of body size vertebrate predator/prey
# - RangeInvert = range of body size predatory invertebrate /prey
#
#Returns :
# a jacobian binary matrix with interaction strength scaled to predator/prey body size ratios accounting a different scaling
# regarding the category of the predator (Fish or not fish)

Scale_JmatBS <- function(Jmat, Nodes, Node_Category, Node_Size, RangeVert, RangeInvert){

  library(scales)

  options(warn = -1)

  JJ <- Jmat
  RangeVertRev = range(1/RangeVert)
  RangeInvertRev = range(1/RangeInvert)
  lower_tri_str_index <- which((lower.tri(JJ) == TRUE) & (JJ != 0), arr.ind = TRUE)
  upper_tri_str_index <- which((upper.tri(JJ) == TRUE) & (JJ != 0), arr.ind = TRUE)

  # Is there Fish in Node_Category & Is there Fish in the non-null Jacobian columns (preds)
  # if(is.element(Node_Category, "Fish") & length(lower_tri_str_index[Node_Category[lower_tri_str_index[, 2]] == "Fish", 2]) != 0){
  if(any(is.element(Node_Category, "Fish")) & length(lower_tri_str_index[Node_Category[lower_tri_str_index[, 2]] == "Fish", 2]) != 0){ 
    # Jacobian matrix filled with of body size ratio
    # Size of preds which are Fish / Size of their prey
    Pred_PreyFish <- Node_Size[lower_tri_str_index[Node_Category[lower_tri_str_index[, 2]] == "Fish", 2]] /
      Node_Size[lower_tri_str_index[Node_Category[lower_tri_str_index[, 2]] == "Fish", 1]]
    # Size of preds which are not Fish / Size of their prey
    Pred_PreyOther <- Node_Size[lower_tri_str_index[Node_Category[lower_tri_str_index[, 2]] != "Fish", 2]] /
      Node_Size[lower_tri_str_index[Node_Category[lower_tri_str_index[, 2]] != "Fish", 1]]

    JJ[lower_tri_str_index] <- c(Pred_PreyFish, Pred_PreyOther) # Pred size/prey size
    JJ[upper_tri_str_index] <- c(1 / Pred_PreyFish, 1 / Pred_PreyOther) # reverse of Pred size/prey size

    # Rescale the body size ratio between 0.25 and 1.25 and differently weither consumer is a vertebrate or not
    ScalePred_PreyFish <- rescale(Pred_PreyFish, to = c(0.25, 1.25),   #range of interaction strength from rnorm = range (JmatPom[JmatPom>0])
                                from = RangeVert)
    ScalePred_PreyOther <- rescale(Pred_PreyFish, to = c(0.25, 1.25),   #range of interaction strength from rnorm = range (JmatPom[JmatPom>0])
                                from = RangeInvert)

    JJ[lower_tri_str_index] <- c(ScalePred_PreyFish, ScalePred_PreyOther)

    # Do the same for the reverse
    RevScalePred_PreyFish <- rescale(1 / Pred_PreyFish, to = c(0.25,1.25),   #range of interaction strength from rnorm = range (JmatPom[JmatPom>0])
                                   from = RangeVertRev)
    RevScalePred_PreyOther<-rescale(1 / Pred_PreyOther, to = c(0.25,1.25),   #range of interaction strength from rnorm = range (JmatPom[JmatPom>0])
                                    from = range(RangeInvertRev))

    JJ[upper_tri_str_index] <- c(RevScalePred_PreyFish, RevScalePred_PreyOther)

  } else{
    # Jacobian matrix filled with body size ratio
    Pred_PreyOther <- Node_Size[lower_tri_str_index[Node_Category[lower_tri_str_index[, 2]] != "Fish", 2]] /
      Node_Size[lower_tri_str_index[Node_Category[lower_tri_str_index[,2]] != "Fish", 1]]

    JJ[lower_tri_str_index] <- Pred_PreyOther # Pred size/prey size
    JJ[upper_tri_str_index] <- 1 / Pred_PreyOther # reverse of Pred size/prey size

    # Rescale the body size ratio between 0.25 and 1.25 and differently weither consumer is a vertebrate or not
    ScalePred_PreyOther <- rescale(Pred_PreyOther, to = c(0.25, 1.25),   #range of interaction strength from rnorm = range (JmatPom[JmatPom>0])
                                 from = RangeInvert)
    JJ[lower_tri_str_index] <- ScalePred_PreyOther

    # Do the same for the reverse
    RevScalePred_PreyOther <- rescale(1 / Pred_PreyOther, to = c(0.25, 1.25),   #range of interaction strength from rnorm = range (JmatPom[JmatPom>0])
                                    from = range(RangeInvertRev))

    JJ[upper_tri_str_index] <- RevScalePred_PreyOther
  }

  JJ2 <- JJ * Jmat
  dimnames(JJ2) <- list(Nodes, Nodes)
  return(JJ2)
}


#########################################   Stability  #########################################################
# Compute stability metric (function from Sauve et al. 2016) based on the qunatity of self-regulation to add on the diagonal
# of the Jacobian matrix to obtain a  system mathematically stable (i.e., real part of the greatest eigenvalue of J-s2*I is negative)
#
# Inputs :
# J = jacobian matrix from the function "Scale_JmatBS"
# s2 = arbitrary parameter to add on the diagonal
#
#
#Returns :
# the self-regulation value (the lowest it is the more stable is the system).

Stability <- function(J, s2){

  test1 <- (dim(J)[1] == dim(J)[2]) # Is J a square matrix?
  test2 <- FALSE
  if (test1 == TRUE){
    S <- dim(J)[1]
    test2 <- which(diag(J) != vector("numeric", S)) # Does J have a null diagonal?
  }

  if ((test1 == TRUE) & (length(test2) == 0)){ # if J is a square matrix with a null diagonal
    S <- dim(J)[1] # S is the number of species in the network
    s1 <- 0
    I <- diag(S)
    E1 <- max(Re(eigen(J-s1*I, only.values = T)$values))
    E2 <- max(Re(eigen(J-s2*I, only.values = T)$values))

    if ((E1 >= 0) & (E2 < 0)){ # if s2 is well chosen and the system is not already stable
      while ((s2-s1) >= 10^-4){
        stab <- (s1+s2)/2
        E1 <- max(Re(eigen(J - stab*I, only.values = T)$values))
        if (E1 >= 0){
          s1 <- stab
        }
        else {
          s2 <- stab
        }
      }
      return(stab)
    }

    if (E1 < 0){
      # stop("J corresponds to a stable system.")
      stab <- 0
      return(stab)
    }
    if (E2 >= 0){
      # stop("s2 is not high enough.")
      stab <- NA
      return(stab)
    }
  }
  else { # if J is not a square matrix with a null diagonal
    if (test1 == FALSE){
      #stop("J is not a square matrix.")
      stab <- "J is not square"
      return(stab)
    }
    if (length(test2) > 0){
      #stop("J does not have a null diagonal.")
      stab <- "diag not null"
      return(stab)
    }
  }
}



########## REFERENCES cited in the script
# - Gravel, D., Poisot, T., Albouy, C., Velez, L., Mouillot, D., and Freckleton, R. (2013). Inferring food web structure
#         from predator-prey body size relationships. Methods Ecol. Evol. 4(11), 1083-1090. doi: 10.1111/2041-210x.12103.

# - Pomeranz, J.P.F., Wesner, J.S., and Harding, J.S. (2020). Changes in stream food-web structure across a gradient of
#         acid mine drainage increase local community stability Ecology 101(9), e03102.

# - Sauve, A.M.C., Th?bault, E., Pocock, M.J.O., and Fontaine, C. (2016). How plants connect pollination and herbivory
#          networks and their contribution to community stability. Ecology 97(4), 908-917. doi: 10.5061/dryad.3s36r118.

# - Vagnon, C., Cattan?o, F., Guillard, J., and Frossard, V. (2022). Inferring the trophic attributes and consequences of
#          co-occurring lake invaders using an allometric niche model. Biological Invasions. doi: 10.1007/s10530-022-02745-2.

# - Vagnon, C., Cattan?o, F., Goulon, C., Grimardias, D., Guillard, J., and Frossard, V. (2021). An allometric niche model
#          for species interactions in temperate freshwater ecosystems. Ecosphere 12(3), e03420.


