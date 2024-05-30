###############################################################################################################
#               Infer trophic links in the metaweb, reconstruct food webs  and compute stability              #
#                                                                                                             #
#VAGNON1,2, POMERANZ3, LOHEAC4, VALLAT4, GUILLARD1,2, RAYMOND5, SENTIS2,6, FROSSARD1,2                        #
#                                                                                                             #
#1Univ. Savoie Mont Blanc, INRAE, UMR CARRTEL, 74200 Thonon-les-Bains, France.                                #
#2Pôle R&D Ecosystèmes Lacustres (ECLA), OFB-INRAE-USMB, 13182 Aix-en-Provence, France.                       #
#3Colorado Mesa Univ., CO 81501, Colorado, United states.                                                     #
#4Fédération de Savoie pour la Pêche et la Protection du Milieu Aquatique (FDPPMA 73),                        #
#  73230 Saint-Alban-Leysse, France                                                                           #
#5Office Français pour la Biodiversité, Unité Spécialisée Milieux Lacustres, 74200 Thonon-les-Bains, France.  #
#6INRAE, Aix Marseille Univ., UMR RECOVER, 13182 Aix-en-Provence, France.                                     #
#2022                                                                                                         #
#                                                                                                             #
# This script was created for guiding the futur users of the framework developped in the study :             #
#   ""               #
# The user needs to download the test data set "DATA.Rdata", the functions in "cv_Functions_AltitudeLakes.R"  #
# and the quantile regression parameters vor vertebrates () and invertebrates ().                             #
#
# This script is devided in 3 main parts and different sub-parts:
#      I) Trophic link inferences
#        It mixes the aNM (Vagnon et al. 2021) for inferring all possible trophic links and weighting trophic #
#        links (i.e. Niche probability) and the method of Pomeranz et al. (2020) for refining trophic links   #
#        according to node densities depending on lakes richness (i.e. Neutral process) and for generating    #
#        100 binary matrices based on Bernouilli trials for considering variability in consumers diet.        #
#                                                                                                             #
#         4 sub-parts:                                                                                        #
#          a) Infer trophic niche for consumers and the corresponding squared binary matrix                   #
#          b) Refine trophic links according to the trophic characteristics of consumers                      #
#          c) Obtain interaction probabilities based on a weighting procedure and species abundances          #
#          d) Generate adjacency matrices based on the interaction probabilities                              #
#                                                                                                             #
#                                                                                                             #
#      II) Jacobian matrix inferences                                                                         #
#          The transformation of binary adjacency matrices into Jacobian matrices is based on a method from   #
#          Pomeranz et al. 2020 and modified for conducting the present study. Interaction strengths are      #
#          retrieved from a half-normal distribution and attributed randomly on the upper triangle in the     #
#          Jacobian matrix. The values on the lower triangle are symetric                                     #
#                                                                                                             #
#      III) Stability assessment                                                                              #
#           The stability measure is based on a modified function from Sauve et al. (2016). Stability is      #
#           computed from an iterative process in which the values on the diagonal                            #
#                                                                                                             #
#                                                                                                             #
#  References:                                                                                                #
# Pomeranz JPF, Wesner JS, Harding JS (2020) Changes in stream food-web structure across a gradient of acid   #
#       mine drainage increase local community stability.Ecology 101:e03102                                   #
#                                                                                                             #
# Sauve, A. M. C., E. Thébault, M. J. O. Pocock, and C. Fontaine. 2016. How plants connect pollination and    #
#       herbivory networks and their contribution to community stability. Ecology 97:908-917.                 #
#                                                                                                             #
# Vagnon C, Cattaneo F, Goulon C, Grimardias D, Guillard J, Frossard V (2021) An allometric niche model for   #
#       species interactions in temperate freshwater ecosystems Ecosphere 12:e03420                           #
#                                                                                                             #
###############################################################################################################

# Set the working directory in which you want to work
setwd("~/Documents/Post_Doc_CARRTEL/R")

# Packages to load for applying functions
library(stringr)
library(Rlab)
library(scales)

# Load functions and quantile regressions for defining vertebrate and invertebrate prey body size range
source("cv_Functions_AltitudeLakes.R") # aNM functions
load("Param_reginvert.Rdata") # aNM parameters for invertebrate consumers
load("Param_regvert.Rdata") # aNM parameters for vertebrate consumers

# Load data
load("DATA.Rdata") # dataset of nodes to use in the metaweb

#### 1) Infer trophic niche for consumers and the corresponding squared binary matrix ####
# /!\ Don't forget to make a log10 on the body size if it is not already the case

# Obtain trophic niche parameters from the function Niche_Attributes
Niche <- Niche_Attributes(species_name = DATA$Node, body_size = log10(DATA$BodySize), DATA$Phylum)

# Obtain the Binary matrix and the trophic link table from the function Binary_Mat
All_Trophic_Links <- Binary_Mat(Niche$name, Niche$n, Niche$c, Niche$low, Niche$high, "YES") # Fill with yes if the table of links is needed
Bmat <- All_Trophic_Links[[1]] # Get the binary matrix
Links <- All_Trophic_Links[[2]] # Get the table of links

#### 2) Refine trophic links according to the trophic characteristics of consumers using the function Ref_L_Diet ####
Ref_Links <- Ref_L_Diet(Bmat = Bmat, diet = DATA$Trophy, Table = "YES", LinksTab = Links)

#### 3) Interaction probability matrix
# A) Weight trophic links according to the position of prey body size in the range of prey body size that a consumer can consume
# e.g. prey body size falling at the center of the range has highest probabilities than prey body sizes at the range bounds
Mat_W <- Weighting(Niche_attributes = Niche, Bmat = Ref_Links[[1]]) # similar to the Niche probability matrix in Pomeranz et al.2020

#### 4)  Get the Binary matrices and the weigthed matrices for each of the 18 lakes ####
BinMat <- Ref_Links[["Bmat_ref"]]

# B) Use species density for refining trophic links based on neutral processes with the function Rel_Ab_Mat
# i.e. 2 rare species with low densities have lowest probability to interact
AbMat <- Rel_Ab_Mat(DATA$Node, DATA$Category, DATA$Density, BinMat)


# C) Calculate the final interaction probability matrix as the product of the weighted links and the densities or abundances
# i.e. Lakes_Res[["Wmat"]] * Lakes_Res[["AbMat"]]
IntMat <- AbMat* Mat_W
IntMatProb <- IntMat
  for (j in 1:ncol(IntMat)){
    for(k in 1:nrow(IntMat)){
      if(IntMatProb[k,j] != 0){

        # Rescale it between 0.01 and 0.99
        IntMatProb[k,j] <- rescale(IntMat[k,j], to = c(0.01,0.99), from = range(IntMat[which(IntMat != 0)], na.rm = TRUE))
      }
    }
  }


#### 6)  Generate adjacency matrices based on Bernouilli trials with the function Make_Bern ####
BernMat <- Make_Bern(n=150, Wmat=IntMatProb)




# Load packages
library(dplyr)
library(stringr)

#Load functions
source("cv_Functions_AltitudeLakes.R")

########################################   1) INTERACTION STRENGTHS  #########################################
# 1) Transform adjacancy binary matrices into jacobian matrix including interaction strengths
Fish <- droplevels(DATA$Node[DATA$Category == "Fish"])
NoFish <- droplevels(DATA$Node[!DATA$Node %in% Fish])

RangeVert <- range((10^Links$Log_Size_Pred[Links$Pred %in% Fish]) / (10^Links$Log_Size_Prey[Links$Pred %in% Fish]))
RangeInvert <- range((10^Links$Log_Size_Pred[!Links$Pred %in% Fish]) / (10^Links$Log_Size_Prey[!Links$Pred %in% Fish]))

# Compute jacobian matrix for each adjacency matrix
Lakes_Stability <- vector("list")
n = 0
set.seed(n + 1)
for (j in 1:dim(BernMat)[1]){
  AdjMat <- BernMat[j,,]

  # Method 2 in Pomeranz et al 2019 (scaled to  pred:pred bodysize ratio)
  # sample interaction strengths from a half normal distribution
  JmatPom <- Jacobian_Mat(m = AdjMat, Nodes = DATA$Node, Method = "Method_2") # Using the same values for + and - interactions
  #JmatPom <- JacobMat # If you want the same randomly distributed interaction strength as in the Method 1

  #SCale interaction strengths by predator:prey body size ratios
  JacobMatBS <- Scale_JmatBS(Jmat = JmatPom,
                           Nodes = DATA$Node,
                           Node_Category = DATA$Category,
                           Node_Size = DATA$BodySize,
                           RangeVert = RangeVert,
                           RangeInvert = RangeInvert)

  JacobMatBSNoSalmo <- JacobMatBS[NoFish, NoFish] # Consider matrices without salmonids and forage fish

  # Store results in the list
  Lakes_Stability[["JacobMatBS"]][[j]] <- JacobMatBS
}


########################################   2) FOOD WEB STABILITY  ############################################
# s = stability metrics = was defined as the minimum amount of intraspecific competition (e.g., the diagonal
# of the Jacobian matrix, Jii) necessary for a food-web iteration fo be stable
# --> stability() function from the supplementary method of Sauve et al. 2016
# lower value of s = more stable than high value of s (only comparison between systems)
# Values on the diagonal of the Jacobian matrices (i.e., intraspecific competition), were varied until the
# individual matrix was stable (e.g., all of the real parts of the eigenvalues were negative) and s is still
# the same for all values in the diagonal of all matrices.

Stab_Metric <- data.frame(Lake = NA, Simu = NA, S = NA)
for (j in 1:dim(BernMat)[1]){
  
  Stab_JacobMatBS <- Stability(Lakes_Stability[["JacobMatBS"]][[j]], s2 = 2)

  Stab <- data.frame(Lake = unique(DATA$Lake_name), Simu = j, S = Stab_JacobMatBS)

  Stab_Metric <- rbind(Stab_Metric, Stab)
}
Lakes_Stability[["Stability"]] <- na.omit(Stab_Metric)




################################################################################
library(igraph) # pour graph de réseau
library(NetIndices) # pour calcul niveau troph
library(cheddar) # pour graph de réseau (pas testé)
library(matlib) # pour calcul inverse de matrice
library(latex2exp)
library(plotrix) # pour plot cercle
library(expm)
library(sparsevar)
source("jacobian_analysis.R")
source("jacob_plot.R")

# Analyse des Jacobienne
Jacobs <- Lakes_Stability[[1]]
D = rep.int(-1, nrow(JacobMatBS))
JacobianAnalysis <- list()
for (i in 1:length(Jacobs)){
  OneJacob <- matrix(unlist(Jacobs[[i]]), nrow = nrow(JacobMatBS), ncol = ncol(JacobMatBS))
  diag(OneJacob) = D 
  JacobianAnalysis[[i]] <- jacobian_analysis(OneJacob, D)
}

# Visualisation résultats
# setEPS()
# postscript("eigen_spectral_radius.eps")

collectivity = c()
for (i in 1:length(Jacobs)){
# for (i in 1:2){
  results <- JacobianAnalysis[[i]]
  ev <- results[[1]]
  A_ij <- results[[2]]
  val <- ev$values
  collectivity[i] <- jacob_plot(A_ij, val, i)
}
title(main = paste("Mean collectivity = ", round(mean(collectivity),2)))
dev.off()




