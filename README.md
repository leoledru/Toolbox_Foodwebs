# Toolbox Food Webs Analysis

| Function | Description |
| --- | --- |
| `Toolbox_output.Rmd` | The main RMarkdown file, combining detailed writing and R scripts. All the analysis are carried out within this file, using the various R functions available as .R files in the github folder.  |
| `CheckInit.R` | Check if the interaction matrix representing the food web is correct and ready for further analysis. |
| `Heatmaps.R` | Heatmaps of direct and net effects, as a reminder an interaction is read from column j to row i. |
| `ChainCount.R` | Identifies all the trophic chains, i.e. number of paths from apex predator(s) (predator species without predator) to basal species, stores them and measures their length. Gives also the mean length and the standard-deviation. |
| `FoodWebAnalysis.R` | Identification of trophic chains and computation of several proxies: 1) Collectivity, 2) Connectance, 3) Average omnivory, 4) a sub-list with measures characterizing the trophic cascade process expressed by each trophic chain in the food web. |
| `ComputeLinks.R` | For each species find all non-null interactions for each order until a maximum order set by the user. Also find the cumulative *orderLim* from which each species has interacted with all others in the food web. |
| `FoodWebGraph.R` | Vizualisation of the food web as a graphic where species are nodes and interactions are links |
| `InvertedChainGraph.R` | Identification and graphical vizualisation of trophic chains with a long-term inversion of trophic cascade |
| `StabilityAlgorithm.R` | Compute stability metric (function from Sauve et al. 2016) based on the quantity of self-regulation to add on the diagonal of the "Jacobian" matrix to obtain a system mathematically "stable" (i.e., real part of the greatest eigenvalue of J-s2*I is negative) |
| `FoodWebEquilibrium.R` | Simulates food web dynamics with a Generalized-Lokta-Volterra model to find its stable configuration. Remove extinct species until find an equilibrium with only persistant species. |
| `FoodWebPerturb.R` | Simulates a stable system until equilibrium, then perturb the chosen species (one or several) and simulates the system dynamics response |
