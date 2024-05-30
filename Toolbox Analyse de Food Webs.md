# Toolbox Analyse de Food Webs



## Inférence de food webs (Chloé)



## Stabilisation du réseau (?)

- Appliquer une dynamique sur le réseau, si au moins une espèce ne persiste pas ajouter une self-regulation et re-run ; si encore extinction augmenter la self-regulation ; jusqu'à atteindre un food web stable. Permet de définir une self-regulation



## Identification et visualisation des chaînes trophiques

- Calcul niveaux trophiques et omnivorie
- Caractérisation du réseau : Collectivité, Connectance, Omnivorie totale, Niveau trophique maximal
- Identifier et stocker dans un data.frame chaque chaîne trophique : top (unique) / n-1 / n-2 (unique)
- Noter dans la colonne suivante si il y a respect d'une cascade n-step (en gros s'il y a de l'omnivorie inhibant au sein de chaque chaîne même la n-step cascade)
- Graph du food web avec possibilité de choisir une espèce et mise en lumière de ses liens dans le réseau par rapport aux autres liens (rouges contre noirs pour le reste par exemple) / Et autre possibilité, si on choisit pas de mettre une espèce en lumière alors les liens sont différenciés entre ceux où il y a cascade n-step et ceux où il n'y a pas (à cause d'omnivorie) / Mise en lumière d'une chaîne trophique avec-sans son intégration dans le réseau 
- **Graph** --> choisir une espèce, partout où cette espèce est présente dans EdgesAll changer la couleur. Puis regarder à qui elle est relié en ordre deux, et pareil dans EdgesAll, là où il y a ces espèces d'ordre deux, avec une espèce d'ordre 1 changer la couleur. Comme ça l'utilisateur choisi une espèce, et sur le graph ses liens directs, et ses liens d'ordre deux apparaissent avec deux couleurs différents du reste du graph.



## Analyse des cascades de chaque chaîne

- Ajout au data.frame le ratio Net-Cascade/Direct-Cascade et Cascade/Attenuation/Amplification/Inversion
- Ajout de la mesure de collectivité expériencée par chaque chaîne
- Upgrade visualisation du réseau : 4 réseaux côte à côte avec pour chacun les liens mis en lumière pour chaînes avec Cascade/Att/Amp/Inv --> graph du food web, identifier les prédateurs, calculer ratio de cascade pour chaque chaîne et classer les chaînes selon le type, dans EdgeAll attribué une couleur différente selon le type de la chaîne (faudrai des vecteurs genre PredCascade, PredAtt ...., ConsoCascade, ConsoAtt .... ; et ensuite pouvoir faire des ifelse sur EdgesAll, genre ifelse un pred et une prey sont chacun dans un vecteur PredAtt et PreyAtt alors mettre une couleur rouge pour le lien).
- Calcul community-cascade inversion : Oui/Non



## Dynamique de perturbation de cascade

- Possibilité de choisir un noeud (ou plusieurs), appliquer une perturbation positive ou négative (scalée pour assurer non-extinction) et voir réponse du système
- Plot de la réponse dynamique du système





# To-do

- igraph ajout paramètre pour layout : par niveau troph ou bien cercle ou ...
- dynamique de perturb ajouter un print qui donne les espèces qui s'éteignent s'il y en a