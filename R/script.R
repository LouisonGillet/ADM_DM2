# En tête ----
# Desc : DM2 ADM
# Date: 17/11/2024
# Auteur : EL MAZZOUJI Wahel & GILLET Louison 

# Chargement des données ----
Datagenus <- read.csv("data/Datagenus.csv", sep=";")

# Préparation des données ----
data <- Datagenus[1:1000, ] # On ne prend pas en compte la ligne 1001
especes <- paste0("gen", 1:27) # Noms des colonnes des espèces
colonnes_selectionnees <- c(especes, "surface", "forest", "geology") # Colonnes à conserver
data <- data[, colonnes_selectionnees]

# Calcul des densités de peuplement par unité de surface ----
densite_peuplement <- as.matrix(data[especes] / data$surface) # Conversion en matrice 

# Standardisation si nécessaire ----
n <- nrow(densite_peuplement) # Nombre de parcelles
p <- ncol(densite_peuplement) # Nombre d'espèces

# Calcul direct de la matrice des distances
# Initialiser une matrice pour stocker les contributions par variable
q = 27   # Nombre de caractéristiques
contributions <- array(0, dim = c(n, n, q)) 
# Calcul des distances et contributions
for (i in 1:n) {
  for (j in i:n) {
    # Différence au carré pour chaque variable
    diff <- (densite_peuplement[i, ] - densite_peuplement[j, ])^2
    # Distance totale (somme des différences au carré)
    dist <- sum(diff)
    # Contribution de chaque variable
    if (dist != 0) {
      contributions[i, j, ] <- round((diff / dist)*100,2)
      
    }
  }
}
print(contributions[1, 800,])

# Centrage et réduction des densités
# Calcul des moyennes et écarts-types par espèce
moyennes_especes <- colMeans(densite_peuplement)
sd_especes <- sqrt(colSums((densite_peuplement - matrix(moyennes_especes, n, p, byrow = TRUE))^2) / (n))
densite_centree_reduite <- (densite_peuplement - matrix(moyennes_especes, n, p, byrow = TRUE)) / 
  matrix(sd_especes, n, p, byrow = TRUE)

datapeuple <- as.data.frame(densite_centree_reduite)

# Calcul des distances euclidiennes ----
dp <- dist(datapeuple, method = "euclidean")  # Matrice des distances

# CAH avec Ward ----
CAHWard <- hclust(d = dp, method = "ward.D")  # CAH avec la méthode de Ward
CAHMax <- hclust(d = dp, method = "complete")  # CAH avec la méthode de Max
CAHMoy <- hclust(d = dp, method = "mcquitty")  # CAH avec la méthode de Moy
CAHMin <- hclust(d = dp, method = "single")  # CAH avec la méthode de Min 

plot(CAHWard, labels = FALSE, main = "Dendrogramme - Méthode de Ward")
plot(CAHMax, labels = FALSE, main = "Dendrogramme - Méthode Complete (Max)")
plot(CAHMoy, labels = FALSE, main = "Dendrogramme - Méthode McQuitty (Moy)")
plot(CAHMin, labels = FALSE, main = "Dendrogramme - Méthode Single (Min)")

#### Nombre de classes #### 
k <- 8

# Partition en k=4 classes ----
PWard = cutree(tree = CAHWard, k) #Ward
PMax = cutree(tree = CAHMax, k) #Sautmax
PMoy = cutree(tree = CAHMoy, k) #Sautmoy
PMin = cutree(tree = CAHMin, k) #Sautmin

# Calcul des R2 pour chaque variable ----
R2_PWard = cbind(rep(0 , ncol(datapeuple)))
R2_PMax = cbind(rep(0 , ncol(datapeuple)))
R2_PMoy = cbind(rep(0 , ncol(datapeuple)))
R2_PMin = cbind(rep(0 , ncol(datapeuple)))

for (i in 1:ncol(datapeuple)) {
  R2_PWard[i] = summary(lm(datapeuple[,i] ~ as.factor(PWard)))$r.squared
  R2_PMax[i] = summary(lm(datapeuple[,i] ~ as.factor(PMax)))$r.squared
  R2_PMoy[i] = summary(lm(datapeuple[,i] ~ as.factor(PMoy)))$r.squared
  R2_PMin[i] = summary(lm(datapeuple[,i] ~ as.factor(PMin)))$r.squared
}

row.names(R2_PWard) = colnames(datapeuple)
row.names(R2_PMax) = colnames(datapeuple)
row.names(R2_PMoy) = colnames(datapeuple)
row.names(R2_PMin) = colnames(datapeuple)

# Calcul du R2 global ----
R2G_PWard = mean(R2_PWard)
R2G_PMax = mean(R2_PMax)
R2G_PMoy = mean(R2_PMoy)
R2G_PMin = mean(R2_PMin)

# Fonction pour calculer l'indice de Rand
rand_index <- function(partition1, partition2) {
  contingency_table <- table(partition1, partition2)
  n <- sum(contingency_table)
  total_pairs <- choose(n, 2)

  same_cluster_pairs <- sum(choose(contingency_table, 2))
  
  different_cluster_pairs <- total_pairs - sum(apply(contingency_table, 1, sum)^2) / 2 -
    sum(apply(contingency_table, 2, sum)^2) / 2 +
    same_cluster_pairs
  
  # Calcul de l'indice de Rand
  rand_index <- (same_cluster_pairs + different_cluster_pairs) / total_pairs
  
  return(rand_index)
}

# Indices de rand 
rand_Ward_Max <- rand_index(PWard, PMax)
rand_Ward_Moy <- rand_index(PWard, PMoy)
rand_Ward_Min <- rand_index(PWard,PMin)
rand_Max_Moy <- rand_index(PMax, PMoy)
rand_Max_Min <- rand_index(PMax,PMin)
rand_Moy_Min <- rand_index(PMoy,PMin)

# Création des matrices d'indicatrices pour chaque partition ----
ICWard = data.frame(model.matrix(~as.factor(PWard) - 1))
ICMax = data.frame(model.matrix(~as.factor(PMax) - 1))
ICMoy = data.frame(model.matrix(~as.factor(PMoy) - 1))
ICMin = data.frame(model.matrix(~as.factor(PMin) - 1))

mICWard = as.matrix(ICWard)
mDP = as.matrix(datapeuple)
CentresC2 = solve(t(mICWard) %*% mICWard) %*% t(mICWard)%*% mDP
KMDP2 = kmeans(datapeuple, CentresC2)
KMDP2$cluster
boxplot(datapeuple[,27]~as.factor(PWard))
boxplot(datapeuple[,27]~as.factor(KMDP2$cluster))

# Sauvegarde pour Rmd ---- 
save.image(file = "ressources/prepa.RData")

