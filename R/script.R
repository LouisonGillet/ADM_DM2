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
CAHDP <- hclust(d = dp, method = "ward.D")  # CAH avec la méthode de Ward
plot(CAHDP,labels=FALSE)  # Dendrogramme de la hiérarchie de la CAH

# Calcul des R2 pour chaque variable ----
R2_PDP4 = cbind(rep(0 , ncol(datapeuple)))
for (i in cbind(1:ncol(datapeuple))) {
  R2_PDP4[i] = summary(lm(datapeuple[,i]~as.factor(PDP4)))$r.squared
}
row.names(R2_PDP4) = colnames(datapeuple)

# Calcul du R2 global ----
R2G_PV4 = mean(R2_PDP4)

# Résultats finaux ----
print("R2 global pour la partition en 4 classes :")
print(R2G_PV4)

# CAH avec Max ----
CAHDP2 <- hclust(d = dp, method = "complete")  # CAH avec la méthode de Max
plot(CAHDP2,labels=FALSE)  # Dendrogramme de la hiérarchie de la CAH

# CAH avec moy ----
CAHDP3 <- hclust(d = dp, method = "mcquitty")  # CAH avec la méthode de Max

plot(CAHDP3,labels=FALSE) # Dendrogramme de la hiérarchie de la CAH

# Partition en k=4 classes ----
PWard = cutree(tree = CAHDP, k=4) #Ward
PMax = cutree(tree = CAHDP2, k=4) #Sautmax
PMoy = cutree(tree = CAHDP3, k=4) #Sautmoy

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

indice_de_rand <- rand_index(PWard, PMax)
cat("Indice de Rand :", indice_de_rand, "\n")


