#set.seed(6)

T=300

# On va essayer de simuler des chaînes de longueur T de manière uniforme. On aura pas à coup sûr 
# des chaînes de longueur T. Plus T est grand moins cela est probable  
start_time <- Sys.time()

myope_algo_et_poids = function (T) 
{ #
    X_coord = cbind(0,0) # marche 
    bloque = 0
    direction <- seq(0,3) 
    direction <- cbind((direction %/% 2)*(1-2*(direction %% 2)), (1-direction %/% 2)*(1-2*(direction %% 2)))
    poids = rep(4,1)  # vecteur des poids à chaque étape, partant de (0,0) les 4 directions sont disponibles 
    t=1
    while ( t <= T & bloque!= 1 )
    # Tant que la chaîne ne revient pas vers un point déjà visité
    {   
        direction_dispo = matrix(nrow=0,ncol=2) # contient les directions qui sont faisables sans intersection
        # On cherche les directions faisables
        for (i in seq(1,4)){
        if (!match(TRUE, direction[i,1]+X_coord[t,1] == X_coord[,1] & direction[i,2]+X_coord[t,2] == X_coord[,2], nomatch=0))
            {direction_dispo <- rbind(direction_dispo, direction[i,])}
        }
        long =  length(direction_dispo)/2  # correspond à  longueur de direction_dispo
        if (long >=1  ){
            # Il y a au moins une direction faisable, on en choisit une uniformément
            poids <- rbind(poids, long) # le poids de l'étape t correspond au nombre de directions disponibles à cet étape
            indice = floor(runif(1,min=1,max=long+1))
            direction_choisie = direction_dispo[indice, ]
            X_coord <- rbind(X_coord, X_coord[t,] + direction_choisie)
            t = t+1
        }
        else{
            # Il n'y a pas de direction disponible/ faisable
            bloque = 1 
            t=t+1
        }
        
    }
    donnees = data.frame(cbind(X_coord,poids))
    colnames(donnees) = c('Abscisses','Ordonnées','poids')
    rownames(donnees) = seq(1,nrow(X_coord)) 
    return (donnees)
        
}

Y = myope_algo_et_poids(T)
print(Y)
nrow(Y)
plot(Y[1:nrow(Y),1:2], type='o')

print(Sys.time() - start_time)

N = 1000 # nombre de simulations
start_time <- Sys.time()
# pour des longueurs de 10, 20,50,80,100,200,100 on va déterminer le pourcentage de simulations pour lesquelles la chaîne 
# ne se coupe pas avant d'atteindre la longueur
liste = c(1,10,20,50,80,100,200,400)
#liste = seq(1,500)
proportion = rep(0, length(liste))
for (k in seq(1,length(liste))){
    t = liste[k]
    for ( n in seq(1,N))
    {
        if (nrow(myope_algo_et_poids(t)[1]) == t+1){
            proportion[k] = proportion[k] + 1}
    }
    proportion[k] = proportion[k]/N
}

print(Sys.time() - start_time)
plot(liste, proportion, type='o', xlab =  "Longueurs de la chaîne", ylab = "proportion" )

start_time <- Sys.time()

Importance_Sampling1 = function(T,N,l){
    # T durée de la chaîne, N nombre de simulations, l seuil , l<=T
     
 
    liste_poids = matrix( nrow=0 , ncol = 1) # les poids de ces chaînes
    indicatrice = matrix( nrow=0 , ncol = 1) # contiendra des 0 et 1: 1 si la chaîne est de longueur >= l, 0 sinon 
    
    for (nbre_chaine in seq(1,N)){
        resultat = myope_algo_et_poids(T)
        X = resultat[1:nrow(resultat),1:2] # chaîne simulée
        #print(nrow(X))
        #print(prod(resultat[,3]))
        # On affecte à la chaîne simulée comme poids le produit de ces poids si la chaîne est de longueur T
        # Sinon on lui affecte un poids nul
        liste_poids <- rbind(liste_poids, (nrow(X)!= T+1) * prod(resultat[,3])) 
        indicatrice <- rbind(indicatrice, sum(nrow(X) >= l))
        #print(indicatrice[nbre_chaine])
        #print(sum(nrow(X)>=l))   
    }
    #print(liste_poids)
    #s = 0
    #sp = 0 # somme des poids
    #poid = 0
    #for (i in seq(1,N)){
        #print(cbind(s,sp,poid))
        #poid = liste_poids[i]
        #s = s + indicatrice[i]*liste_poids[i]
        #sp = sp + liste_poids[i]
    #}
    valeur_estimee = sum( liste_poids * indicatrice) / sum(liste_poids)
    #print(sum(liste_poids))
    ecart_type = sqrt( N * sum( ( liste_poids / sum( liste_poids ) )^2 * (indicatrice - valeur_estimee)^2))
    
    borne_inf_intervalle_confiance = valeur_estimee - qnorm(0.975)  * sqrt(1/N) * ecart_type
    borne_sup_intervalle_confiance = valeur_estimee + qnorm(0.975) * sqrt(1/N) * ecart_type
    
    tableau = cbind(N, valeur_estimee, ecart_type, borne_inf_intervalle_confiance, borne_sup_intervalle_confiance)
    tableau = data.frame(tableau)
    colnames(tableau) = c('NbreSimulation', 'Proba estimée', 'Ecart_type', 'Borne Inf', "Borne Sup")
    return (tableau)
}

T = 200
seuil = 195
resultat1 = Importance_Sampling1(T,1000,seuil)
resultat1
print(Sys.time() - start_time)

start_time <- Sys.time()

Importance_Sampling2 = function(T, N, l, Nbre_iter_max){
    # T durée de la chaîne, N nombre de simulations, l seuil , l<=T, Nbre_iter_max nombre maximal d'itérations
    # On impose  ici la condition d'avoir N chaînes de longueur T
    n = 0
    liste_poids = matrix( nrow=0 , ncol = 1) # les poids de ces chaînes
    indicatrice = matrix( nrow=0 , ncol = 1) # contiendra des 0 et 1: 1 si la chaîne est de longueur >= l, 0 sinon 
    
    while (n < N){
        resultat = myope_algo_et_poids(T)
        X = resultat[1:nrow(resultat),1:2] # chaîne simulée
        #print(nrow(X))
        #print(prod(resultat[,3]))
        n_iter = 1
        while (nrow(X)!= T+1 & n_iter < Nbre_iter_max){ 
        # tant qu'on n'a pas une chaîne de longueur désirée ( T ici ) et qu'on n'a pas atteint 
        # le nombre max de simulations
        resultat = myope_algo_et_poids(T)
        X = resultat[1:nrow(resultat),1:2]
        n_iter = n_iter + 1
        } 
        
        if (nrow(X)!= T+1){
        # On affecte à la chaîne simulée comme poids le produit de ces poids
        liste_poids <- rbind(liste_poids, prod(resultat[,3])) 
        indicatrice <- rbind(indicatrice, sum(nrow(X) >= l))
        n = n + 1
        # print(cbind(sum(nrow(X) >= l),prod(resultat[,3]))) 
        }
        
        #print(indicatrice[nbre_chaine])
        #print(s)
    }
    #print(indicatrice)
    #print(cbind(sum( liste_poids * indicatrice),sum(liste_poids)))
    valeur_estimee = sum( liste_poids * indicatrice) / sum(liste_poids)
    # print(valeur_estimee)
    ecart_type = sqrt( N * sum( ( liste_poids * (indicatrice - valeur_estimee)/ sum( liste_poids ) )^2))
    
    borne_inf_intervalle_confiance = valeur_estimee - qnorm(0.975) * sqrt(1/N) * ecart_type
    borne_sup_intervalle_confiance = valeur_estimee + qnorm(0.975)  * sqrt(1/N) * ecart_type
    
    tableau = cbind(N, valeur_estimee, ecart_type, borne_inf_intervalle_confiance, borne_sup_intervalle_confiance)
    tableau = data.frame(tableau)
    colnames(tableau) = c('NbreSimulation', 'Proba estimée', 'Ecart_type', 'Borne Inf', "Borne Sup")
    return (tableau)
}

T = 80
N = 1000
seuil = 70
Nbre_iter_max = 50

resultat2 = Importance_Sampling2(T, N, seuil, Nbre_iter_max)
resultat2

print(Sys.time() - start_time)


