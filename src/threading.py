#!/usr/bin/env python3

from math import sqrt
import numpy as np

class carbon_alpha:
    """
        Classe représentant un carbone alpha de manière géométrique 
        (coordonnées dans le plan et facteurs physico-chimiques.
    """
    def __init__(self, num, x, y, z):
        self.num = num
        self.x = x
        self.y = y
        self.z = z
        
def pdb_parser(path_to_file):
    """
        Fonction qui va parser un fichier pdb afin d'extraire les positions
        x,y,et z de chaque carbone alpha. Renvoie une liste d'objets carbon_alpha
    """
    with open(path_to_file,"r") as filin:
        liste_ca = [carbon_alpha(line[6:11].strip(),float(line[30:38].strip()),float(line[38:46].strip()), 
        float(line[46:54].strip())) for line in filin if line[12:16].strip() == "CA"]
    return liste_ca
        
def fasta_parser(path_to_file):
    """
        Fonction qui va parcourir un fichier fasta pour extraire tous les acides aminés.
        Elle renvoie une chaine de caractère.
    """
    with open(path_to_file,"r") as filin:
        seq = ""
        for line in filin:
            if line[0] == ">":
                continue
            else:
                seq += line[:-1]
    return seq
    
def euclidean_dist(ca_1,ca_2):
    """
        Fonction qui renvoie la distance euclidienne entre deux objets carbon_alpha.
    """
    return  sqrt((ca_1.x - ca_2.x)**2 + (ca_1.y - ca_2.y)**2 + (ca_1.z- ca_2.z)**2)

#Dans un premier temps il s'agit de calculer une matrice de distance entre tous les acides aminés de la protéine, deux à deux.

liste_ca = pdb_parser("../data/pdb/2ai9.pdb")
seq = fasta_parser("../data/fasta/6p4y.fasta.txt")

dist_matrix = np.zeros((len(liste_ca),len(liste_ca)))

for i in dist_matrix:
    for j in dist_matrix:
        dist_matrix[i,j] = euclidean_dist(liste_ca[i],liste_ca[j])
    
print(dist_matrix)
