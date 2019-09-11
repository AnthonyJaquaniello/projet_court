from math import sqrt
import numpy as np
import re

class carbon_alpha:
    """
        Classe représentant un carbone alpha de manière géométrique 
        (coordonnées dans le plan et facteurs physico-chimiques.
    """
    def __init__(self, name, x, y, z):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        
def pdb_parser(path_to_file):
    """
        Fonction qui va parser un fichier pdb afin d'extraire les positions
        x,y,et z de chaque carbone alpha. Renvoie une liste d'objets carbon_alpha
    """
    with open(path_to_file,"r") as filin:
        liste_ca = [carbon_alpha(line[17:20].strip(),float(line[30:38].strip()),float(line[38:46].strip()), 
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
    
def dope_limitator(dope):
    """
        N'extrait du fichier dope initial, que les lignes décrivant des interactions entre carbones alphas.
        Elle crée un fichier doep condensé ne contenant que ces interactions.
    """
    with open(dope,"r") as filin, open("dope_limited.txt","w") as filout:
        for line in filin:
            name = re.findall(catch_name,line)
            if name[1] == 'CA':
                if name[3] == 'CA':
                    filout.write("{}\n".format(line))
                   
def dope_parser(dope):
    """
        Fonction qui parse le fichier dope en input et renvoie une liste de dictionnaire.
    """
    with open(dope,"r") as filin:
        liste_dope= []
        for line in filin:
            d = {"AA_1":re.findall(catch_name,line)[0],"AA_2":re.findall(catch_name,line)[2],
                 "Energy":re.findall(catch_number,line)}
            liste_dope.append(d)
    return liste_dope

#Dans un premier temps il s'agit de calculer une matrice de distance entre tous les acides aminés de la protéine, deux à deux.

liste_ca = pdb_parser("../data/pdb/2ai9.pdb")
seq = fasta_parser("../data/fasta/6p4y.fasta.txt")

dist_matrix = np.ones((len(liste_ca),len(liste_ca)))

for i in range(0,dist_matrix.shape[0]):
    for j in range(0,dist_matrix.shape[1]):
        dist_matrix[i,j] = euclidean_dist(liste_ca[i],liste_ca[j])
    
#Dans un second temps on va devoir relier l'information de distance avec le potentiel DOPE
catch_number = re.compile("-{0,}[0-9]{1,}.[0-9]{1,}")
catch_name = re.compile("[A-Z]{1,3}")

dope_limitator("../data/dope.par.txt")



