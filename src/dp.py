from math import sqrt
from math import floor
import numpy as np
import re

## \mainpage Projet court de M2 Bioinformatique
# \section intro_sec Introduction
# Réalisation d'un programme reprenant la méthode décrite dans l'article 3
#basé sur la double programmation dynamique (pour plus d’information 4). Le
#threading (enfilage) [1,2,3] est une stratégie pour rechercher des séquences
#compatibles avec une structure. Seul les carbones α de la protéine seront
#considérés. Vous utiliserez les potentiels statistiques DOPE [5].
# \section ressources Ressources
# 1) Jones, D.T., Taylor, W.R. & Thornton, J.M. (1992) A new approach to protein fold recognition. Nature. 358, 86-89.\n
# 2) Jones, D.T., Miller, R.T. & Thornton, J.M. (1995) Successful protein fold recognition by optimal sequence threading validated by
#rigorous blind testing. Proteins. 23, 387-397.\n
# 3) Jones, D.T. (1998) THREADER : Protein Sequence Threading by Double Dynamic Programming. (in) Computational Methods in
#Molecular Biology. Steven Salzberg, David Searls, and Simon Kasif, Eds. Elsevier Science. Chapter 13.\n
# 4) Protein Structure Comparison Using SAP - Springer\n
# 5) http://www.dsimb.inserm.fr/~gelly/doc/dope.par

class CarbonAlpha:
    """Classe représentant un carbone alpha de manière géométrique
    (coordonnées dans le plan et facteurs physico-chimiques."""
    def __init__(self, name, x, y, z):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
    def euclidean_dist(self, ca):
        """On passe a la méthode notre objet et un autre objet de la classe carbon_alpha."""
        return sqrt((self.x - ca.x)**2 + (self.y - ca.y)**2 + (self.z- ca.z)**2)
def pdb_parser(path_to_file):
    """Fonction qui va parser un fichier pdb afin d'extraire les positions
    x,y,et z de chaque carbone alpha. Renvoie une liste d'objets carbon_alpha"""
    with open(path_to_file, "r") as filin:
        liste_ca = [CarbonAlpha(line[17:20].strip(),
                                float(line[30:38].strip()),
                                float(line[38:46].strip()),
                                float(line[46:54].strip()))
                    for line in filin if line[12:16].strip() == "CA"]
    return liste_ca

def fasta_parser(path_to_file):
    """Fonction qui va parcourir un fichier fasta pour extraire tous les acides aminés.
    Elle renvoie une chaine de caractère."""
    with open(path_to_file, "r") as filin:
        seq = ""
        for line in filin:
            if line[0] == ">":
                continue
            else:
                seq += line[:-1]
    return seq

def dope_limitator(dopin, dopout):
    """N'extrait du fichier dope initial, que les lignes décrivant des
    interactions entre carbones alphas. Elle crée un fichier dope
    condensé ne contenant que ces interactions."""
    with open(dopin, "r") as filin, open(dopout, "w") as filout:
        for line in filin:
            name = re.findall(catch_name, line)
            if name[1] == 'CA':
                if name[3] == 'CA':
                    filout.write("{}".format(line))

def dope_parser(dope):
    """Fonction qui parse le fichier dope en input et renvoie une liste de dictionnaire."""
    with open(dope, "r") as filin:
        liste_dope = []
        for line in filin:
            d = {
                "AA_1": re.findall(catch_name, line)[0],
                "AA_2": re.findall(catch_name, line)[2],
                "Energy": re.findall(catch_number, line)
            }
            liste_dope.append(d)
    return liste_dope

def dist_to_dope(dist, name_1, name_2, liste):
    """Fais la conversion entre une distance entre deux résidus, et un
     potentiel DOPE dans la structure de donnée 'liste'.
    Nécessite le nom des deux résidus dont la distance a été évaluée."""
    for d in liste:
        if d["AA_1"] == name_1:
            if d["AA_2"] == name_2:
                if dist != 0:
                    return float(d["Energy"][floor(dist*4) - 1])
            return float(d["Energy"][0])

def min_finder(matrix):
    """Renvoie le minimum d'une matrice."""
    liste = list(matrix.reshape(matrix.size))
    return min(liste)
#==================================================

#Dans un premier temps il s'agit de calculer une matrice de distance,
#entre tous les acides aminés de la protéine, deux à deux.
NB_AA = 25
NB_CA = 25
print("Veuillez entrer le chemin vers un fichier pdb: ")
PATH_PDB = input()
print("Veuillez entrer le chemin vers un fichier fasta: ")
PATH_FASTA = input()
print("Veuillez entrer le chemin vers un fichier dope: ")
IN_DOPE = input()
print("Veuillez la destination du nouveau fichier dope: ")
OUT_DOPE = input()

liste_ca = pdb_parser(PATH_PDB)
seq = fasta_parser(PATH_FASTA)
dist_matrix = np.zeros((NB_AA, NB_CA))

for i in range(0, dist_matrix.shape[0]):
    for j in range(0, dist_matrix.shape[1]):
        dist_matrix[i, j] = liste_ca[i].euclidean_dist(liste_ca[j])
#Dans un second temps on va devoir relier l'information de distance avec le potentiel DOPE
catch_number = re.compile("-{0,}[0-9]{1,}.[0-9]{1,}")
catch_name = re.compile("[A-Z]{1,3}")
dope_limitator(IN_DOPE, OUT_DOPE)
liste_dope = dope_parser(OUT_DOPE)
abreviation = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET",
    "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR",
    "V": "VAL", "W": "TRP", "Y": "TYR"
}

prime_matrix = np.zeros((NB_AA, NB_CA, dist_matrix.shape[0] + 1, dist_matrix.shape[1] + 1))
#on stocke dans chaque case de cette matrice, une low_matrix

#construction des low_lvl_matrix:

for i in range(0, NB_AA): #résidu
    for j in range(0, NB_CA): #ca
        #pour chaque couple on doit créer une low level matrix
        low_level_matrix = np.zeros((NB_AA + 1, NB_CA + 1))
        for k in range(1, low_level_matrix.shape[0]):
            for l in range(1, low_level_matrix.shape[1]):
                choice = min(low_level_matrix[k-1, l-1], low_level_matrix[k-1, l],
                             low_level_matrix[k, l-1])
                distance = liste_ca[j].euclidean_dist(liste_ca[l-1])
                try:
                    energy = dist_to_dope(distance, abreviation[seq[i]],
                                          abreviation[seq[k-1]], liste_dope)
                except IndexError:
                    energy = 0
                low_level_matrix[k, l] = choice + energy
        prime_matrix[i, j, :, :] = low_level_matrix

#remplissage de la high_lvl_matrix:

high_level_matrix = np.zeros((NB_AA + 1, NB_CA + 1))

for i in range(1, high_level_matrix.shape[0]):
    for j in range(1, high_level_matrix.shape[1]):
        choice = min([high_level_matrix[i-1, j-1], high_level_matrix[i-1, j],
                      high_level_matrix[i, j-1]])
        high_level_matrix[i, j] = choice + min_finder(prime_matrix[i-1, j-1, :, :])
#resultat final:
print("Score d'adéquation entre la séquence '{}' et la structure '{}' = {:.2f}"
      .format(PATH_FASTA, PATH_PDB, min_finder(high_level_matrix)))
