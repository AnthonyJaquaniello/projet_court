## Information préalable:

* La complexité de l'algorithme de double programmation dynamique est n²*m² où n est la longueur de la séquence fasta et m le nombre de carbones alpha du template. Afin d'éviter une explosion de la complexité, éviter de dépasser n = 100 et m = 100.
* Le réglage de n et m est possible dans le code (NB_AA, et NB_CA)
* Pour fonctionner, le programme a besoin d'un fichier de potentiel statistique (DOPE) (cf ci dessous)

## Lancement:

python3 dp.py <\arg1> <\arg2> <\arg3> <\arg4>

## Arguments:

1. Chemin vers le fichier .pdb
2. Chemin vers le fichier .fasta.txt
3. Chemin vers le fichier .par (fichier dope)
4. Chemin vers le fichier dope de sortie (nouveau fichier qui sera crée à la volée).

## Exemple:

**python3 dp.py ../data/pdb/1lry.pdb ../data/fasta/1lry.fasta.txt ../data/dope.par.txt ../data/dope_limited.txt**

### Jeux de donnés:

Pour effectuer un test du programme nous avons mit à disposition quelques fichiers fasta et pdb sur le dépôt github associé:
**https://github.com/AnthonyJaquaniello/projet_court**

### Signification du résultat:

Le programme renvoie une valeur négative qui est un potentiel énergétique issu d'une double programmation dynamique. Cette valeur
indique l'adéquation de la structure à la séquence proposée. Plus cette valeur est négative, plus l'adéquation est forte
(énergie plus basse donc stabilité du modèle plus grand).

## Complément:

Documentation disponible dans le dossier doc, au format html.


