# Intégration du cycle cellulaire dans un modèle à base de contrainte

  L'objectif de ce stage est de faire une modélisation hybride dans le but de pouvoir expliquer les comportements différents observés entre cellules saines et  cancéreuses. Pour cela, des méthodes à base de contraintes itératives (Dynamic Flux Balance Analysis, DFBA) permettent de calculer des voies métaboliques optimales à chaque étape du cycle et ainsi obtenir un modèle dynamique. Dans ce rapport, deux objectifs sont abordés : Définir des paramètres pour retrouver les comportements de la cellule pour chacune de ses phases de façon statique, puis prendre en main les méthodes python DFBA. Ainsi, en optimisant la biomasse pour maximiser les besoins métaboliques de la cellule et en simulant quelques régulations enzymatiques, les résultats obtenus concordent avec les résultats expérimentaux connus [1](#Référence) : une production d’énergie, d’acides aminés et nucléotides passant par la voie de la glycolyse principalement pour la phase G1, une production de nucléotides par la voie des pentoses phosphate pour la phase et S et la production de lipides pour la phase G2. Ce sont ces paramètres qui pourront ensuite être utilisés pour définir le modèle dynamique avec la méthode DFBA.

## Pré-requis

Librairies python:

* Cobrapy (0.24)
*	DFBA (Tourigny) [2](#Référence)
*	Matplotlib
*	Numpy
*	Pandas
*	Seaborn

## Description du contenu utilisé

| Fichier| Contenu |
|-----------|-----------|
| `./small_model/model_mario2013.xml` | Modèle CHO [3](Référence) au format sbml |
| `./small_model/model_mario2013G1.xml ` | Modèle CHO [3](Référence) au format sbml modifié pour la phase G1 du cycle cellulaire |
| `./small_model/model_mario2013S.xml ` | Modèle CHO [3](Référence) au format sbml modifié pour la phase S du cycle cellulaire |
| `./small_model/model_mario2013G2.xml` | Modèle CHO [3](Référence) au format sbml modifié pour la phase G2 du cycle cellulaire |
| `HMRdatabase2_00.xml` | Modèle HMR [4](#Référence) au format sbml |
| `fba_small_new.ipynb` | Notebook réalisant des FBA sur le modèle CHO pour chaque phase du cylce cellulaire dans le but de valider le modèle |
| `fba_cellcycle2.ipynb` | Notebook réalisant des FBA sur le modèle HMR pour chaque phase du cylce cellulaire dans le but de valider le modèle |
| `dfba_small.py` | Scipt simulant un dFBA pour le modèle CHO sans condition particulière (pas intégré au cycle celulaire) |

## Références <a name="Références"></a>

[1]	Moulin C. Analyse des voies métaboliques au cours du cycle cellulaire : application au métabolisme du cancer. Bio-Informatique, Biologie Systémique [q-bio.QM]. Université Paris-Saclay, 2020. Français. ⟨NNT : 2020UPASG022⟩. ⟨tel-03164566⟩

[2]	Tourigny et al., (2020). dfba: Software for efficient simulation of dynamic flux-balance analysis models in Python. Journal of Open Source Software, 5(52), 2342, https://doi.org/10.21105/joss.02342

[3] Ghorbaniaghdam, A., Henry, O. & Jolicoeur, M. A kinetic-metabolic model based on cell energetic state: study of CHO cell behavior under Na-butyrate stimulation. Bioprocess Biosyst Eng 36, 469–487 (2013). https://doi.org/10.1007/s00449-012-0804-3

[4] Mardinoglu, A., Agren, R., Kampf, C. et al. Genome-scale metabolic modelling of hepatocytes reveals serine deficiency in patients with non-alcoholic fatty liver disease. Nat Commun 5, 3083 (2014). https://doi.org/10.1038/ncomms4083
