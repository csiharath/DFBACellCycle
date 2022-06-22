# Intégration du cycle cellulaire dans un modèle à base de contrainte
L'objectif de ce stage est de faire une modélisation hybride dans le but de pouvoir expliquer les comportements différents observés entre cellules saines et cancéreuses. Pour cela, des méthodes à base de contraintes itératives (Dynamic Flux Balance Analysis, DFBA) permettent de calculer des voies métaboliques optimales à chaque étape du cycle et ainsi obtenir un modèle dynamique. Dans ce rapport, deux objectifs sont abordés : Définir des paramètres pour retrouver les comportements de la cellule pour chacune de ses phases de façon statique, puis prendre en main les méthodes python DFBA. Ainsi, en optimisant la biomasse pour maximiser les besoins métaboliques de la cellule et en simulant quelques régulations enzymatiques, les résultats obtenus concordent avec les résultats expérimentaux connus [1](#Référence) : une production d’énergie, d’acides aminés et nucléotides passant par la voie de la glycolyse principalement pour la phase G1, une production de nucléotides par la voie des pentoses phosphate pour la phase et S et la production de lipides pour la phase G2. Ce sont ces paramètres qui pourront ensuite être utilisés pour définir le modèle dynamique avec la méthode DFBA.

## Pré-requis

Librairies python:

* Cobrapy (0.24)
*	DFBA (Tourigny) [2](#Référence)
*	Matplotlib
*	Numpy
*	Pandas
*	Seaborn

## Références <a name="Références"></a>

[1]	Moulin C. Analyse des voies métaboliques au cours du cycle cellulaire : application au métabolisme du cancer. Bio-Informatique, Biologie Systémique [q-bio.QM]. Université Paris-Saclay, 2020. Français. ⟨NNT : 2020UPASG022⟩. ⟨tel-03164566⟩

[2]	Tourigny et al., (2020). dfba: Software for efficient simulation of dynamic flux-balance analysis models in Python. Journal of Open Source Software, 5(52), 2342, https://doi.org/10.21105/joss.02342
