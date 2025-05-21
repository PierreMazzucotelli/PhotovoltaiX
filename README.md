# PhotovoltaiX

<img src="https://github.com/PierreMazzucotelli/PhotovoltaiX/blob/main/logo.png" alt="Logo du Projet" width="400" height="400">



Ce dépôt contient un module de simulation destiné à modéliser les caractéristiques électriques et les profils énergétiques d'une jonction p-n, tels qu'utilisés dans les dispositifs photovoltaïques. Il a été réalisé durant le MODAL de physique de la matière condensée en deuxième année à l'Ecole Polytechnique.

## 📌 Aperçu

Le fichier `Simulation Modal.py` propose diverses fonctions permettant de :

- Lire les données expérimentales (tension et courant) à partir de fichiers texte.
- Calculer les propriétés physiques de la jonction (barrière de potentiel, largeurs de déplétion, profil de potentiel).
- Déterminer le courant théorique de la jonction via la loi diode, en intégrant notamment l'impact du courant de saturation (`jcc`), de la résistance série (`Rs`) et de la résistance shunt (`Rp`).
- Tracer les caractéristiques courant-tension (I-V) en échelle linéaire ou logarithmique.
- Afficher plusieurs courbes I-V issues de la variation d'un paramètre (ex. `Rs`, `Rp` ou `jcc`) tout en maintenant les autres constantes.

---

## ⚙️ Fonctionnalités

- **Lecture et traitement de données**  
  La fonction `lire_donnees_fichier(file_path)` permet de charger les mesures expérimentales à partir d'un fichier texte.

- **Modélisation de la jonction p-n**  
  Des fonctions telles que `potential_barrier`, `depletion_widths` et `potential_profile` calculent les propriétés essentielles de la jonction en fonction des paramètres de dopage et de la tension appliquée.

- **Calcul et affichage de la caractéristique I-V**  
  Les fonctions `j_th` et `f` calculent le courant en intégrant la loi exponentielle de la jonction et les effets résistifs.  
  Des fonctions de tracé comme `plot_current_voltage_custom` et `plot_current_voltage_log` permettent d'afficher ces courbes.

- **Variations paramétriques**  
  Plusieurs fonctions (`plot_multiple_iv_curves`, `plot_multiple_iv_curves_Rp`, `plot_multiple_iv_curves_jcc`) permettent d'observer l'impact de la variation de `Rs`, `Rp` ou `jcc` sur la caractéristique I-V.

---
## Licence
Ce projet est sous licence MIT.

---

## 🧰 Prérequis

- Python 3.x  
- Les packages suivants sont nécessaires :

```bash
pip install numpy matplotlib scipy
