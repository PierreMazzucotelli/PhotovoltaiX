import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# Physical constants
q = 1.602176634e-19      # Elementary charge (C)
k_B = 1.380649e-23       # Boltzmann constant (J/K)
epsilon_0 = 8.854187817e-12  # Vacuum permittivity (F/m)
epsilon_si = 11.7 * epsilon_0  # Permittivity of silicon (F/m)

# Material and doping parameters
T = 300.0  # Temperature in Kelvin

#paramètres de la jonction (homojonciton)
Eg = 1.1*1.6*10**(-19)  # Bandgap energy of silicon (eV)
NA=1e12  # Acceptor concentration (m^-3)
ND=1e18  # Donor concentration (m^-3)
Nc=1e20  # Effective density of states in the conduction band (m^-3)
Nv=1e20  # Effective density of states in the valence band (m^-3)


def lire_donnees_fichier(file_path):
    """Lit les données de tension et de courant à partir d'un fichier texte."""
    voltages = []
    currents = []

    with open(file_path, 'r') as file:
        # Lire et ignorer la première ligne (en-tête)
        next(file)

        # Parcourir les lignes du fichier
        for line in file:
            # Diviser la ligne en colonnes
            columns = line.split()

            # Vérifier si la ligne contient des données valides
            if len(columns) >= 2:
                # Extraire la tension et le courant
                voltage = float(columns[0].replace(',', '.'))
                current = float(columns[1].replace(',', '.'))

                # Ajouter les valeurs aux listes
                voltages.append(voltage)
                currents.append(current)

    return voltages, currents

def potential_barrier(NA, ND, nC, nV, Eg):
    return (Eg-k_B*T*np.log((nC*nV)/(NA*ND)))/q

def depletion_widths(V):
    """
    Calcule les largeurs de déplétion modifiées par la tension appliquée V.
    Pour une tension appliquée, le potentiel effectif devient Veff = Vbi - V.
    """
    Veff = np.abs(potential_barrier(NA, ND, Nc, Nv, Eg) - V)
    Wn = np.sqrt((Veff * epsilon_0 * epsilon_si * 2 * NA) / (q * ND * (ND + NA)))
    Wp = np.sqrt((Veff * epsilon_0 * epsilon_si * 2 * ND) / (q * NA * (ND + NA)))
    return Wn, Wp

def potential_profile(x, V, beta=2.0):
    """
   Calcule le profil continu de potentiel le long de x, uniquement dans la région de déplétion [-Wp, Wn].
    Pour x entre -Wp et 0 (côté p) : profil quadratique décroissant de Vbi à une valeur commune à x = 0.
    Pour x entre 0 et Wn (côté n) : profil quadratique montant de 0 jusqu'à cette même valeur à x = 0.
    """
    Vbi = potential_barrier(NA, ND, Nc, Nv, Eg) - V
    Wn, Wp = depletion_widths(V)
    # Centre de la transition
    xm = (-Wp + Wn) / 2.0
    phi = np.where(x <= -Wp, Vbi,
                   np.where(x >= Wn, 0, Vbi * 0.5 * (1 - np.tanh(beta * (x - xm)))))
    return phi

def plot_energy_profile(num=500):
    """
    Trace l'évolution du profil de potentiel (énergie) en fonction de la distance x
    pour une tension appliquée V et affiche les bornes de déplétion (-Wp et Wn).
    """
    Vlist = [-1, 0, 1]
    Wn, Wp = depletion_widths(V)
    x = np.linspace(-Wp - 3, Wn + 3, num)
    phi = potential_profile(x, V)

    plt.figure()
    plt.plot(x, phi, label=f"V = {V:.2f} V")

    for V in Vlist:
        # Calcul du potentiel de jonction Vbi (pour x = -Wp) et de 0 à x = Wn
        Vbi = potential_barrier(NA, ND, Nc, Nv, Eg) - V

        # Ajout des marqueurs aux positions -Wp et Wn
        plt.plot([-Wp, Wn], [Vbi, 0], 'kx', markersize=8, label="Bornes de déplétion")

        # Ajout des lignes verticales indiquant -Wp et Wn
        plt.axvline(x=-Wp, linestyle="--", color="k", alpha=0.5)
        plt.axvline(x=Wn, linestyle="--", color="k", alpha=0.5)

        plt.xlabel("Distance x (m)")
        plt.ylabel("Potentiel (V)")
        plt.title("Comparaison des profils d'énergie de la jonction p/n")
        plt.grid(True)
        plt.legend()
        plt.show()


Dn = 1.0e-3  # Diffusion coefficient for electrons (m^2/s)
ni = np.exp(-Eg/(2*k_B*T))*(Nc*Nv)**0.5  # Intrinsic carrier concentration (m^-3)
Dp = 1.0e-3  # Diffusion coefficient for holes (m^2/s)
taun = 1.0e-3  # Electron lifetime (s)
taup = 1.0e-3  # Hole lifetime (s)
# Diffusion lengths
Ln = np.sqrt(Dn * taun)  # Electron diffusion length (m)
Lp = np.sqrt(Dp * taup)  # Hole diffusion length (m)

dp = 200*10**-6 # epaisseur de p
dn = 1*10**-6 # epaisseur de n

def coth(x):
    return 1/np.tanh(x)

def js(V):
    coef1 =q*Dn*ni**2/(Ln*NA)
    coef2 = q*Dp*ni**2/(Lp*ND)
    Wn, Wp = depletion_widths(V)
    js = coef1*coth((dp-Wp)/Ln) + coef2*coth((dn-Wn)/Lp)
    return js

def j_th(V, Rp, Rs):
    """
    Calcule le courant de la jonction p-n en fonction de la tension appliquée V.
    """
    j1 = js(V)
    j0 = j1 * (np.exp(q * V / (k_B * T)) - 1)
    return j0*10**(-4)-jcc + 0*V/Rp # Convertit en A/cm^2 et Rp en Ohm.cm2

def f(j, V, Rs):
    return j_th(V - Rs * j, Rp, Rs) - j * (1 + Rs / Rp) + V / Rp

def plot_current_voltage_custom(Rs_value):
    """
    Trace la caractéristique I-V pour une valeur de Rs donnée.
    """
    V0 = np.linspace(-1, 1, 1000)
    j0 = []
    for V in V0:
        def f2(j):
            return f(j, V, Rs_value)
        j0.append(brentq(f2, -2/Rp, 100))
    j0 = np.array(j0)
    plt.figure()
    plt.plot(V0, j0, label=f"Rs = {Rs_value}")
    plt.xlabel("Tension (V)")
    plt.ylabel("Courant (A/cm^2)")
    plt.title("Caractéristique courant-tension de la jonction p-n")
    plt.grid(True)
    plt.legend()
    plt.show()

def plot_current_voltage_log():
    """
    Trace la caractéristique courant-tension de la jonction p-n sur une échelle logarithmique.
    """
    V0 = np.linspace(-1, 1, 1000)
    j0 = []
    for V in V0 :
        def f2(j) :
            return f(j, V, Rs)
        j0.append(brentq(f2, -2/Rp, 100))

    j0 = np.array(j0)

    plt.plot(V0, np.abs(j0), label="Caractéristique I-V (échelle log)")
    plt.scatter(voltages, np.abs(currents), marker='+', linestyle='', label="HJ exp")
    
"""
# Obscurité
file_path = "2505071 HJ 4 cm2_avec_masque0.txt"
voltages, currents = lire_donnees_fichier(file_path)
jcc = 0.00016 # Saturation current density (A/m^2)
Rp = 6*10**1     # Ohm*cm2
Rs = 2.5
plot_current_voltage_log()

# sous lumière
file_path = "2505071 HJ 4 cm2_avec_masque2.txt"
voltages, currents = lire_donnees_fichier(file_path)
jcc = 0.0228 # Saturation current density (A/m^2)
Rp = 4.7*10**1     # Ohm*cm2
Rs = 2.8
plot_current_voltage_log()

plt.xlabel("Tension (V)")
plt.ylabel("log(Courant)")
plt.yscale('log')
plt.title("Caractéristique courant-tension de la jonction p-n (échelle log)")
plt.grid(True)
plt.legend()
plt.show()
"""""
'''
def plot_multiple_iv_curves(rs_values):
    """
    Affiche plusieurs courbes I-V théoriques pour différentes valeurs de Rs.
    
    Parameters:
        rs_values (list): Liste des valeurs de Rs à tester.
    """
    V0 = np.linspace(-1, 1, 1000)
    plt.figure()
    for Rs_temp in rs_values:
        j_curve = []
        for V in V0:
            def f2(j):
                return f(j, V, Rs_temp)
            j_curve.append(brentq(f2, -2/Rp, 100))
        j_curve = np.array(j_curve)
        plt.plot(V0, np.abs(j_curve), label=f"Rs = {Rs_temp}")
    plt.xlabel("Tension (V)")
    plt.ylabel("Courant (A/cm^2)")
    plt.title("Caractéristiques I-V pour différentes valeurs de Rs à Rp = 100 Ohm.cm2 et aucun photocourant")
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.show()

# Exemple d'appel avec différentes valeurs de Rs
Rp = 10000
jcc = 0
rs_list = [0.1, 1.0,10, 100, 1000]
plot_multiple_iv_curves(rs_list)'''
'''
def plot_multiple_iv_curves_Rp(rp_values, Rs_value):
    """
    Affiche plusieurs courbes I-V théoriques pour différentes valeurs de Rp
    avec Rs et jcc constants.

    Parameters:
        rp_values (list): Liste des valeurs de Rp à tester.
        Rs_value (float): Valeur constante de Rs.
    """
    V0 = np.linspace(-1, 1, 1000)
    plt.figure()
    for Rp_val in rp_values:
        j_curve = []
        for V in V0:
            def f2(j):
                # On utilise ici Rp_val pour la courbe en cours et Rs_value constant.
                return j_th(V - Rs_value * j, Rp_val, Rs_value) - j*(1 + Rs_value/Rp_val) + V/Rp_val
            j_curve.append(brentq(f2, -2/Rp_val, 500000))
        j_curve = np.array(j_curve)
        plt.plot(V0, np.abs(j_curve), label=f"Rp = {Rp_val}")
    plt.xlabel("Tension (V)")
    plt.ylabel("Courant (A/cm^2)")
    plt.title("Caractéristiques I-V pour différentes valeurs de Rp\navec Rs nul et aucun photocourant")
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.show()

# Exemple d'appel avec différentes valeurs de Rp et Rs constant
Rs_value = 0.000001   # valeur fixe de Rs
jcc = 0          # jcc constant
rp_list = [10, 100, 1000, 10000, 100000]  # par exemple, différentes valeurs de Rp
plot_multiple_iv_curves_Rp(rp_list, Rs_value)'''

def j_th_modified(V, Rp, Rs, jcc_val):
    """
    Calcule le courant de la jonction p-n en fonction de la tension V,
    avec une valeur de jcc (courant de saturation) donnée.
    Retourne le courant en A/cm^2.
    """
    j1 = js(V)
    j0 = j1 * (np.exp(q * V / (k_B * T)) - 1)
    return j0 * 1e-4 - jcc_val + 0 * V / Rp

def f_jcc(j, V, Rs, Rp, jcc_val):
    """
    Fonction intermédiaire intégrant la chute de tension,
    pour une valeur de jcc constante.
    """
    return j_th_modified(V - Rs * j, Rp, Rs, jcc_val) - j * (1 + Rs / Rp) + V / Rp

def plot_multiple_iv_curves_jcc(jcc_values, Rp_const, Rs_const):
    """
    Affiche plusieurs courbes I-V théoriques pour différentes valeurs de jcc,
    avec Rp et Rs constants.
    
    Parameters:
        jcc_values (list): Liste des valeurs de jcc à tester.
        Rp_const (float): Valeur constante de Rp (Ohm.cm2).
        Rs_const (float): Valeur constante de Rs.
    """
    V0 = np.linspace(-1, 1, 3000)
    plt.figure()
    for jcc_val in jcc_values:
        j_curve = []
        for V in V0:
            def f2(j):
                return f_jcc(j, V, Rs_const, Rp_const, jcc_val)
            j_curve.append(brentq(f2, -2000/Rp_const, 50000))
        j_curve = np.array(j_curve)
        plt.plot(V0, np.abs(j_curve), label=f"jcc = {jcc_val} A/cm^2")
    plt.xlabel("Tension (V)")
    plt.ylabel("Courant (A/cm^2)")
    plt.title("Caractéristiques I-V pour différentes valeurs de jcc\navec Rp=1000 Ohm.cm2 et Rs=0.0001 Ohm.cm2 constants")
    plt.yscale("log")
    plt.grid(True)
    plt.legend()
    plt.show()

# Exemple d'appel avec différentes valeurs de jcc
Rp = 1000      # Valeur constante de Rp (Ohm.cm2)
Rs = 0.0001        # Valeur constante de Rs
jcc_values = [0, 0.00001, 0.0001,0.001, 0.01, 0.1, 1]  # Exemple de valeurs de jcc (A/m^2)
plot_multiple_iv_curves_jcc(jcc_values, Rp, Rs)