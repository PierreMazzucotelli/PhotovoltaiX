import numpy as np
import matplotlib.pyplot as plt

# Physical constants
q = 1.602176634e-19      # Elementary charge (C)
k_B = 1.380649e-23       # Boltzmann constant (J/K)
epsilon_0 = 8.854187817e-12  # Vacuum permittivity (F/m)
epsilon_si = 11.7 * epsilon_0  # Permittivity of silicon (F/m)

# Material and doping parameters
T = 300.0  # Temperature in Kelvin
Na = 1e24  # Acceptor concentration (m^-3)
Nd = 1e22  # Donor concentration (m^-3)
ni = 1.5e16  # Intrinsic carrier concentration (m^-3) at 300K

#paramètres de la jonction (homojonciton)
Eg = 1.12  # Bandgap energy of silicon (eV)
NA=1e16  # Acceptor concentration (m^-3)
ND=1e16  # Donor concentration (m^-3)
Nc=1e19  # Effective density of states in the conduction band (m^-3)
Nv=1e19  # Effective density of states in the valence band (m^-3)

# On calcule deltaV
def potential_barrier(NA, ND, nC, nV, Eg):
    return (Eg-k_B*T*np.log((nC*nV)/(NA*ND)))/q

def depletion_widths(V):
    """
    Calcule les largeurs de déplétion modifiées par la tension appliquée V.
    Pour une tension appliquée, le potentiel effectif devient Veff = Vbi - V.
    """
    Veff = potential_barrier(NA, ND, Nc, Nv, Eg) - V
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

if __name__ == "__main__":
    plot_energy_profile()