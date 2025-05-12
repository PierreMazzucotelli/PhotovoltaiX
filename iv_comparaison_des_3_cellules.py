import matplotlib.pyplot as plt
import numpy as np

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

def MPP(V, I) :
    power = 0
    imax = 0
    for i in range(len(V)) :
        if(abs(V[i]*I[i])>power and I[i]<0 and V[i]>0) :
            imax = i
            power = abs(V[i]*I[i])
    return 1000*power, imax     # power en mW/cm^2


file_path = "2505071 HJ 4 cm2_avec_masque.txt"
voltages, currents = lire_donnees_fichier(file_path)

# Tracer le courant en fonction de la tension
plt.plot(voltages, currents, marker='+', linestyle='', label="Hétérojonction 2505071 - FF = 28.6")
power, imax = MPP(voltages, currents)
plt.scatter(voltages[imax], currents[imax], s = 100)
plt.text(voltages[imax], currents[imax]-0.002, "   Pmax = " + f"{power:.2g}" + " mW/cm²", fontsize=10)

file_path = "tandem_photo1.txt"
voltages, currents = lire_donnees_fichier(file_path)

# Tracer le courant en fonction de la tension
plt.plot(voltages, currents, marker='+', linestyle='', label="Tandem - FF = 65.7")
power, imax = MPP(voltages, currents)
plt.scatter(voltages[imax], currents[imax], s = 100)
plt.text(voltages[imax], currents[imax]-0.002, "   Pmax = " + f"{power:.2g}" + " mW/cm²", fontsize=10)

file_path = "2505022-2-sous iradiation1000.txt"
voltages, currents = lire_donnees_fichier(file_path)

# Tracer le courant en fonction de la tension
plt.plot(voltages, currents, marker='+', linestyle='', label="PIN 2505022 - FF = 59.5")
power, imax = MPP(voltages, currents)
plt.scatter(voltages[imax], currents[imax], s = 100)
plt.text(voltages[imax], currents[imax]-0.002, "   Pmax = " + f"{power:.2g}" + " mW/cm²", fontsize=10)

# Tracer le courant en fonction de la tension
plt.title('Courbes I(V) de 3 cellules différentes')
plt.xlabel('Tension (V)')
plt.ylabel('Courant (A)')
#plt.yscale('log')
plt.grid(True)
plt.legend()
plt.show()
