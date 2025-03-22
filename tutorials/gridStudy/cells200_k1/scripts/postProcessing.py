import numpy as np
import matplotlib.pyplot as plt

def read_file_to_numpy_array(file_path, delimiter=None):
    """
    Liest Werte aus einer Datei ein und speichert sie in einem NumPy-Array.

    :param file_path: Pfad zur Eingabedatei
    :param delimiter: Trennzeichen zwischen den Werten (z. B. ',' für CSV)
                      Wenn None, wird automatisch versucht, die Datei zu parsen.
    :return: NumPy-Array mit den eingelesenen Werten
    """
    try:
        # Laden der Datei als NumPy-Array
        data = np.loadtxt(file_path, delimiter=delimiter)
        return data
    except Exception as e:
        print(f"Fehler beim Einlesen der Datei: {e}")
        return None

# Beispiel: Dateipfad eingeben
file_path = "./postProcessing/sample/500/x4_U.xy"  # Dateiname anpassen

# Datei einlesen und in NumPy-Array speichern
u = read_file_to_numpy_array(file_path)
u_x = u[:, 0:2]
x = u_x[:,0]
y = u_x[:,1]

u_x_pol_List = []
for i in range(u[:,1].size):
    if u[i,1]>1e-3:
        u_x_pol_List.append(u[i, 0:2])
u_x_pol = np.array(u_x_pol_List)
x_pol = u_x_pol[:,0]
y_pol = u_x_pol[:,1]

if u is not None:
    print("Werte eingelesen.")
    #print(x, y)
else:
    print("Die Datei konnte nicht gelesen werden.")

# Interpolation: Polynom 2. Grades (Parabel) finden
coefficients = np.polyfit(x_pol, y_pol, deg=2)  # Grad 2 für eine Parabel
polynomial = np.poly1d(coefficients)

# Nullstellen berechnen
roots = np.roots(coefficients)

# Werte für die Parabel plotten
x_new = np.linspace(roots[0], roots[1], 500)  # Erstelle fein aufgelöste x-Werte
y_new = polynomial(x_new)  # Berechne die y-Werte der Parabel

# Plot der Punkte und der interpolierten Parabel
plt.scatter(x, y, color='red', marker = ".", linewidth=0.5, label='Werte der Zellmittelpunkte')  # Originalpunkte
plt.plot(x_new, y_new, color='blue', linewidth=1, label='Interpolierte Parabel')  # Parabel
k = 1

# Analytische Lösung
analyticSolution_pol = 0.015*(x_pol-0.1-k*0.02)/0.4*(2-(x_pol-0.1-k*0.02)/0.4)
# L2-Fehler berechnen
l2_error = np.sqrt(np.sum((analyticSolution_pol - polynomial(x_pol)) ** 2))
print("L2-Error: ", l2_error)

plt.plot(x_new, 0.015*(x_new-0.1-k*0.02)/0.4*(2-(x_new-0.1-k*0.02)/0.4), color='green', linewidth=1, linestyle = "dotted", label='Analytische Lösung')

# Nullstellen in den Plot einfügen
for root in roots:
    plt.scatter(root, 0, color='green', marker = "x", zorder=5)  # Nullstelle bei y=0
    print(root)

# Achsenbeschriftung und Titel
plt.xlabel('x')
plt.ylabel('y')
#plt.xlim([0.05,0.15])
#plt.ylim([0.01,0.016])
plt.title('Interpolation einer Parabel mit Nullstellen')
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Nullachse
plt.legend()
plt.grid()

# Speichere das Bild ab
output_file = "images/interpolierte_parabel.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"Das Bild wurde gespeichert als: {output_file}")
