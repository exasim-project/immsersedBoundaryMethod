# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 13:32:21 2025

@author: Jannik Weber (ms8248)
         Karlsruher Insitut für Technologie (KIT)
         Institut für Strömungsmechanik (ISTM)
"""

import numpy as np

#case settings
DATA_DIR = 'postProcessing/controlVolumeSampling/2986/'
rho = 1
lines = ['frontLine', 'backLine', 'upperLine', 'lowerLine']
fields = ["p", "U", "grad(U)"]
lineNormals = [np.array([-1, 0, 0]), np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, -1, 0])]
angle = 0 #degree
Re = 8000
mu = rho/Re
Area = 1
dynamicPressure = 0.5*rho #1/2*rho(1)*u(1)^2
pointsPerLenght = 100
deltaS = 1/pointsPerLenght

#Functions:
def getForce(lineName, normal, fieldName):
    file = lineName + '_' + fieldName + '.xy'
    F = 0
    m = 0
    with open(DATA_DIR + file, "r") as f:

        for line in f:

            line = line.strip()
            if not line:
                continue  # leere Zeilen ignorieren
            
            # Zeile in Werte zerlegen (Leerzeichen getrennt)
            values = line.split()
            # in Float konvertieren
            values = [float(v) for v in values]

            if fieldName == 'p':
                pValue = np.array(values[1])
                F += pValue*normal*deltaS

            if fieldName == 'U':
                UValue = np.array((values[1], values[2], values[3]))
                F += rho*UValue*np.dot(UValue, normal)*deltaS
                m += rho*np.dot(UValue, normal)*deltaS

            if fieldName == 'grad(U)':
                gradUValue = np.array((values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9]))
                gradUValue = np.reshape(gradUValue, (3,3))
                shearStressTensor = mu*(gradUValue + np.transpose(gradUValue))
                F -= np.dot(shearStressTensor, normal)*deltaS

    return F, m


def main():
    F = 0
    m = 0
    for i in range(len(lines)):
        for field in range(len(fields)):
            lineName = lines[i]
            normal = lineNormals[i]
            fieldName = fields[field]

            f = getForce(lineName, normal, fieldName)
            F -= f[0]   #Force
            m -= f[1]   #Mass

    # Berücksichtigung des Massenfehlers:

    F_m = -m * np.array((np.cos(np.radians(angle)), np.sin(np.radians(angle)), 0))
    F_L_m = np.dot(F + F_m, [-np.sin(np.radians(angle)), np.cos(np.radians(angle)), 0])
    F_D_m = np.dot(F + F_m, [np.cos(np.radians(angle)), np.sin(np.radians(angle)), 0])
    c_L_m = F_L_m/(dynamicPressure*Area)
    c_D_m = F_D_m/(dynamicPressure*Area)

    F_L = np.dot(F, [-np.sin(np.radians(angle)), np.cos(np.radians(angle)), 0])
    F_D = np.dot(F, [np.cos(np.radians(angle)), np.sin(np.radians(angle)), 0])

    c_L = F_L/(dynamicPressure*Area)
    c_D = F_D/(dynamicPressure*Area)

    print("F =", F, ", m =", m)
    #print("F_L = ", F_L, ", F_D = ", F_D)
    print("c_L =", c_L, ", c_D =", c_D)
    print("c_L_m =", c_L_m, ", c_D_m =", c_D_m)
    
if __name__ == "__main__":
    main()