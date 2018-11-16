#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Homework B004
Numerical methods for non-stationary Partial Differential Equations (PDEs):
Finite Differences and Finite Volumes

Mathematics and Modeling
Master of Mathematics and Applications

@author: Andre ALMEIDA ROCHA (3701739) <rocha.matcomp@gmail.com>
"""

import numpy as np
import matplotlib.pyplot as plt

def main():
    points = [10, 20, 40, 80, 160, 320]

    errors = []
    for NX in points:
        norme = difference_finites(NX)
        errors.append(norme)

    print(errors)

    order1 = [1/point for point in points]
    order2 = [1/(point**2) for point in points]
    order4 = [1/(point**4) for point in points]

    plt.figure(0)
    line1 = plt.plot(points, order1, label = 'First order convergence')
    line2 = plt.plot(points, order2, label = 'Second order convergence')
    line3 = plt.plot(points, order4, label = 'Fourth order convergence')
    line4 = plt.plot(points, errors, "r", label = 'Errors')
    plt.legend(handles = [line1[0], line2[0], line3[0], line4[0]])

    plt.figure(1)
    linelog1 = plt.loglog(points, order1, label = 'Log of first order convergence')
    linelog2 = plt.loglog(points, order2, label = 'Log of second order convergence')
    linelog3 = plt.loglog(points, order4, label = 'Log of fourth order convergence')
    linelog4 = plt.loglog(points, errors, "r", label = 'Log of Errors')
    plt.legend(handles = [linelog1[0], linelog2[0], linelog3[0], linelog4[0]])

    plt.show()

def difference_finites(NX):
    """
    Finite differences method
    NX : nombre de points de grille
    """
    # PARAMETRES PHYSIQUES
    Lx = 1.0  #taille du domaine
    #T = 0.4 #temps d'integration
    T = 0.0085 #temps d'integration

    print("T final=",T)

    # PARAMETRES NUMERIQUES
    dx = Lx/NX #pas de grille (espace)
    CFL = 0.42
    dt = CFL * dx**2    #pas de grille (temps)
    NT = int(T/dt)  #nombre de pas de temps

    print("Nombre pas de espace= ",NX)
    print("dx= ", dx)
    print("Nombre pas de temps= ",NT)
    print("dt= ", dt)

    # PARAMETRES FONCTION INITIAL
    alpha = 1
    beta = 1

    # Pour la figure
    xx = np.zeros(NX+1)
    for i in np.arange(0,NX+1):
          xx[i] = i*dx

    #Initialisation
    ddU   = np.zeros(NX+1)
    U_old =  np.zeros(NX+1)
    U_new =  np.zeros(NX+1)

    # Condition initial
    U_data = np.zeros(NX+1)
    U_data = alpha * np.sin(2*np.pi*xx) + beta * np.cos(2*np.pi*xx)

    # Solution exacte
    U_sol =  np.zeros(NX+1)
    U_sol = np.exp(-4*(np.pi**2) * T) * U_data

    for i in np.arange(0, NX+1):
        U_old[i] = U_data[i]

    plt.figure(-1)
    line = plt.plot(xx, U_data, label = 'best')
    plt.legend(handles = [line[0]])
    plt.xlabel('x')
    plt.ylabel('u')

    # Boucle en temps
    time=0.
    for n in np.arange(0, NT):

        time=time+dt
        if (n%100 == 0):
            print ("t=",time)

        # schemas explicites
        for j in np.arange(2, NX-1):
            ddU[j] = 4/3*(U_old[j+1] - 2*U_old[j] + U_old[j-1]) - 1/12*(U_old[j+2] - 2*U_old[j] + U_old[j-2]) + 0.5*CFL*(U_old[j+2] - 4*U_old[j+1] + 6*U_old[j] - 4*U_old[j-1] + U_old[j-2])

        ddU[NX-1] = 4/3*(U_old[0] - 2*U_old[NX-1] + U_old[NX-2]) - 1/12*(U_old[1] - 2*U_old[NX-1] + U_old[NX-3]) + 0.5*CFL*(U_old[1] - 4*U_old[0] + 6*U_old[NX-1] - 4*U_old[NX-2] + U_old[NX-3])
        ddU[1] = 4/3 * (U_old[2] - 2 * U_old[1] + U_old[0]) - 1/12 * (U_old[3] - 2 * U_old[1] + U_old[NX-1]) + 0.5*CFL*(U_old[3] - 4*U_old[2] + 6*U_old[1] - 4*U_old[0] + U_old[NX-1])
        ddU[0] = 4/3 * (U_old[1] - 2 * U_old[0] + U_old[NX-1]) - 1/12 * (U_old[2] - 2 * U_old[0] + U_old[NX-2]) + 0.5*CFL*(U_old[2] - 4 * U_old[1] + 6 * U_old[0] - 4 * U_old[NX-1] + U_old[NX-2])


        # Actualisation schemas explicites
        for j in np.arange(0, NX):
            U_new[j] = U_old[j] - CFL * ddU[j]

        # Actualisation pour la 1-periodicite
        U_new[NX] = U_new[0]

        # Actualisation en temps
        for j in np.arange(0, NX):
            U_new[j] = U_old[j]


    print ("tFinal=", time)


    # caclul erreur numerique
    norme = np.sqrt(dx * np.sum(np.fabs(U_new - U_sol)**2))

    #for j in np.arange(0, NX-1):
    #    norme = norme + dx*np.fabs(U_new[j] - U_sol[j])**2
    #norme = np.sqrt(norme)

    print ("Erreur/norme L2",norme)
    plt.figure(0)

    line1 = plt.plot(xx, U_new, "r", marker='x', label = 'discret')
    line2 = plt.plot(xx, U_sol, "b", marker='o', label = 'exact')
    plt.legend(handles = [line1[0], line2[0]])
    plt.xlabel('x')
    plt.ylabel('u')

    plt.show()

    return norme

if __name__ == '__main__':
    main()
