"""
TP1 Numerical analysis
"""   
# 1D ADVECTION EQUATION with FD #
# 1D ADVECTION EQUATION with FD #
# 1D ADVECTION EQUATION with FD #

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import random


def initial_data():
    """
    Initialization of data.
    """
    

def lax_wendroff_scheme(U_new, U_old, ddU, NX, CFL):
    """
    Lax-Wendroff_scheme.
    """
    for j in np.arange(1,NX-1):
        ddU[j] = (U_old[j+1] - U_old[j-1] + CFL*(2*U_old[j] - U_old[j+1] - U_old[j-1]))/2.
    
    ddU[0] = (U_old[1] - U_old[NX-2] + CFL*(2*U_old[0] - U_old[1] - U_old[NX-2]))/2.
       
    update_scheme(U_new, U_old, ddU, NX, CFL)
    
def centered_explicit_scheme(U_new, U_old, ddU, NX, CFL):
    """
    Centered explicit scheme.
    """
    for j in np.arange(1, NX-1):
        ddU[j] = (U_old[j+1] - U_old[j-1])/2.0

    ddU[0] = (U_old[1] - U_old[NX-2])/2.0
       
    update_scheme(U_new, U_old, ddU, NX, CFL)
    
def centered_implicit_scheme(U_new, U_old, ddU, NX, CFL):
    """
    Centered implicit scheme.
    """
    for j in np.arange(0,NX-1):
        second_membre[j] = U_old[j]        
        
    res = np.dot(B, second_membre)    
    for j in np.arange(0,NX-1):
        U_new[j] = res[j]
        
    update_scheme(U_new, U_old, ddU, NX, CFL)

def upwind_scheme(U_new, U_old, ddU, NX, CFL):
    """
    Upwind scheme
    """
    for j in np.arange(1, NX-1):
        ddU[j] = (U_old[j] - U_old[j-1])

    ddU[0] = (U_old[1] - U_old[NX-2])
    
    update_scheme(U_new, U_old, ddU, NX, CFL)
    
def error():
    """
    Error Estimation
    """
    pass

def montecarlo_scheme(U_new, U_int, U_old, CFL, NX):
    """
    Monte Carlo scheme.
    """        
    for j in np.arange(1, NX):
        random1 = random.uniform(0.,1.)
        
        if random1 < CFL:
            U_int[j] = U_old[j-1]
        else:
            U_int[j]=U_old[j]
    
    random1 = random.uniform(0.,1.)
    
    if random1 < CFL:
        U_int[0] = U_old[NX-2]
    else:
        U_int[0] = U_old[0]
        
    update_probabilist_scheme(U_new, U_int, U_old, NX)
    
def glimm_scheme(U_new, U_int, U_old, CFL, NX):
    """
    Glimm scheme.
    """
    random1 = random.uniform(0.,1.)
    
    for j in np.arange(1, NX):    
        if random1 < CFL:
            U_int[j] = U_old[j-1]
        else:
            U_int[j]=U_old[j]    
    
    if random1 < CFL:
        U_int[0] = U_old[NX-2]
    else:
        U_int[0] = U_old[0]
        
    update_probabilist_scheme(U_new, U_int, U_old, NX)

def update_scheme(U_new, U_old, ddU, NX, CFL):
    """
    
    """
    for j in np.arange(0,NX-1):
        U_new[j] = U_old[j] - CFL * ddU[j]
     
    U_new[NX-1] = U_new[0]
    U_old = U_new
    


def update_probabilist_scheme(U_new, U_int, U_old, NX):
    """
    Update
    """    
    U_new[NX-1] = U_new[0]
    
    for j in np.arange(0,NX):
        U_new[j]=U_int[j]
        U_old[j]=U_new[j]

def execute_scheme():
    """
    Execute the scheme
    """
    time = 0.
    for n in np.arange(0, NT):
        
        time = time + dt    
        if (n % 100==0):
            print ("n=",n,": t=",time)
        
        #lax_wendroff_scheme(U_new, U_old, ddU, NX, CFL)
        #centered_explicit_scheme(U_new, U_old, ddU, NX, CFL)
        #centered_implicit_scheme(U_new, U_old, ddU, NX, CFL)
        upwind_scheme(U_new, U_old, ddU, NX, CFL)
                        
        #montecarlo_scheme(U_new, U_int, U_old, CFL, NX)
        #glimm_scheme(U_new, U_int, U_old, CFL, NX)    
    
        #Error        
        
    print ("tFinal=",time)

def main():
    """
    Definitions of constants and schemes
    """
    #pi = 3.14159265358979323846
    pi = np.pi
    
    # PARAMETRES PHYSIQUES
    Coeff = 1.  #coefficient physique
    Lx = 1.0  #taille du domaine
    
    T = 1.3 #temps d'integration
    print("T final=",T)
    
    
    # PARAMETRES NUMERIQUES
    NX = 101  #nombre de points de grille
    dx = Lx/(NX-1) #pas de grille (espace)
    
    # Courant index
    CFL = 0.2
    
    dt = CFL * dx    #pas de grille (temps)
    
    NT = int(T/dt)  #nombre de pas de temps
    
    print("Espace steps: ", NX)
    print("Time steps: ", NT)
    print("dx= ",dx)
    print("dt= ",dt)
    
    
    # Pour la figure
    xx = np.zeros(NX)
    for i in np.arange(0, NX):
          xx[i] = i*dx
    
    #Initialisation
    ddU   = np.zeros(NX)
    
    # Initial data u_0 = sin
    U_data = np.zeros(NX)
    U_old = np.zeros(NX)
    U_int = np.zeros(NX)
    U_new = np.zeros(NX)
    
    U_data = np.sin(2*pi*xx)
    
    
    # Creneau
    """"
    for j in np.arange(0,NX):
        U_data[j]=0.
        if (0.25<xx[j]<0.75): 
            U_data[j]=1.
    """
    
    for i in np.arange(0, NX):
        U_old[i]= U_data[i]
        U_int[i]= U_data[i]
    
    U_new = np.zeros(NX)
       
    # Remplissage de la matrice pour le schema implicite
    MX = NX
    
    second_membre = np.ones(MX-1)
    res = np.zeros(MX-1)
    A = np.zeros((MX-1,MX-1),float)
    
    for j in np.arange(0,MX-1):
        A[j,j]=1.
    for j in np.arange(0,MX-2):
        A[j,j+1]=+CFL/2.0
        A[j+1,j]=-CFL/2.0
         
    A[MX-2,0] = +CFL/2.0
    A[0,MX-2] = -CFL/2.0
    
    B = np.linalg.inv(A)
         
    # Boucle en temps
    random.seed()
    
    execute_scheme()
    
    plt.figure(-1)
    plt.plot(xx,U_old)
    plt.legend(("t=0."), 'best')
    plt.xlabel('x')
    plt.ylabel('u')

    plt.figure(0)
    
    plt.plot(xx,U_new,"r",marker='x')
    plt.legend(("t=T"), 'best')
    plt.xlabel('x')
    plt.ylabel('u')
    
    plt.show()

if __name__ == '__main__':
    main()