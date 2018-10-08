#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Grégoire Allaire, Flore Nabet et Thomas Wick
# Ecole Polytechnique
# MAP 411
# Hiver 2017/2018

# Execution of *.py files
# Possiblity 1: in terminal
#  terminal> python3 file.py  # here file.py = convection.py

# Possiblity 2: executing in an python environment
# such as spyder.

# Résolution numérique de l'équation de la chaleur
# u,t - u,xx = 0
# schéma implicite
########################################################
# Load packages (need to be installed first if
# not yet done - but is not difficult)
import numpy as np
import matplotlib.pyplot as plt # pour plot functons
#plt.switch_backend('tkAgg')  # necessary for OS SUSE 13.1 version,
# otherwise, the plt.show() function will not display any window

import pylab
import scipy as sp
from scipy import sparse
from scipy.sparse import linalg


########################################################
# Paramètres du probleme
lg = 10.        # intervalle en x=[-lg,lg]
nx = 201       # nombre de points du maillage
dx = (2*lg)/(nx-1)       # dx = pas d'espace
cfl = 2.0     # cfl = dt/dx^2
dt = dx*dx*cfl # dt = pas de temps
Tfinal = 1.   # Temps final souhaité

print("parametres de la discretisation : domaine (", -lg , lg , ") discretise avec nx=", nx, " points de discretisation et une taille de maille dx=", dx)

x = np.linspace(-lg,lg,nx)

# Initialize u0
u0 = np.zeros(len(x))
#print(len(u0))

# Set specific u0 values (same as in the scilab program)
for k in range (len(x)):
    if (1.0 - x[k]**2) < 0:
       u0[k] = 0
    else:
       u0[k] = 1.0 - x[k]**2 # donnée initiale


## Plot initial condition and save it to a file
plt.figure(1)
plt.clf()
plt.plot(x,u0) # trace u0 en fonction de x
plt.title('condition initiale')
plt.xlabel('x')
plt.ylabel('u(x)')
#plt.savefig('u0.png')
plt.show()
plt.pause(1.)

########################################################
# Schemas numeriques

# Initialize u by the initial data u0
u = u0.copy() # il faut faire une copie sinon on va changer u0 en changeant u

# Construction de la matrice (creuse) pour le schema implicite
#-nu*cfl*u_{j-1}^{n+1}+(1+2*nu*cfl)u_j^{n+1}-nu*cfl*u_{j+1}^{n+1}=u_j^n
Aimp=sp.sparse.diags([-cfl, 1+2*cfl, -cfl], [-1, 0, 1], shape=(nx, nx),format = 'csc')
lu_Aimp = linalg.splu(Aimp, permc_spec = 'NATURAL')

#lu_Aimp est un objet de type scipy.sparse.linalg.SuperLU qui contient:
# - une matrice triangulaire inferieure lu_Aimp.L (lower)
# - une matrice triangulaire superieure lu_Aimp.U (upper)
# - et eventuellement des matrices de permutation mais on a bloque cette possibilite ici avec l'option permc_spec = 'NATURAL'
#   Le resultat est la factorisation Aimp = lu_Aimp.L * lu_Aimp.U
# l'interet est qu'il est peu couteux de resoudre des systemes lineaires pour des matrices triangulaires



# Nombre de pas de temps effectues
nt = int(Tfinal/dt)
Tfinal = nt*dt # on corrige le temps final (si Tfinal/dt n'est pas entier)

# Time loop
for n in range(1,nt+1):

  # Schéma implicite en temps
    u = lu_Aimp.solve(u)
  # resout les deux systemes lineaires avec matrices triangulaires necessaires pour calculer uimp = Aimp^-1 uimp
  #  uimp = sp.sparse.linalg.spsolve(Aimp,uimp)

 # Print solution
    if n%1 == 0:
      plt.figure(1)
      plt.clf()
      plt.plot(x,u0,'b',x,u,'r')
      plt.xlabel('$x$')
      plt.title('Schema implicite, $t=$%s' %(n*dt))
      plt.pause(0.1)

        # Print solution into seperate files
#        fig, ax = plt.subplots(nrows=1,ncols=1)
#        ax.plot(x,uexp,'b',x,uimp,'r') # trace uexp et uimp en fonction de x
#        # On peux utiliser *.png ou *.pdf
#        fig.savefig(str(n) + 'u.png')
#        plt.close(fig)

####################################################################
# comparaison solution exacte avec solution numerique au temps final

uexacte = np.zeros(len(uinit))

def noyauc(x,t):
    return np.exp(-x**2/(4*t))/np.sqrt(4*np.pi*t)

# on calcule la solution exacte en utilisant la formule des rectangles
# uexacte(xi,Tfinal)=\int u0(y) noyauc(xi-y,Tfinal)
# uexacte(xi,Tfinal)~sum_{j=0}^{2*lg/dx} dx*u0(xj) noyauc(xi-xj,Tfinal)
# et xi=-lg+i*dx donc xi-xj=(i-j)*dx
for i in range(int(nx)):
    for j in range(int(nx)-1):
        uexacte[i] = uexacte[i] + u0[j]*dx*noyauc((i-j)*dx,Tfinal)


plt.figure(2)
plt.clf()
plt.suptitle('Comparaison entre solutions exacte et approchee au temps Tfinal$=$%s' %(Tfinal))
#plt.suptitle('Comparaison entre la solution exacte et la solution numerique au temps final $T=$ %s', %(Tfinal))
plt.plot(x,u0,'b',x,u,'or',x,uexacte,'k')
plt.legend(['Donnee initiale','Schema implicite','Solution exacte'],loc='best')
plt.show()
