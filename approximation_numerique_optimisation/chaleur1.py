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
# schéma centré instable
########################################################
# Load packages (need to be installed first if
# not yet done - but is not difficult)
import numpy as np
import matplotlib.pyplot as plt # pour plot functons
#plt.switch_backend('tkAgg')  # necessary for OS SUSE 13.1 version,
# otherwise, the plt.show() function will not display any window

import pylab
import scipy as sp
#from scipy import sparse
from scipy.sparse import linalg


########################################################
# Paramètres du probleme
lg = 10.        # intervalle en x=[-lg,lg]
nx = 201       # nombre de points du maillage
dx = (2*lg)/(nx-1)       # dx = pas d'espace
cfl = 0.1     # cfl = dt/dx^2
dt = dx*dx*cfl # dt = pas de temps
Tfinal = 0.03   # Temps final souhaité

print("schema centre en temps")
print("parametres : domaine (", -lg , lg , ") discretise avec nx=", nx, " points et une taille de maille dx=", dx)

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


# Set specific values at time step -1
u1 = u0.copy() # il faut faire une copie sinon on va changer u0 en changeant u1
uinit = u0.copy()

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

# Construction de la matrice (creuse) pour le schema explicite
# u(-10,t)=u(10,t)=0 on peut donc enlever ces deux coordonnees du systeme
#u_j^{n+1}=nu*cfl*u_{j-1}^n+(1-2*nu*cfl)u_j^n+nu*cfl*u_{j+1}^n
A=sp.sparse.diags([2.*cfl, -4.*cfl, 2.*cfl], [-1, 0, 1], shape=(nx-2, nx-2))


# Nombre de pas de temps effectues
nt = int(Tfinal/dt)
Tfinal = nt*dt # on corrige le temps final (si Tfinal/dt n'est pas entier)

# Time loop
for n in range(1,nt+1):

  # Schéma centré en temps
  u[1:len(u)-1]=u1[1:len(u)-1] + A*u0[1:len(u)-1]
  u1 = u0.copy()
  u0 = u.copy()

 # Print solution
  if n%2 == 0:
      plt.figure(1)
      plt.clf()
      plt.plot(x,uinit,'b',x,u,'r')
      plt.xlabel('$x$')
      plt.title('Schema centre, $t=$%s' %(n*dt))
      plt.pause(1.)

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
        uexacte[i] = uexacte[i] + uinit[j]*dx*noyauc((i-j)*dx,Tfinal)


plt.pause(1.)
plt.figure(2)
plt.suptitle('Comparaison entre solutions exacte et approchee au temps Tfinal$=$%s' %(Tfinal))
#plt.suptitle('Comparaison entre la solution exacte et la solution numerique au temps final $T=$ %s', %(Tfinal))
plt.plot(x,u0,'b',x,u,'or',x,uexacte,'k')
plt.legend(['Donnee initiale','Schema centre','Solution exacte'],loc='best')
plt.show()
