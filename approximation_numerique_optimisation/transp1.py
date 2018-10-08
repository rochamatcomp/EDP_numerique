#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Grégoire Allaire, Flore Nabet et Thomas Wick
# Ecole Polytechnique
# MAP 411
# Hiver 2017/2018

# Résolution numérique de l'équation de transport
# u,t + v u,x = 0
# schéma explicite centré
########################################################
# Load packages 
import numpy as np
import matplotlib.pyplot as plt 

import pylab
import scipy as sp
from scipy import sparse
from scipy.sparse import linalg


########################################################
# Paramètres du probleme
v = 1.0        # v = vitesse
lg = 10.        # intervalle en x=[-lg,lg]
nx = 201       # nombre de points du maillage
dx = (2*lg)/(nx-1)       # dx = pas d'espace
cfl = 0.9     # cfl = dt/dx
dt = dx*cfl/v # dt = pas de temps
Tfinal = 2.   # Temps final souhaité

print("parametres de la discretisation : domaine (", -lg , lg , ") discretise avec nx=", nx, " points de discretisation et une taille de maille dx=", dx)

x = np.linspace(-lg,lg,nx)

# Initialize u0
u0 = np.zeros(len(x))
#print(len(u0))

# Set specific u0 values (same as in the scilab program)
for i in range (len(x)):
    if (1.0 - x[i]**2) < 0:
       u0[i] = 0
    else:
       u0[i] = 1.0 - x[i]**2 # donnée initiale


## Plot initial condition and save it to a file
plt.figure(1)
plt.clf()
plt.plot(x,u0) # trace u0 en fonction de x
plt.title('condition initiale, transport')
plt.xlabel('x')
plt.ylabel('u(x)')
#plt.savefig('u0.png')
plt.show()
plt.pause(1.)

########################################################
# Schéma numérique

# Initialize u by the initial data u0
u = u0.copy() # il faut faire une copie sinon on va changer u0 en changeant u

# Construction de la matrice (creuse) pour le schema explicite
# u(-10,t)=u(10,t)=0 on peut donc enlever ces deux coordonnees du systeme
#u_j^{n+1}=nu*cfl*u_{j-1}^n+(1-2*nu*cfl)u_j^n+nu*cfl*u_{j+1}^n
A=sp.sparse.diags([v*dt/(2*dx), 1., -v*dt/(2*dx)], [-1, 0, 1], shape=(nx-2, nx-2))

# Nombre de pas de temps effectues
nt = int(Tfinal/dt)
Tfinal = nt*dt # on corrige le temps final (si Tfinal/dt n'est pas entier)

# Time loop
for n in range(1,nt+1):

  # Schéma explicite centré
    u[1:len(u)-1]=A*u[1:len(u)-1]

 # Print solution
    if n%2 == 0:
      plt.figure(1)
      plt.clf()
      plt.plot(x,u0,'b',x,u,'r')
      plt.xlabel('$x$')
      plt.title('Schema explicite centre, $t=$%s' %(n*dt))
      plt.pause(0.1)

####################################################################
# comparaison solution exacte avec solution numerique au temps final

uexacte = np.zeros(len(u))

# on calcule la solution exacte qui est la donnée initiale transportée
for i in range(int(nx)):
    j = i - int(Tfinal*v/dx) ;
    j = max(j,1)
    uexacte[i] = u0[j]

plt.figure(2)
plt.clf()
plt.suptitle('Transport: comparaison au temps Tfinal$=$%s' %(Tfinal))
#plt.suptitle('Comparaison entre la solution exacte et la solution numerique au temps final $T=$ %s', %(Tfinal))
plt.plot(x,u0,'b',x,u,'or',x,uexacte,'k')
plt.legend(['Donnee initiale','Schema explicite centre','Solution exacte'],loc='best')
plt.show()
