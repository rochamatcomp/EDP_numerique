# 1D ADVECTION EQUATION with FD #
# 1D ADVECTION EQUATION with FD #
# 1D ADVECTION EQUATION with FD #

import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import random


pi=3.14159265358979323846


# PARAMETRES PHYSIQUES
Coeff = 1.  #coefficient physique
Lx = 1.0  #taille du domaine

T = 1. #temps d'integration
print("T final=",T)

#print (np.log(2))

# PARAMETRES NUMERIQUES
NX = 11  #nombre de points de grille
dx = Lx/(NX-1) #pas de grille (espace)
CFL=0.5

dt = CFL*dx    #pas de grille (temps)

NT = int(T/dt)  #nombre de pas de temps
#print("Nombre pas de temps= ",NT)
#print("dt= ",dt)



# Pour la figure
xx = np.zeros(NX)
for i in np.arange(0,NX):
      xx[i]=i*dx

#plt.ion()
#ax = plt.gca(projection='3d')


#Initialisation
ddU   = np.zeros(NX)

# u_0=sin
U_data = np.zeros(NX)
U_old = np.zeros(NX)
U_int = np.zeros(NX)
U_new = np.zeros(NX)

U_data = np.cos(2*pi*xx)
# Creneau
""""
for j in np.arange(0,NX):
    U_data[j]=0.
    if (0.25<xx[j]<0.75): 
        U_data[j]=1.
"""
for i in np.arange(0,NX):
    U_old[i]= U_data[i]
    U_int[i]= U_data[i]

U_new = np.zeros(NX)

#print("init, U_new",U_new)
#print("init, U_old",U_old)
#print (xx)
#print (xx[5])


#plt.figure(-1)
#plt.plot(xx,U_old)
#plt.legend(("t=0."), 'best')
#plt.xlabel('x')
#plt.ylabel('u')
#plt.ylabel('YY')
#plt.show()

# Remplissage de la matrice pour le schema implicite
MX=NX
#MX=5

second_membre=np.ones(MX-1)
res=np.zeros(MX-1)
A=np.zeros((MX-1,MX-1),float)

for j in np.arange(0,MX-1):
    A[j,j]=1.
for j in np.arange(0,MX-2):
    A[j,j+1]=+CFL/2.
    A[j+1,j]=-CFL/2.
A[MX-2,0]=+CFL/2.
A[0,MX-2]=-CFL/2.

    
B=np.linalg.inv(A)
    
#print(A)
#print(B)
#print(second_membre)
#res=np.dot(B,second_membre)
#print(res)


# Boucle en temps

random.seed()


time=0.
for n in np.arange(0,NT):
    
    time=time+dt    
    if (n%10==0):
 #       print("")
 #       print("")
        print ("t=",time)
 #       print ("U_old, begin0",U_old)
        
 

    
# Upwind        
    for j in np.arange(1,NX-1):
        ddU[j] = U_old[j]-U_old[j-1]       
    ddU[0] = U_old[0]-U_old[NX-2]
    
#    print ("U_old, begin1",U_old)
    
    
# Centered
#    for j in np.arange(1,NX-1):
#        ddU[j] = (U_old[j+1]-U_old[j-1])/2.
#    ddU[0] = (U_old[1]-U_old[NX-2])/2.
    
# Lax-Wendroff
    """
    for j in np.arange(1,NX-1):
        ddU[j] = (U_old[j+1]-U_old[j-1]+CFL*(2*U_old[j]-U_old[j+1]-U_old[j-1]))/2.
    ddU[0] = (U_old[1]-U_old[NX-2]+CFL*(2*U_old[0]-U_old[1]-U_old[NX-2]))/2.    
    """
    
#    print ("U_old, begin2",U_old)
    
    
    
   
# Acutalisation   scehams explicites
    for j in np.arange(0,NX-1):
        U_new[j]=U_old[j]-(dt/dx)*ddU[j]
 #   print ("U_old, begin3",U_old)
    
    U_new[NX-1]=U_new[0]
    
#    print ("U_old, begin4",U_old)
    
    U_old=U_new
    

# Implicit
    """
    for j in np.arange(0,NX-1):
        second_membre[j]=U_old[j]        
    res=np.dot(B,second_membre)    
    for j in np.arange(0,NX-1):
        U_new[j]=res[j]        
        
    U_new[NX-1]=U_new[0]
    U_old=U_new

    """

# Monte Carlo
    """"
    random2=random.uniform(0.,1.)
    
    for j in np.arange(1,NX-1):
        random1=random.uniform(0.,1.)
        if (random1<CFL):
            U_int[j]=U_old[j-1]
        else:
            U_int[j]=U_old[j]
    
            
    random1=random.uniform(0.,1.)
    if (random1<CFL):
        U_int[0]=U_old[NX-2]
    else:
        U_int[0]=U_old[0]

        
    U_new[NX-1]=U_new[0]
    
    for j in np.arange(0,NX):
        U_new[j]=U_int[j]
        U_old[j]=U_new[j]
    """
    
# Glimm:  a faire


    

    
    

    
            
        
            
            

 
   
    #toto=U_old
#    U_data=U_new
   # U_new=toto

# FIXED BC
#   U_data[0,:]=0.0



# PERIODIC BC
#   U_data[0,:]=U_data[NX-2,:]
#   U_data[NX-1,:]=U_data[1,:]
#   U_data[:,0]=U_data[:,NY-2]
#   U_data[:,NY-1]=U_data[:,1]

   # if (n%10==0):
       
 #       print ("A3",n)
     
#plt.figure(n)

# caclul erreur numerique
norme_L2=0.
for j in np.arange(0,NX-1):
    norme_L2=norme_L2+ dx*np.fabs(U_new[j]-U_int[j])**2
norme_L2=np.sqrt(norme_L2)
print ("Points", NX-1, "Erreur/norme L2",norme_L2)



print ("tFinal=",time)
#print(U_new)   
plt.figure(0)

plt.plot(xx,U_new,"r",marker='x')
plt.plot(xx,U_data,"b",marker='o')
plt.legend(("t=T"), 'best')
plt.xlabel('x')
plt.ylabel('u')
      #  plt.draw()
    # plt.show()
     
 #       print ("A34",n)
     
       
  #   plotlabel= "N = " + str(n+1)

     #ax.cla()
     #ax.plot_surface(xx,yy,U_data,vmin=-0.1,vmax=0.1,cmap=cm.jet,antialiased=False,linewidth=0,rstride=1,cstride=1)
     #ax.set_zlim3d(-0.1,0.1)

 #    plt.pcolormesh(xx,yy,U_data)
#     plt.axis('image')
 #    plt.clim(-0.1,0.1)

#     plt.title(plotlabel)
 #    plt.draw()
 #    if 'qt' in plt.get_backend().lower():
 #       QtGui.qApp.processEvents()

plt.show()

