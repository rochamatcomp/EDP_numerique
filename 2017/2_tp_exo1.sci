clear
// champ de vitesse
deff('y=a(x)','y=x')

// position initiale des trajectoires 
dx=0.2
x0=[-2:dx:2]
J=length(x0)

// echelle de temps
Tend = 1
N=30
dt=Tend/N
T=[0:dt:Tend]

// stockage des trajectoires en espace et temps
xT=zeros(J,N+1)

// position init des trajectoires
xT(:,1)=x0(:)

// ouverture fenetre graphique
xset("window",0);
xtitle('trajectoires','x','t')

// boucle en temps
for n=1:N
    // boucle sur les positions
    for j=1:J
        xT(j,n+1)=xT(j,n)+a(xT(j,n))*dt
    end
    // plot des nouvelles positions 
    Tplot = T(n).*ones(1,J)
    plot2d(xT(:,n),Tplot,-1)
    sleep(300)
end

// plot des caracterisitiques
for j=1:J
    plot2d(xT(j,:),T,2)
end





