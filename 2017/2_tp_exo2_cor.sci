clear
// schema VF decentre ammont pour eq transport a coeff variable 1D

// *** vitesse ***
deff('y=a(x)','y=((x))')

// *** space mesh size ***
J=200 
dx = 4/J
X=[-2:dx:2]

// *** vitesse aux nodes ***
aj=a(X)
A = max(abs(aj))

// *** time mesh size - CFL ***
Tend=0.8
dt=dx/(2*A)
T=0:dt:Tend
N=length(T)

// *** centres VF ***
Xvf=X+dx/2
Xvf=Xvf(1:$-1)

// *** fonction initiale ***
deff('y=u0(x)','if x < -0.5 then, y=0, elseif x>0.5 then, y=0, else, y=1, end,')


deff('y=ufin(x)','if x < -0.5*exp(.8) then, y=0, elseif x>0.5*exp(.8) then, y=0, else, y=exp(-.8), end,')
// *** solution discrete ***
U=zeros(J,N)
Ufin=zeros(J)

// initialisation
for j=1:J, U(j,1)=u0(Xvf(j)); Ufin(j)=ufin(Xvf(j)); end;

// scheme
for n=1:N-1
    //calcul du flux
    for j=1:J-1
        flux(j)=max(aj(j+1),0)*U(j,n)-max(-aj(j+1),0)*U(j+1,n);
        fluxd(j)=max(-aj(j+1),0)*(U(j,n)-U(j+1,n));   
        fluxg(j)=max(aj(j+1),0)*(U(j,n)-U(j+1,n));   
        end
    U(:,n+1)=dx*U(:,n);
    for j=1:J-1
        U(j,n+1)=U(j,n+1)-fluxd(j)*dt;

        U(j+1,n+1)=U(j+1,n+1)+fluxg(j)*dt;
    end
    U(:,n+1)=U(:,n+1)/dx;
end

// plot a differents temps
t1=floor(N/4) ; t2=floor(N/2) ; t3=floor(3*N/4)

Xvf=Xvf'

xset("window",0) ;
plot2d([Xvf,Xvf],[U(:,1),U(:,1)],[1,-1])
xtitle('time init','x','u')

xset("window",1) ;
plot2d([Xvf,Xvf],[U(:,t1),U(:,t1)],[1,-1])
xtitle('time 1','x','u')

xset("window",2) ;
plot2d([Xvf,Xvf],[U(:,t2),U(:,t2)],[1,-1])
xtitle('time 2','x','u')

xset("window",3) ;
plot2d([Xvf,Xvf],[U(:,t3),U(:,t3)],[1,-1])
xtitle('time 3','x','u')

xset("window",4) ; 
plot2d([Xvf,Xvf],[U(:,$),U(:,$)],[1,-1])
xtitle('time final','x','u')






