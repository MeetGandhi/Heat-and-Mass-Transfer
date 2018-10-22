%explicit considering discretization only in x direction

L = 220E-3;             % total length of domain, m
d = 1.25E-3;            % total thickness of domain, m
w = 22E-3;              % width of plastic film, m
rho = 7850;             % density of material, kg/m3
cp = 435;               % heat capacity, J/kg-K
k = 60;                 % thermal conductivity, W/m-K

T_inf = 25+273.15;             % free stream temperature, degC
T_init = T_inf;         % initial temperature, degC
h = 100;                % convective heat transfer coefficient, W/m2-K

delta_t_ON = 10.0;      % time for which laser is switched on, s
q0_laser = 85E+3;       % laser irradiation, W/m2;
t_final = 30.0;         % total simulation time, s

nx = 51;               % number of points in x-direction
ny = 4;                 % number of points in y-direction

nTimes = 61;            % total number of time steps

nTemps = nx * ny;       % number of temperatures at each time step
dt = t_final/(nTimes-1);% time step, s

dx = L/(nx-1);          % discretization length in x-direction, m
dy = d/(ny-1);          % discretization length in y-direction, m

x = 0:dx:L;             % vector of x locations (for plotting), m
y = 0:dy:d;             % vector of y-locations (for plotting), m
t = 0:dt:t_final;       % vector of times (for plotting), m

TE=zeros(nx,nTimes);
TE(:,1)=T_init;

for i=1:(nTimes-1)
    if dt*i<=delta_t_ON
       TE(1,i+1)=(((k*d)/(dx))*(TE(2,i)-TE(1,i))+h*(T_inf-TE(1,i))*dx/2+(q0_laser*dx/2)+(rho*cp*dx*d*TE(1,i))/2)/(rho*cp*dx*d)/2;
       for j=2:int8((0.022/dx))+1
           TE(j,i+1)=(((k*d)/(dx))*(TE(j+1,i)-TE(j,i))+((k*d)/(dx))*(TE(j-1,i)-TE(j,i))+h*dx*(T_inf-TE(j,i))+(q0_laser*dx)+(rho*cp*dx*d*TE(j,i)))/rho*cp*dx*d;
       end
       for j=int8((0.022/dx))+2:nx-1
           TE(j,i+1)=(((k*d)/(dx))*(TE(j+1,i)-TE(j,i))+((k*d)/(dx))*(TE(j-1,i)-TE(j,i))+h*dx*(T_inf-TE(j,i))+(rho*cp*dx*d*TE(j,i)))/rho*cp*dx*d;
       end
       TE(nx,i+1)=(((k*d)/dx)*(TE(nx-1,i)-TE(nx,i))+h*dx/2*(T_inf-TE(nx,i))+(rho*cp*dx*d*TE(nx,i))/2)/(rho*cp*dx*d)/2;
       
    else
       TE(1,i+1)=(((k*d)/(dx))*(TE(2,i)-TE(1,i))+h*(T_inf-TE(1,i))*dx/2+(rho*cp*dx*d*TE(1,i))/2)/(rho*cp*dx*d)/2;
       for j=2:int8((0.022/dx))+1
           TE(j,i+1)=(((k*d)/(dx))*(TE(j+1,i)-TE(j,i))+((k*d)/(dx))*(TE(j-1,i)-TE(j,i))+h*dx*(T_inf-TE(j,i))+(rho*cp*dx*d*TE(j,i)))/rho*cp*dx*d;
       end
       for j=int8((0.022/dx))+2:nx-1
           TE(j,i+1)=(((k*d)/(dx))*(TE(j+1,i)-TE(j,i))+((k*d)/(dx))*(TE(j-1,i)-TE(j,i))+h*dx*(T_inf-TE(j,i))+(rho*cp*dx*d*TE(j,i)))/rho*cp*dx*d;
       end
       TE(nx,i+1)=(((k*d)/dx)*(TE(nx-1,i)-TE(nx,i))+h*dx/2*(T_inf-TE(nx,i))+(rho*cp*dx*d*TE(nx,i))/2)/(rho*cp*dx*d)/2;
    end
end
