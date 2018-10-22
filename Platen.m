clear all;
clc;
n=2.0;
L=0.02;
W=n*L;
LD=6.0*(W+L);
d=5.0*(L/2.0);
w=L;
rho=8000.0;
cp=500.0;
k=15.0;

T_inf_o=25.0;
T_init=25.0;
T_inf_i=200.0;
h_o=250.0;
h_i=500.0;

t_final = 30.0;
nTimes=61.0;
dt=t_final/(nTimes-1);

dx=0.004;
dy=0.005;

nx=(30.0*(n+1))+1;
ny=11;
nTemps=nx*ny;
x = 0:dx:L;             
y = 0:dy:d;             
t = 0:dt:t_final;       

T = zeros(nTemps, nTimes);
A = zeros(nTemps, nTemps);
B = A;
C = zeros(nTemps, 1);
T(:, 1) = T_init;

start=((W/(2.0*dx))+1)+((L/dy)*nx);

v1=[];
for e=start+nx:nx:start+((L/dy)*(nx))-nx
    v1=[v1,e];
end

v2=[];
for e=start+(L/dx)+nx:nx:start+(L/dx)+((L/dy)*(nx))-nx
    v2=[v2,e];
end

for i=1:nTemps
    
    if i==1
        A(i,i+1)=(k*dy)/(2.0*dx);
        A(i,i+nx)=(k*dx)/(2.0*dy);
        B(i,i)=-1*(rho*cp*(dy/2.0)*(dx/2.0))/dt;
        C(i)=-1.0*h_o*(dx/2)*T_inf_o;
        A(i,i)= -(k*dy)/(2.0*dx) - (k*dx)/(2.0*dy) -(rho*cp*(dy/2.0)*(dx/2.0))/dt - h_o*(dx/2);
    
    elseif i<nx
        A(i,i-1)=(k*dy/(2.0*dx));
        A(i,i+1)=(k*dy/(2.0*dx));
        A(i,i+nx)=(k * dx/ dy);
        B(i,i)=(-1.0* rho * cp * dx * (dy/2.0)/ dt);
        C(i)=(- h_o * dx * T_inf_o);
        A(i,i)= -1.0*( (k*dy/(2.0*dx)) + (k*dy/(2.0*dx)) + (k * dx/ dy) + (rho * cp * dx * (dy/2.0)/ dt) +  h_o * dx  );
        
    elseif i==nx
        A(i,i-1)=(k*dy)/(2.0*dx);
        A(i,i+nx)=(k*dx)/(2.0*dy);
        B(i,i)=(-1.0* rho * cp * (dx/2.0) * (dy/2.0)/ dt);
        C(i) = - h_o * T_inf_o * ((dx/2.0) + (dy/2.0));
        A(i,i) = -1.0 * ( (k*dy)/(2.0*dx) + (k*dx)/(2.0*dy) +  (rho * cp * (dx/2.0) * (dy/2.0)/ dt)+ h_o *  ((dx/2.0) + (dy/2.0))) ;
        
    elseif rem(i-1,nx)==0 && i~=1 && i~=nTemps-nx+1
        A(i,i-nx)=(k*dx)/(2.0*dy);
        A(i,i+1)=(k*dy)/dx;
        A(i,i+nx)=(k*dx)/(2.0*dy);
        B(i,i)=(- rho * cp * (dx/2.0) * (dy)) / dt ;
        A(i,i)=-(k*dx)/(2.0*dy)-(k*dx)/(2.0*dy)-(k*dy)/dx+(- rho * cp * (dx/2.0) * (dy)) / dt ;       
    
    elseif i==nTemps-nx+1
        A(i,i-nx)=(k*dx)/(2.0*dy);
        A(i,i+1)=(k*dy)/(2.0*dx);
        B(i,i)=(- rho * cp * (dx/2.0) * (dy/2.0)) / dt ;
        A(i,i)=-(k*dy)/(2.0*dx)-(k*dx)/(2.0*dy)+(- rho * cp * (dx/2.0) * (dy/2.0)) / dt ;
        
    elseif rem(i,nx)==0 && i~=91 && i~=nTemps
        A(i,i-nx)=(k*dx)/(2.0*dy);
        A(i,i+nx)=(k*dx)/(2.0*dy);
        A(i,i-1)=(k*dy)/(dx);
        B(i,i)=(- rho * cp * (dx/2.0) * (dy)) / dt;
        C(i)=-h_o*(dy)*T_inf_o;
        A(i,i)=(- rho * cp * (dx/2.0) * (dy)) / dt - (k*dx)/(2.0*dy) - (k*dx)/(2.0*dy) - (k*dy)/(dx) -h_o*(dy);
        
    elseif i==nTemps
        A(i,i-1)=(k*dy)/(2.0*dx);
        A(i,i-nx)=(k*dx)/(2.0*dy);
        B(i,i)=(- rho * cp * (dx/2.0) * (dy/2.0)) / dt;
        C(i)=-h_o*(dy/2.0)*T_inf_o;
        A(i,i)=-(k*dy)/(2.0*dx)-(k*dx)/(2.0*dy)+(- rho * cp * (dx/2.0) * (dy/2.0)) / dt -h_o*(dy/2.0) ;
        
    elseif i<nTemps && i>nTemps-nx+1
        A(i,i-1)=(k*dy)/(2.0*dx);
        A(i,i+1)=(k*dy)/(2.0*dx);
        A(i,i-nx)=(k*dx)/dy;
        B(i,i)=(- rho * cp * (dx) * (dy/2.0)) / dt;
        A(i,i)=-(k*dy)/(2.0*dx)-(k*dy)/(2.0*dx)-(k*dx)/dy + (- rho * cp * (dx) * (dy/2.0)) / dt;
        
    elseif i==start
        A(i,i-1)=(k*dy)/dx;
        A(i,i-nx)=(k*dx)/dy;
        A(i,i+1)=(k*dy)/(2.0*dx);
        A(i,i+nx)=(k*dx)/(2.0*dy);
        C(i)=-h_i*((dx/2.0)+(dy/2.0))*T_inf_i;
        B(i,i)=(- rho * cp * (0.7500 * dx * dy)) / dt;
        A(i,i)=-((k*dy)/dx)-((k*dx)/dy)-((k*dx)/(2.0*dy))-((k*dy)/(2.0*dx))-h_i*((dx/2.0)+(dy/2.0))+(- rho * cp * (0.7500 * dx * dy)) / dt;
        
    elseif i==start+(L/dx)
        A(i,i-1)=(k*dy)/(2.0*dx);
        A(i,i+1)=(k*dy)/dx;
        A(i,i-nx)=(k*dx)/dy;
        A(i,i+nx)=(k*dx)/(2.0*dy);
        C(i)=-h_i*((dx/2.0)+(dy/2.0))*T_inf_i;
        B(i,i)=(- rho * cp * (0.7500 * dx * dy)) / dt;
        A(i,i)=-((k*dy)/dx)-((k*dx)/dy)-((k*dx)/(2.0*dy))-((k*dy)/(2.0*dx))-h_i*((dx/2.0)+(dy/2.0))+(- rho * cp * (0.7500 * dx * dy)) / dt;
        
    elseif i==start+((L/dy)*(nx))
        A(i,i-1)=(k*dy)/dx;
        A(i,i+nx)=(k*dx)/dy;
        A(i,i-nx)=(k*dx)/(2.0*dy);
        A(i,i+1)=(k*dy)/(2.0*dx);
        C(i)=-h_i*((dx/2.0)+(dy/2.0))*T_inf_i;
        B(i,i)=(- rho * cp * (0.7500 * dx * dy)) / dt;
        A(i,i)=-((k*dy)/dx)-((k*dx)/dy)-((k*dx)/(2.0*dy))-((k*dy)/(2.0*dx))-h_i*((dx/2.0)+(dy/2.0))+(- rho * cp * (0.7500 * dx * dy)) / dt;
        
    elseif i==start+(L/dx)+((L/dy)*(nx))
        A(i,i+1)=(k*dy)/dx;
        A(i,i+nx)=(k*dx)/dy;
        A(i,i-nx)=(k*dx)/(2.0*dy);
        A(i,i-1)=(k*dy)/(2.0*dx);
        C(i)=-h_i*((dx/2.0)+(dy/2.0))*T_inf_i;
        B(i,i)=(- rho * cp * (0.7500 * dx * dy)) / dt;
        A(i,i)=-((k*dy)/dx)-((k*dx)/dy)-((k*dx)/(2.0*dy))-((k*dy)/(2.0*dx))-h_i*((dx/2.0)+(dy/2.0))+(- rho * cp * (0.7500 * dx * dy)) / dt;
       
    elseif i>start && i<start+(L/dx)
        A(i,i-1)=(k*dy)/(2.0*dx);
        A(i,i+1)=(k*dy)/(2.0*dx);
        A(i,i-nx)=(k*dx)/(dy);
        C(i)=-h_i*dx*T_inf_i;
        B(i,i)=(- rho * cp * (0.500 * dx * dy)) / dt;
        A(i,i)=-(k*dy)/(2.0*dx)-(k*dy)/(2.0*dx)-(k*dx)/(dy)-h_i*dx+(- rho * cp * (0.500 * dx * dy)) / dt;
        
    elseif i>start+((L/dy)*(nx)) && i<start+(L/dx)+((L/dy)*(nx))
        A(i,i-1)=(k*dy)/(2.0*dx);
        A(i,i+1)=(k*dy)/(2.0*dx);
        A(i,i+nx)=(k*dx)/(dy);
        C(i)=-h_i*dx*T_inf_i;
        B(i,i)=(- rho * cp * (0.500 * dx * dy)) / dt;
        A(i,i)=-(k*dy)/(2.0*dx)-(k*dy)/(2.0*dx)-(k*dx)/(dy)-h_i*dx+(- rho * cp * (0.500 * dx * dy)) / dt;
       
    elseif ismember(i,v1)
        A(i,i-1)=(k*dy)/dx;
        A(i,i-nx)=(k*dx)/(2.0*dy);
        A(i,i+nx)=(k*dx)/(2.0*dy);
        C(i)=-h_i*dy*T_inf_i;
        B(i,i)=(- rho * cp * (0.500 * dx * dy)) / dt;
        A(i,i)=-((k*dy)/dx)-((k*dx)/(2.0*dy))-((k*dx)/(2.0*dy))-(h_i*dy)+((- rho * cp * (0.500 * dx * dy)) / dt);
        
    elseif ismember(i,v2)
        A(i,i+1)=(k*dy)/dx;
        A(i,i-nx)=(k*dx)/(2.0*dy);
        A(i,i+nx)=(k*dx)/(2.0*dy);
        C(i)=-h_i*dy*T_inf_i;
        B(i,i)=(- rho * cp * (0.500 * dx * dy)) / dt;
        A(i,i)=-((k*dy)/dx)-((k*dx)/(2.0*dy))-((k*dx)/(2.0*dy))-(h_i*dy)+((- rho * cp * (0.500 * dx * dy)) / dt);
        
    else
        A(i,i-1)=(k*dy)/dx;
        A(i,i+1)=(k*dy)/dx;
        A(i,i-nx)=(k*dx)/dy;
        A(i,i+nx)=(k*dx)/dy;
        B(i,i)=(- rho * cp * (dx * dy)) / dt;
        A(i,i)=-(2*(k*dy)/dx)-(2*(k*dx)/dy)+(- rho * cp * (dx * dy)) / dt;
        
    end
end

for j = 1:(nTimes-1)
    T(:,j+1)=(A\B*T(:,j))+(A\C);
end

figure
% plot(T(start,:),'r');
% hold on;
% 
% plot(T(start+(L/dx),:),'g');
% hold on;
% 
% plot(T(start+((L/dy)*(nx)),:),'b');
% hold on;
% 
% plot(T(start+(L/dx)+((L/dy)*(nx)),:),'r*');
% hold on;
% 
% plot(T(nTemps,:));
% hold on;

for o=1:nTemps
    plot(T(o,:));
    hold on;
end

uu=[];
for u=1:nTimes
    if abs(T(1,u)-T(floor(start-((0.020/dy)*nx)),u))>5
        uu=[uu,u];
    end
    
end
Tmat = reshape(T(:,45), [nx,ny]);
Tplot_r = Tmat.';

% Do post-processing/ plotting here
figure;
imagesc(Tplot_r);   
colorbar;
% At steady state, verify the energy balance equation, and determine the
% total heat transferred to the process fluid=s
s=0;
for r=370:374
    s=s+(h_i*(T_inf_i-T(r,61)));
end
for r=734:739
    s=s+(h_i*(T_inf_i-T(r,61)));
end
for d=1:nTemps
   if ismember(d,v1) || ismember(d,v2)
       s=s+(h_i*(T_inf_i-T(d,61)));
       
   end
end