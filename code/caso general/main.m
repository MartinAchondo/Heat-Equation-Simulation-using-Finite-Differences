
% Programa principal EDP Calor
% Caso convectivo en interior

format compact

% Modificacion aislante interior o exterior
aisl_in = true;

if aisl_in
    lambda_1 = 2;
    pC_1 = 2*10^6;
    lambda_2 = 0.1;
    pC_2 = 10^6;
    L1 = 0.45;
    L2 = 0.05;
else
    lambda_2 = 2;
    pC_2 = 2*10^6;
    lambda_1 = 0.1;
    pC_1 = 10^6;
    L2 = 0.45;
    L1 = 0.05;
end

alpha_1 = lambda_1/pC_1;
alpha_2 = lambda_2/pC_2;
h = 4;
prop = [alpha_1,alpha_2,lambda_1,lambda_2,h];

Tw = 10;
dTw = 5;
T_inf = 20;
prob = [Tw,dTw,T_inf];

% parametros numericos

dx = 0.01;  %% 0.005   0.01
dt = 10;   %1   10 15

t = 0;
hrs = 49;
hrsi = 25;
tf = hrs*3600;

N1 = floor(L1/dx);
N2 = floor(L2/dx);
N = N1 + N2 + 1;
Nt = hrs + 1;

Ti = ones(1,N)*Tw;
L = linspace(0,L1+L2,N);
k = N1 + 1;

% outputs para figuras en pared interior

dt2 = 30;
if dt2<dt
    dt2 = dt;
end
Nt2f = hrs*3600/dt2;
Nt2i = (hrs-24)*3600/dt2+1;


%%%%% EULER EXPLICITO %%%%%%%

[T_out1,TN] = Euler_Explicito(Ti,N,Nt,t,k,dx,dt,tf,prob,prop,dt2,Nt2f,Nt2i);
f1 = figure(1);
plot_hrs(L,T_out1,hrsi,hrs)
hold on

[q,qj] = Heat_flux_int(TN,dt2,prob,prop);

T_out = Stationary(prob,prop,L1,L1+L2,k,N);
plot(L,T_out,'DisplayName','Stationary')
hold on
legend('Location', 'Best')
set(f1, 'Position', [150 90 650 530]);

%%%%% EULER IMPLICITO %%%%%%

[T_out2,TN] = Euler_Implicito(Ti,N,Nt,t,k,dx,dt,tf,prob,prop,dt2,Nt2f,Nt2i);
f2 = figure(2);
plot_hrs(L,T_out2,hrsi,hrs)
hold on

[q,qj] = Heat_flux_int(TN,dt2,prob,prop);

T_out = Stationary(prob,prop,L1,L1+L2,k,N);
plot(L,T_out,'DisplayName','Stationary')
hold on
legend('Location', 'Best')
set(f2, 'Position', [150 90 650 530]);

%%%%% CRANK NICOLSON %%%%%%%


[T_out3,TN] = Crank_Nicolson(Ti,N,Nt,t,k,dx,dt,tf,prob,prop,dt2,Nt2f,Nt2i);
f3 = figure(3);
plot_hrs(L,T_out3,hrsi,hrs)
hold on

[q,qj] = Heat_flux_int(TN,dt2,prob,prop);

T_out = Stationary(prob,prop,L1,L1+L2,k,N);
plot(L,T_out,'DisplayName','Stationary')
hold on
legend('Location', 'Best')
set(f3, 'Position', [150 90 650 530]);

%%%%%% GRAFICO CALOR EN CARA INTERIOR %%%%%

figure(4)
t = linspace(0,24,length(qj));
plot(t,qj,'DisplayName','Flujo calor cara interna 24 hrs') 
ylabel('q [W/m2]')
xlabel('t [hrs]')
xlim([0 24])
legend

%%%%% GRAFICO TEMPERATURA EN  INTERIOR %%%%%

figure(5)
plot(t,TN,'DisplayName','Temperatura cara iterna') 
ylabel('T [°C]')
xlabel('t [hrs]')
xlim([0 24])
legend

%%%%% AMPLITUD TMPERATURA OSCILANTE %%%%

figure(6)
Tsin = T_out3(hrs-24:hrs,:);
T_tot = zeros(1,N);
for j=1:N
    xx = max(Tsin(:,j));
    T_tot(j) = xx;
end
plot(L,T_tot-T_out)
xlabel('L [m]')
ylabel('T [°C]')
legend('Amplitud oscilante A(x)') 



%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Funcion plots

function plot_hrs(L,T,ti,tf)
    if ti==0
        ti=1;
        tf = tf + 1;
    end
    for j=ti:tf
        %a = mod(j,24);
        % if a==0
        %     a=24
        % elseif a==1
        %     a=25
        % end
        plot(L,T(j,:),'DisplayName',strcat(int2str(j-1),'hrs'))
        hold on
    end
    xlabel('L [m]')
    ylabel('T [°C]')
    legend
    ylim([4 21])
end
