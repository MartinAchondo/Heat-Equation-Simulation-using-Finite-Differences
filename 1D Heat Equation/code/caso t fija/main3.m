
% Programa principal
% Caso temperatura fija en interior

format compact

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

dx = 0.01;  %% 0.005   0.01
dt = 1;   %1   10

t = 0;
hrs = 24;
hrsi = 0;
tf = hrs*3600;


N1 = floor(L1/dx);
N2 = floor(L2/dx);
N = N1 + N2 + 1;
Nt = hrs + 1;

Ti = ones(1,N)*Tw;
L = linspace(0,L1+L2,N);
k = N1 + 1;

dt2 = 30;
if dt2<dt
    dt2 = dt;
end
Nt2f = hrs*3600/dt2;
Nt2i = (hrs-24)*3600/dt2+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%


[T_out1,TN] = Euler_Explicito(Ti,N,Nt,t,k,dx,dt,tf,prob,prop,dt2,Nt2f,Nt2i);
f1 = figure(1);
plot_hrs(L,T_out1,hrsi,hrs)
hold on

T_out = Stationary(prob,prop,L1,L1+L2,k,N);
plot(L,T_out,'DisplayName','Stationary')
hold on
legend('Location', 'Best')
set(f1, 'Position', [150 90 650 530]);



%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T_out] = Stationary(prob,prop,L1,L,k,N)
    
    r = num2cell(prop);
    [~,~,lambda_1,lambda_2,h] = deal(r{:});
    r = num2cell(prob);
    [Tw,~,T_inf] = deal(r{:});

    x = linspace(0,L,N);

    phi = h*(L*lambda_1-L1*(lambda_1-lambda_2))+lambda_1*lambda_2;
    c1 = h*lambda_2*(T_inf-Tw)/phi;
    c2 = h*lambda_1*(T_inf-Tw)/phi;
    c3 = (h*(L*lambda_1*Tw-L1*(lambda_1-lambda_2)*T_inf)+lambda_1*lambda_2*Tw)/phi;

    T_out = zeros(1,N);
    T_out(1:k) = c1*x(1:k) + Tw;
    T_out(k+1:N) = c2*x(k+1:N) + c3;


end


function [T_out,TN_out] = Euler_Explicito(Ti,N,Nt,t,k,dx,dt,tf,prob,prop,dt2,Nt2f,Nt2i)

    r = num2cell(prop);
    [alpha_1,alpha_2,lambda_1,lambda_2,h] = deal(r{:});
    r = num2cell(prob);
    [Tw,dTw,T_inf] = deal(r{:});

    T = zeros(1,N);
    T_old = Ti;
    T_out = zeros(Nt,N);
    TN = zeros(1,Nt2f);
    z = 1;
    zz = 1;
    T_out(z,:) = Ti(:);

    while t<tf
        t = t + dt;
        for i=1:N
            if i==1
                T(i) = Tw + dTw*sin((2*pi*t/3600)/(24));
                %T(i) = Tw;
            elseif i>1 && i<k
                T(i) = T_old(i) + (dt/(dx^2))*(T_old(i+1)-2*T_old(i)+T_old(i-1))*alpha_1;
            elseif i>k && i<N
                T(i) = T_old(i) + (dt/(dx^2))*(T_old(i+1)-2*T_old(i)+T_old(i-1))*alpha_2;
            elseif i==N
                T(i) = T_inf;
                %T(i) = (lambda_2*T(i-1)/dx + h*T_inf)/(h+lambda_2/dx);
            end
        end
        T(k) = (lambda_2*T(k+1)+lambda_1*T(k-1))/(lambda_1+lambda_2);
        T_old(:) = T(:);
        if mod(t,3600)==0
            z = z + 1;
            T_out(z,:) = T(:);
        end
        if mod(t,dt2) == 0
            zz = zz + 1;
            TN(zz) = T(N);
        end
    end
    TN_out = TN(Nt2i:Nt2f);
end


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
