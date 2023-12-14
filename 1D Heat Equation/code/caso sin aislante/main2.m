
% Programa principal
% Caso sin aislante

format compact

aisl_in = true;


lambda_1 = 2;
pC_1 = 2*10^6;
L1 = 0.45;



alpha_1 = lambda_1/pC_1;
h = 4;
prop = [alpha_1,lambda_1,h];

Tw = 10;
dTw = 5;
T_inf = 20;
prob = [Tw,dTw,T_inf];

dx = 0.01;  %% 0.005   0.01
dt = 1;

t = 0;
hrs = 394;
hrsi = 370;
tf = hrs*3600;


N1 = floor(L1/dx);
N = N1 + 1;
Nt = hrs + 1;

Ti = ones(1,N)*Tw;
L = linspace(0,L1,N);

dt2 = 30;
Nt2f = hrs*3600/dt2;
Nt2i = (hrs-24)*3600/dt2+1;


%%%%%%%%%%%%%%%%%

[T_out1,TN] = Euler_Explicito(Ti,N,Nt,t,dx,dt,tf,prob,prop,dt2,Nt2f,Nt2i);
plot_hrs(L,T_out1,hrsi,hrs)
hold on

T_out = Stationary(prob,prop,L1,N);
plot(L,T_out)
hold on


[q,qj] = Heat_flux_int(TN,dt2,prob,prop);


% Amplitud oscilante
figure(5)
Tsin = T_out1(hrs-24:hrs,:);
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

function [q,qj] = Heat_flux_int(T,dt,prob,prop)
    
    r = num2cell(prop);
    [~,~,h] = deal(r{:});
    r = num2cell(prob);
    [~,~,T_inf] = deal(r{:});
    qj(:) = -h*(T(:)-T_inf);

    q = Total_heat(qj,dt);

    function [qt] = Total_heat(qj,dt)
        sum = 0;
        for j=1:length(qj)-1
            sum = sum + qj(j) + qj(j+1);
        end
        qt = sum*dt/2;
    end

end

function [T_out] = Stationary(prob,prop,L,N)
    
    r = num2cell(prop);
    [~,lambda_1,h] = deal(r{:});
    r = num2cell(prob);
    [Tw,~,T_inf] = deal(r{:});

    x = linspace(0,L,N);

    c1 = h*(T_inf-Tw)/(lambda_1+h*L);

    T_out = zeros(1,N);
    T_out(:) = c1*x(:) + Tw;


end


function [T_out,TN_out] = Euler_Explicito(Ti,N,Nt,t,dx,dt,tf,prob,prop,dt2,Nt2f,Nt2i)

    r = num2cell(prop);
    [alpha_1,lambda_1,h] = deal(r{:});
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
            elseif i>1 && i<N
                T(i) = T_old(i) + (dt/(dx^2))*(T_old(i+1)-2*T_old(i)+T_old(i-1))*alpha_1;
            elseif i==N
                %T(i) = T_inf;
                T(i) = (lambda_1*T(i-1)/dx + h*T_inf)/(h+lambda_1/dx);
            end
        end
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