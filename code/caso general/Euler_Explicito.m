

% Euler explicito

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
                %T(i) = T_inf;
                T(i) = (lambda_2*T(i-1)/dx + h*T_inf)/(h+lambda_2/dx);
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

