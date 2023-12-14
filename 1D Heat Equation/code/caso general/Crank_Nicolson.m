
% Crank Nicolson

function [T_out,TN_out] = Crank_Nicolson(Ti,N,Nt,t,k,dx,dt,tf,prob,prop,dt2,Nt2f,Nt2i)

    r = num2cell(prop);
    [alpha_1,alpha_2,lambda_1,lambda_2,h] = deal(r{:});
    r = num2cell(prob);
    [Tw,dTw,T_inf] = deal(r{:});

    T_old(:) = Ti(:);
    T_out = zeros(Nt,N);
    TN = zeros(1,Nt2f);
    z = 1;
    zzz = 1;
    T_out(z,:) = Ti(:);
    
    sigma_1 = alpha_1*dt/(dx^2);
    sigma_2 = alpha_2*dt/(dx^2);
    theta = lambda_2/(2*dx);
    zz = 0;
    if mod(N,3) == 2
        zz = 1;
        M = floor(N/3) + zz;
    elseif mod(N,3) == 0
        M = floor(N/3);
    end
    kk = floor(k/3);
    if mod(k,3)==0
        kk = kk - 1;
    end

    A = zeros(M,3,3);
    B = zeros(M,3,3);
    C = zeros(M,3,3);
    F = zeros(M,3);
    for i=1:M
        if i==1
            C(i,:,:) = [1 0 0; -1 2+2/sigma_1 -1; 0 -1 2+2/sigma_1];
            B(i,:,:) = [0 0 0; 0 0 0; -1 0 0];            
        elseif i<=kk
            A(i,:,:) = [0 0 -1; 0 0 0; 0 0 0]; 
            C(i,:,:) = [2+2/sigma_1 -1 0; -1 2+2/sigma_1 -1; 0 -1 2+2/sigma_1];
            B(i,:,:) = [0 0 0; 0 0 0; -1 0 0]; 
        elseif i==kk+1
            if mod(k,3)==0
                A(i,:,:) = [0 0 -1; 0 0 0; 0 0 0]; 
                C(i,:,:) = [2+2/sigma_1 -1 0; -1 2+2/sigma_1 -1; lambda_1 -4*lambda_1 3*(lambda_1+lambda_2)];
                B(i,:,:) = [0 0 0; 0 0 0; -4*lambda_2 lambda_2 0];
            else
                A(i,:,:) = [0 lambda_1 -4*lambda_1; 0 0 0; 0 0 0]; 
                C(i,:,:) = [3*(lambda_1+lambda_2) -4*lambda_2 lambda_2; -1 2+2/sigma_2 -1; 0 -1 2+2/sigma_2];
                B(i,:,:) = [0 0 0; 0 0 0; -1 0 0];
            end
        elseif i>kk+1 && i<M
            A(i,:,:) = [0 0 -1; 0 0 0; 0 0 0];
            C(i,:,:) = [2+2/sigma_2 -1 0; -1 2+2/sigma_2 -1; 0 -1 2+2/sigma_2];
            B(i,:,:) = [0 0 0; 0 0 0; -1 0 0];  
        elseif i==M
            if mod(N,3) == 2
                A(i,:,:) = [0 0 -1; 0 0 theta; 0 0 0];
                C(i,:,:) = [2+2/sigma_2 -1 0; -4*theta 3*theta+h 0; 0 0 1];
            elseif mod(N,3) == 0
                A(i,:,:) = [0 0 -1; 0 0 0; 0 0 0];
                C(i,:,:) = [2+2/sigma_2 -1 0; -1 2+2/sigma_2 -1; theta -4*theta 3*theta+h];
            end
        end 
    end
    
    while t<tf
        t = t + dt;
        for i=1:M
            if i==1
                F(i,1) = Tw + dTw*sin((2*pi*t/3600)/(24));
                F(i,2) = T_old(i*3-2)+(2/sigma_1-2)*T_old(i*3-1)+T_old(i*3);
                F(i,3) = T_old(i*3-1)+(2/sigma_1-2)*T_old(i*3)+T_old(i*3+1);               
            elseif i<=kk
                F(i,1) = T_old(i*3-3)+(2/sigma_1-2)*T_old(i*3-2)+T_old(i*3-1);
                F(i,2) = T_old(i*3-2)+(2/sigma_1-2)*T_old(i*3-1)+T_old(i*3);
                F(i,3) = T_old(i*3-1)+(2/sigma_1-2)*T_old(i*3)+T_old(i*3+1);
            elseif i==kk+1
                if mod(k,3)==0
                    F(i,1) = T_old(i*3-3)+(2/sigma_1-2)*T_old(i*3-2)+T_old(i*3-1);
                    F(i,2) = T_old(i*3-2)+(2/sigma_1-2)*T_old(i*3-1)+T_old(i*3);
                    F(i,3) = 0;
                else
                    F(i,1) = 0;
                    F(i,2) = T_old(i*3-2)+(2/sigma_2-2)*T_old(i*3-1)+T_old(i*3);
                    F(i,3) = T_old(i*3-1)+(2/sigma_2-2)*T_old(i*3)+T_old(i*3+1);
                end
            elseif i>kk+1 && i<M
                F(i,1) = T_old(i*3-3)+(2/sigma_2-2)*T_old(i*3-2)+T_old(i*3-1);
                F(i,2) = T_old(i*3-2)+(2/sigma_2-2)*T_old(i*3-1)+T_old(i*3);
                F(i,3) = T_old(i*3-1)+(2/sigma_2-2)*T_old(i*3)+T_old(i*3+1);   
            elseif i==M
                if mod(N,3) == 2
                    F(i,1) = T_old(i*3-3)+(2/sigma_2-2)*T_old(i*3-2)+T_old(i*3-1);
                    F(i,2) =  h*T_inf;
                    F(i,3) = 1;
                elseif mod(N,3) == 0
                    F(i,1) = T_old(i*3-3)+(2/sigma_2-2)*T_old(i*3-2)+T_old(i*3-1);
                    F(i,2) = T_old(i*3-2)+(2/sigma_2-2)*T_old(i*3-1)+T_old(i*3);
                    F(i,3) = h*T_inf;
                end
            end 
        end
        T_star = Thomas_algorithm(A,B,C,F,M);
        T_old(:) = T_star(1:M*3-zz);
        if mod(t,3600)==0
            z = z + 1;
            T_out(z,:) = T_old(:);
        end
        if mod(t,dt2) == 0
            zzz = zzz + 1;
            TN(zzz) = T_old(N);
        end
    end
    TN_out = TN(Nt2i:Nt2f);
end
