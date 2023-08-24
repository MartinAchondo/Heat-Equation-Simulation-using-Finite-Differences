% Perfil de temperatura estacionaria teorica

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