
% Algoritmo de thomas por bloques

function [X_out] = Thomas(A,B,C,F,n)

    psi_1 = zeros(n,3,3);
    psi_2 = zeros(n,3);
    X = zeros(n,3);
    X_out = zeros(1,n*3);
    
    psi_1(1,:,:) = (to_mat(C,1))\to_mat(B,1);
    psi_2(1,:) = transpose((to_mat(C,1))\to_vec(F,1)');
    

    for i=2:n-1
        psi_1(i,:,:) = (to_mat(C,i)-to_mat(A,i)*to_mat(psi_1,i-1))\to_mat(B,i);
        psi_2(i,:) = transpose((to_mat(C,i)-to_mat(A,i)*to_mat(psi_1,i-1))\transpose(to_vec(F,i)-transpose(to_mat(A,i)*to_vec(psi_2,i-1)')));
    end
    psi_2(n,:) = transpose((to_mat(C,n)-to_mat(A,n)*to_mat(psi_1,n-1))\transpose(to_vec(F,n)-transpose(to_mat(A,n)*to_vec(psi_2,n-1)')));

    X(n,:) = to_vec(psi_2,n);
    for ss=1:n-1
        i = n-ss;
        X(i,:) = to_vec(psi_2,i) - transpose(to_mat(psi_1,i)*to_vec(X,i+1)'); 
    end

    for i=1:n
        X_out(i*3-2) = X(i,1);
        X_out(i*3-1) = X(i,2);
        X_out(i*3) = X(i,3);
    end


    function [Z] = to_mat(X,q)
        Z = zeros(3,3);
        for l=1:3
            for j=1:3
                Z(l,j) = X(q,l,j);
            end
        end
    end

    function [Z] = to_vec(X,q)
        Z = zeros(1,3);
        for l=1:3
            Z(l) = X(q,l);
        end
    end

end
