

% Calculo del flujo de calor en interior 

function [q,qj] = Heat_flux_int(T,dt,prob,prop)
    
    r = num2cell(prop);
    [~,~,~,~,h] = deal(r{:});
    r = num2cell(prob);
    [~,~,T_inf] = deal(r{:});
    qj(:) = h*(T(:)-T_inf);

    q = Total_heat(qj,dt);

    function [qt] = Total_heat(qj,dt)
        sum = 0;
        for j=1:length(qj)-1
            sum = sum + qj(j) + qj(j+1);
        end
        qt = sum*dt/2;
    end

end
