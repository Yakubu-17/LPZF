function [ power] = PowerControl(gamma_matrix, p_max)
[rows, columns]=size(gamma_matrix);
LM=rows;
K=columns;
power = zeros(LM,K);
for m=1:LM
    for k=1:K
        power(m,k) = (gamma_matrix(m,k)/sum(gamma_matrix(m,:)))*p_max; 
    end
end

%power = p_max * ones(LM, K);    
    %Normalize gamma for each AP to get the power control values
    %power1= gamma_matrix ./ sum(gamma_matrix, 2) * p_max;
    %size(power1)
end