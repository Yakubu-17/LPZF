%This function returns the sum rate for MRT and Full pilot zero forcing.
function[sumrate_MRT, sumrate_FZF] = functionSimulationPU(L,M,K,p_maxdB,ppdB,beta_matrix,Phii,tau,Pilot_indices,pilot_book,experiment)
%Converting p_maxdB  and pilot power to watt
p_max=10^(p_maxdB/10); pp=10^(ppdB/10);


%Solve for gamma,theta and c matrix.
nu=zeros(L,K);gamma_matrix=zeros(L,K);
theta_matrix=zeros(L,K); c_matrix=zeros(L,K);
for l=1:L  
    for kx=1:K
        d=0;
        for iidx=1:K
            d=d+ beta_matrix(l,iidx)*norm(Phii(:,kx)'*Phii(:,iidx));
        end
        nu(l,kx)=1 + tau*pp*d;
        gamma_matrix(l,kx)=pp*tau*beta_matrix(l,kx)^2*(1/nu(l,kx));
        c_matrix(l,kx) = sqrt(pp)*beta_matrix(l,kx)*(1/nu(l,kx));
        theta_matrix(l,kx) = gamma_matrix(l,kx)/(c_matrix(l,kx)^2);
    end    
end


%normalization factor MRT
p_norm_MRT=zeros(L,K);
for n=1:L
   for kk=1:K
      p_norm_MRT(n,kk)= 1/(M*theta_matrix(n,kk)); 
   end
end


%normalization factor FZF
p_norm_FZF=zeros(L,K);
for nn=1:L
   for kc=1:K
      p_norm_FZF(nn,kc)= ((M-tau)*theta_matrix(nn,kc)); 
   end
end

%Power control fraction
p_control = zeros(L,K);
for mm=1:L
    for k=1:K
        p_control(mm,k) = (gamma_matrix(mm,k)/sum(gamma_matrix(mm,:))); 
        %p_control(mm,k) = (gamma_matrix(mm,k) / sum(gamma_matrix(:, k)));
    end
end





P_norm_MRT=repelem(p_norm_MRT,M,1);
Beta = repelem(beta_matrix,M,1);
%Gamma = repelem(gamma_matrix,M,1);
P_control = repelem(p_control,M,1);
P_norm_FZF=repelem(p_norm_FZF,M,1);

H = (randn(L*M,K,experiment) + 1i*randn(L*M,K,experiment))*sqrt(0.5); % channel

%Storage Terms
Desired_MRT = zeros(K,experiment); %desired signal 
MUI_MRT = zeros(K,experiment); %MUG + BUG
Desired_FZF = zeros(K,experiment); %desired signal for FZF
MUI_FZF = zeros(K,experiment); %MUG + BUG for FZF


%power =PowerControl(Gamma, p_max); %Estimating power after applying power control.

for expr=1:experiment
H_matrix = H(:,:,expr).*sqrt(Beta); 
[Hest_MRT,Hest_FZF] = Estimator(pp,tau,Phii,H_matrix,Beta,pilot_book,Pilot_indices); %Estimate H_bar matrix
Wprecoder_MRT=Hest_MRT.*sqrt(P_norm_MRT);
Wprecoder_FZF=Hest_FZF.*sqrt(P_norm_FZF);
    for jx = 1:K
        power1 = diag(P_control(:,jx));%power control fraction.
        Desired_MRT(jx,expr) = abs(H_matrix(:,jx)'*sqrt(power1 )*Wprecoder_MRT(:,jx));%Desired signal for MRT. 
        Desired_FZF(jx,expr) = abs(H_matrix(:,jx)'*sqrt(power1 )*Wprecoder_FZF(:,jx));%Desired signal for FZF. 


        %Computation for power MUI,BUG and pilot contamination.
        for jxx=1:K
            power2 = diag(P_control(:,jxx)); %power control fraction.
            MUI_MRT(jx,expr) = MUI_MRT(jx,expr) + abs(H_matrix(:,jx)'*power2*Wprecoder_MRT(:,jxx)).^2;
            MUI_FZF(jx,expr) = MUI_FZF(jx,expr) + abs(H_matrix(:,jx)'*power2*Wprecoder_FZF(:,jxx)).^2;
        end
       
    end
    
end

Desiredx_MRT = real(mean(Desired_MRT,2));
MUIx_MRT = mean(MUI_MRT,2);
Desiredx_FZF = real(mean(Desired_FZF,2));
MUIx_FZF = mean(MUI_FZF,2);

%Computing for MRT rate
rate_MRT=zeros(1,K);

for kj=1:K

    sinr_MRT  = p_max*Desiredx_MRT(kj)^2/( p_max*(MUIx_MRT(kj)-Desiredx_MRT(kj)^2) + 1 );
    rate_MRT(kj)= log2(1 + sinr_MRT );
end


%Computing for FZF rate
rate_FZF=zeros(1,K);
for dj=1:K

    sinr_FZF = p_max*Desiredx_FZF(dj)^2/( p_max*(MUIx_FZF(dj)-Desiredx_FZF(dj)^2) + 1 ); 
    rate_FZF(dj)= log2(1 + sinr_FZF);
end



sumrate_MRT = sum(rate_MRT)/K; 
sumrate_FZF = sum(rate_FZF)/K; 
end