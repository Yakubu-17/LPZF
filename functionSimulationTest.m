%This function returns the sum rate for MRT and Full pilot zero forcing.
function[sumrate_MRT,sumrate_FZF] = functionSimulationTest(L,M,K,p_maxdB,ppdB,beta_matrix,Phii,tau,Pilot_indices,pilot_book,experiment)
%Converting p_maxdB  and pilot power to watt
p_max=10^(p_maxdB/10); pp=10^(ppdB/10);


nu=zeros(L,K);gamma_matrix=zeros(L,K);eta=zeros(L,K);etazf=zeros(L,K);
% Solve for lambda and nu: the variance of the channel estimates
for m=1:L
   
    for idx=1:K
        d=0;
        for iidx=1:K
            d=d+beta_matrix(m,iidx)*norm(Phii(:,idx)'*Phii(:,iidx));
        end
        nu(m,idx)=1 + tau*pp*d;
        gamma_matrix(m,idx)=pp*tau*beta_matrix(m,idx)^2*(1/nu(m,idx));
    end    
end
%power normalization factor
for n=1:L
   for kk=1:K
      eta(n,kk)= 1/(M*sum(gamma_matrix(n,:)));  
   end
end

for n=1:L
   for kk=1:K
      etazf(n,kk)= ((M-tau)*sum(gamma_matrix(n,:)));  
   end
end


%Power control fraction
% p_control = zeros(L,K);
% for mm=1:L
%     for k=1:K
%         p_control(mm,k) = (gamma_matrix(mm,k)/sum(gamma_matrix(mm,:))); 
%         %p_control(mm,k) = (gamma_matrix(mm,k) / sum(gamma_matrix(:, k)));
%     end
% end

Eta=repelem(eta,M,1);
Etazf=repelem(etazf,M,1);
Beta = repelem(beta_matrix,M,1);
%P_control = repelem(p_control,M,1);
%Gamma = repelem(gamma_matrix,M,1);
H = (randn(L*M,K,experiment) + 1i*randn(L*M,K,experiment))*sqrt(0.5); % channel


%power =PowerControl(Gamma, p_max);

%Storage Terms
Desired_MRT = zeros(K,experiment); %desired signal 
MUI_MRT = zeros(K,experiment); %MUG + BUG
Desired_FZF = zeros(K,experiment); %desired signal 
MUI_FZF = zeros(K,experiment); %MUG + BUG
alpha = 1;
for expr=1:experiment
    H_matrix = H(:,:,expr).*sqrt(Beta); 
    %[~,~,Hhat] = Estimator(pp,tau,Phii,H_matrix,Beta,pilot_book,Pilot_indices); %Estimate H_bar matrix and Hhat
    %[~,~,Hhat] = test_CE(pp,tau,Phii,H_matrix,Beta,pilot_book,Pilot_indices); %Estimate H_bar matrix and Hhat
    [Hhat,Hzf] = ChannelEstimator(tau,pp,H_matrix,Beta,Phii,alpha);%estimate the channel
     for jx = 1:K
        power1 = diag(Eta(:,jx));
        power1zf = diag(Etazf(:,jx));
        %p1 = diag(power(:,jx));%power control fraction.
        Desired_MRT(jx,expr) = (H_matrix(:,jx)'*sqrt(power1)*Hhat(:,jx));%Desired signal for MRT. 
        Desired_FZF(jx,expr) = (H_matrix(:,jx)'*sqrt(power1zf)*Hzf(:,jx));%Desired signal for MRT.
        


        %Computation for power MUI,BUG and pilot contamination.
        for jxx=1:K
            power2 = diag(Eta(:,jxx));
            power2zf = diag(Etazf(:,jxx)); 
            %p2 = diag(power(:,jxx));%power control fraction.
            MUI_MRT(jx,expr) = MUI_MRT(jx,expr) + abs(H_matrix(:,jx)'*sqrt(power2)*Hhat(:,jxx)).^2;
            MUI_FZF(jx,expr) = MUI_FZF(jx,expr) + abs(H_matrix(:,jx)'*sqrt(power2zf)*Hzf(:,jxx)).^2;
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

    sinr_MRT  = p_max*Desiredx_MRT(kj)^2/(p_max*(MUIx_MRT(kj)-Desiredx_MRT(kj)^2) + 1 );
    rate_MRT(kj)= log2(1 + sinr_MRT );
end

rate_FZF=zeros(1,K);

for kj=1:K

    sinr_FZF  = p_max*Desiredx_FZF(kj)^2/(p_max*(MUIx_FZF(kj)-Desiredx_FZF(kj)^2) + 1 );
    rate_FZF(kj)= log2(1 + sinr_FZF );
end



sumrate_MRT = sum(rate_MRT);
sumrate_FZF = sum(rate_FZF); 


end