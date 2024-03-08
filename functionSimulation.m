%This function returns the sum rate for MRT and Full pilot zero forcing.
function[sumrate_MRT, sumrate_FZF,sumrate_PZF,sumrate_PPZF] = functionSimulation(L,M,K,p_maxdB,ppdB,beta_matrix,Phii,tau,Pilot_indices,pilot_book,S,W,Z,X,E_S_temp,tau_S_l,experiment)


p_max=10^(p_maxdB/10); pp=10^(ppdB/10);


%Solve for gamma,theta and c matrix.
nu=zeros(L,K);gamma_matrix=zeros(L,K);
theta_matrix=zeros(L,K); c_matrix=zeros(L,K);
for l=1:L
    for kx=1:K
        d=0;
        for iidx=1:K
            d=d + beta_matrix(l,iidx)*norm(Phii(:,kx)'*Phii(:,iidx));
        end
        nu(l,kx)= tau*pp*d + 1;
        gamma_matrix(l,kx)=pp*tau*beta_matrix(l,kx)^2*(1/nu(l,kx));
        c_matrix(l,kx) = sqrt(pp)*beta_matrix(l,kx)*(1/nu(l,kx));
        theta_matrix(l,kx) = gamma_matrix(l,kx)/(c_matrix(l,kx)^2);
    end
end





%normalization factor for MRT
p_norm_MRT=zeros(L,K);
for n=1:L
    for kk=1:K
        p_norm_MRT(n,kk)= sqrt(M*theta_matrix(n,kk));
    end
end


%normalization factor FZF
p_norm_FZF=zeros(L,K);
for nn=1:L
    for kc=1:K
        p_norm_FZF(nn,kc)= sqrt(1/((M-tau)*theta_matrix(nn,kc)));
    end
end

%normalization factor PZF
p_norm_PZF=zeros(L,K);
for n=1:L
    for k=1:K
        p_norm_PZF(n,k)= sqrt(1/((M-tau_S_l)*theta_matrix(n,k)));
    end
end



%normalization factor PPZF_MRT
p_norm_PPZF_MRT=zeros(L,K);
for n=1:L
    for kk=1:K
        p_norm_PPZF_MRT(n,kk)= sqrt((M-tau_S_l)*theta_matrix(n,kk));
    end
end




% p_control = ones(L,K);

P_norm_MRT=repelem(p_norm_MRT,M,1);
Beta = repelem(beta_matrix,M,1);
Gamma = repelem(gamma_matrix,M,1);
P_norm_FZF=repelem(p_norm_FZF,M,1);
P_norm_PZF=repelem(p_norm_PZF,M,1);
P_norm_PPZF_MRT=repelem(p_norm_PPZF_MRT,M,1);


H = (randn(L*M,K,experiment) + 1i*randn(L*M,K,experiment))*sqrt(0.5); % channel

%Storage Terms
Desired_MRT = zeros(K,experiment); %desired signal
MUI_MRT = zeros(K,experiment); %MUG + BUG
Desired_FZF = zeros(K,experiment); %desired signal for FZF
MUI_FZF = zeros(K,experiment); %MUG + BUG for FZF
Desired_PZF = zeros(K,experiment); %desired signal for PZF
MUI_PZF = zeros(K,experiment); %MUG + BUG for PZF
Desired_PZF_MRT = zeros(K,experiment); %desired signal for PZF_MRT
MUI_PZF_MRT = zeros(K,experiment); %MUG + BUG for PZF_MRT
Desired_PPZF_MRT = zeros(K,experiment); %desired signal for PPZF_MRT
MUI_PPZF_MRT = zeros(K,experiment); %MUG + BUG for PPZF_MRT





power =PowerControl(Gamma, p_max); %Estimating power after applying power control.


for expr=1:experiment
    H_matrix = H(:,:,expr).*sqrt(Beta);
    [Hest_MRT,Hest_FZF,Hest_PZF,Hest_PZF_MRT,Hest_PPZF_MRT] = ChannelEstimator1(tau,pp,H_matrix,M,Phii,pilot_book,Pilot_indices,S,W,Z,X,E_S_temp);%estimate the channel

    %Forming precoders.
    Wprecoder_MRT=Hest_MRT./P_norm_MRT;
    Wprecoder_FZF=Hest_FZF./P_norm_FZF;
    Wprecoder_PZF=Hest_PZF./P_norm_PZF;
    Wprecoder_PZF_MRT = Hest_PZF_MRT./P_norm_MRT;
    Wprecoder_PPZF_MRT = Hest_PPZF_MRT./P_norm_PPZF_MRT;
    for jx = 1:K

        desired_MRT_temp = 0;
        desired_FZF_temp = 0;
        desired_PZF_temp = 0;
        desired_PZF_MRT_temp = 0;
        desired_PPZF_MRT_temp = 0;
        for l=1:L
            power1 = diag(power((l-1)*M+1:M*l,jx));%power control fraction.
            desired_MRT_temp =desired_MRT_temp + (H_matrix((l-1)*M+1:M*l,jx)'*sqrt(power1 )*Wprecoder_MRT((l-1)*M+1:M*l,jx));%Desired signal for MRT.
            desired_FZF_temp = desired_FZF_temp + (H_matrix((l-1)*M+1:M*l,jx)'*sqrt(power1 )*Wprecoder_FZF((l-1)*M+1:M*l,jx));%Desired signal for FZF.
            desired_PZF_temp =desired_PZF_temp + (H_matrix((l-1)*M+1:M*l,jx)'*sqrt(power1 )*Wprecoder_PZF((l-1)*M+1:M*l,jx));%Desired signal for PZF.
            desired_PZF_MRT_temp =desired_PZF_MRT_temp + (H_matrix((l-1)*M+1:M*l,jx)'*sqrt(power1 )*Wprecoder_PZF_MRT((l-1)*M+1:M*l,jx));%Desired signal for PZF_MRT.
            desired_PPZF_MRT_temp = desired_PPZF_MRT_temp + (H_matrix((l-1)*M+1:M*l,jx)'*sqrt(power1 )*Wprecoder_PPZF_MRT((l-1)*M+1:M*l,jx));%Desired signal for PPZF_MRT.
        end
        Desired_MRT(jx,expr)=desired_MRT_temp;
        Desired_FZF(jx,expr)=desired_FZF_temp;
        Desired_PZF(jx,expr)=desired_PZF_temp;
        Desired_PZF_MRT(jx,expr)=desired_PZF_MRT_temp;
        Desired_PPZF_MRT(jx,expr)=desired_PPZF_MRT_temp;

        %Computation for power MUI,BUG and pilot contamination.
        for jxx=1:K
            mui_MRT_temp=0;
            mui_FZF_temp=0;
            mui_PZF_temp=0;
            mui_PZF_MRT_temp=0;
            mui_PPZF_MRT_temp=0;

            for l=1:L
                power2 = diag(power((l-1)*M+1:M*l,jxx)); %power control fraction.
                mui_MRT_temp = mui_MRT_temp + (H_matrix((l-1)*M+1:M*l,jx)'*sqrt(power2)*Wprecoder_MRT((l-1)*M+1:M*l,jxx));
                mui_FZF_temp = mui_FZF_temp + (H_matrix((l-1)*M+1:M*l,jx)'*sqrt(power2)*Wprecoder_FZF((l-1)*M+1:M*l,jxx));
                mui_PZF_temp = mui_PZF_temp + (H_matrix((l-1)*M+1:M*l,jx)'*sqrt(power2)*Wprecoder_PZF((l-1)*M+1:M*l,jxx));
                mui_PZF_MRT_temp = mui_PZF_MRT_temp + (H_matrix((l-1)*M+1:M*l,jx)'*sqrt(power2)*Wprecoder_PZF_MRT((l-1)*M+1:M*l,jxx));
                mui_PPZF_MRT_temp = mui_PPZF_MRT_temp + (H_matrix((l-1)*M+1:M*l,jx)'*sqrt(power2)*Wprecoder_PPZF_MRT((l-1)*M+1:M*l,jxx));
            end
            MUI_MRT(jx,expr) = MUI_MRT(jx,expr) + abs(mui_MRT_temp)^2;
            MUI_FZF(jx,expr) = MUI_FZF(jx,expr) + abs(mui_FZF_temp)^2;
            MUI_PZF(jx,expr) = MUI_PZF(jx,expr) + abs(mui_PZF_temp)^2;
            MUI_PZF_MRT(jx,expr) = MUI_PZF_MRT(jx,expr) + abs(mui_PZF_MRT_temp)^2;
            MUI_PPZF_MRT(jx,expr) = MUI_PPZF_MRT(jx,expr) + abs(mui_PPZF_MRT_temp)^2;
        end

    end

end

Desiredx_MRT = mean(real(Desired_MRT),2);
MUIx_MRT = mean(MUI_MRT,2);
Desiredx_FZF = mean(real(Desired_FZF),2);
MUIx_FZF = mean(MUI_FZF,2);
Desiredx_PZF = mean(real(Desired_PZF),2);
MUIx_PZF = mean(MUI_PZF,2);
Desiredx_PZF_MRT = mean(real(Desired_PZF_MRT),2);
MUIx_PZF_MRT = mean(MUI_PZF_MRT,2);
Desiredx_PPZF_MRT = mean(real(Desired_PPZF_MRT),2);
MUIx_PPZF_MRT = mean(MUI_PPZF_MRT,2);




%Computing for MRT rate
rate_MRT=zeros(1,K);

for kj=1:K
    sinr_MRT  = Desiredx_MRT(kj)^2/( (MUIx_MRT(kj)-Desiredx_MRT(kj)^2) + 1 );
    rate_MRT(kj)= log2(1 + sinr_MRT );
end


%Computing for FZF rate
rate_FZF=zeros(1,K);
for dj=1:K

    sinr_FZF = Desiredx_FZF(dj)^2/( (MUIx_FZF(dj)-Desiredx_FZF(dj)^2) + 1 );
    rate_FZF(dj)= log2(1 + sinr_FZF);
end

%Computing for PZF rate
rate_PZF=zeros(1,K);
for pj=1:K
    sinr_PZF = Desiredx_PZF(pj)^2/( (MUIx_PZF(pj)-Desiredx_PZF(pj)^2) + 1 );
    rate_PZF(pj)= log2(1 + sinr_PZF);
end

%Computing for PZF_MRT rate
rate_PZF_MRT=zeros(1,K);
for fj=1:K
    sinr_PZF_MRT = Desiredx_PZF_MRT(fj)^2/( (MUIx_PZF_MRT(fj)-Desiredx_PZF_MRT(fj)^2) + 1 );
    rate_PZF_MRT(fj)= log2(1 + sinr_PZF_MRT);
end

rate_PZF_total = rate_PZF+rate_PZF_MRT;

%Computing for PPZF_MRT rate
rate_PPZF_MRT=zeros(1,K);
for qj=1:K
    sinr_PPZF_MRT = Desiredx_PPZF_MRT(qj)^2/( (MUIx_PPZF_MRT(qj)-Desiredx_PPZF_MRT(qj)^2) + 1 );
    rate_PPZF_MRT(qj)= log2(1 + sinr_PPZF_MRT);
end

rate_PPZF_total = rate_PPZF_MRT + rate_PZF;










sumrate_MRT = sum(rate_MRT);
sumrate_FZF = sum(rate_FZF);
sumrate_PZF = sum(rate_PZF_total);
sumrate_PPZF = sum(rate_PPZF_total);
end