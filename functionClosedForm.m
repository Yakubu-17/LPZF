%Closed form compuation
function[sumrate_MRT, sumrate_FZF, sumrate_PZF,sumrate_PPZF] = functionClosedForm(L,M,K,p_maxdB,ppdB,beta_matrix,Phii,tau,delta,tau_S_l)
%Converting p_maxdB  and pilot power to watt
p_max=10^(p_maxdB/10); pp=10^(ppdB/10);


%Solve for gamma,theta and c matrix.
nu=zeros(L,K);gamma_matrix=zeros(L,K);
%theta_matrix=zeros(L,K); c_matrix=zeros(L,K);
for l=1:L  
    for k=1:K
        d=0;
        for iidx=1:K
            d=d+ beta_matrix(l,iidx)*norm(Phii(:,k)'*Phii(:,iidx));
        end
        nu(l,k)=1 + tau*pp*d;
        gamma_matrix(l,k)=pp*tau*beta_matrix(l,k)^2*(1/nu(l,k));
%         c_matrix(l,k) = sqrt(pp)*beta_matrix(l,k)*(1/nu(l,k));
%         theta_matrix(l,k) = gamma_matrix(l,k)/(c_matrix(l,k)^2);
    end    
end


%p_control = p_max*ones(L,K);
p_control=PowerControl(gamma_matrix,p_max);

%Storing terms 
desired_signal = zeros(L,K);
desired_signal_pzf = zeros(L,K);
desired_signal_ppzf = zeros(L,K);
MUI=zeros(L,K);
PC=zeros(L,K,K);
PC_pzf=zeros(L,K,K);
PC_ppzf=zeros(L,K,K);
MUI_FZF=zeros(L,K);
MUI_PZF=zeros(L,K);
MUI_PPZF=zeros(L,K);

for l=1:L
    for k=1:K

        desired_signal(l,k)=sqrt(p_control(l,k)*gamma_matrix(l,k));
        desired_signal_pzf(l,k)=sqrt((M-delta(l,k)*tau_S_l)*p_control(l,k)*gamma_matrix(l,k));
        desired_signal_ppzf(l,k)=sqrt((M-tau_S_l)*p_control(l,k)*gamma_matrix(l,k));


        for t=1:K
            if Phii(:,t)==Phii(:,k)
                PC(l,k,t) = sqrt(p_control(l,t)*gamma_matrix(l,k));
                PC_pzf(l,k,t) = sqrt((M-delta(l,t)*tau_S_l)*p_control(l,t)*gamma_matrix(l,k));
                PC_ppzf(l,k,t) = sqrt((M-tau_S_l)*p_control(l,t)*gamma_matrix(l,k));
            end 
           MUI(l,k)= MUI(l,k) + p_control(l,t)*beta_matrix(l,k);
           MUI_FZF(l,k) = MUI_FZF(l,k) + p_control(l,t)*(beta_matrix(l,k)-gamma_matrix(l,k));
           MUI_PZF(l,k) = MUI_PZF(l,k) + p_control(l,t)*(beta_matrix(l,k)- (delta(l,t)*delta(l,k)*gamma_matrix(l,k)));
           MUI_PPZF(l,k) = MUI_PPZF(l,k) + p_control(l,t)*(beta_matrix(l,k)- (delta(l,k)*gamma_matrix(l,k)));
        end
       
    end
end




%Rate computation
rate=zeros(1,K); %initialize for the per UE SE
pck=sum(PC,1).^2; PiCont=sum(pck,3);

%desired_signal = ones(L,K);
for kj=1:K
    num = M*sum(desired_signal(:,kj))^2; 
    sinr = num/(M*(PiCont(kj)-sum(desired_signal(:,kj))^2) + sum(MUI(:,kj)) + 1);
    rate(kj) = log2(1 + sinr);
    

end

%Rate computation for FZF
rate_fzf=zeros(1,K);
for kj=1:K
    num = (M-tau)*sum(desired_signal(:,kj))^2; 
    sinr_fzf = num/(((M-tau)*PiCont(kj)-num) + sum(MUI_FZF(:,kj)) + 1);
    rate_fzf(kj) = log2(1 + sinr_fzf);
end


%Rate computation for PZF
rate_pzf=zeros(1,K);
pck_pzf=sum(PC_pzf,1).^2; PiCont_pzf=sum(pck_pzf,3);

for kx=1:K
    num = sum(desired_signal_pzf(:,kx))^2;
    sinr_pzf = num/( PiCont_pzf(kx) - sum(desired_signal_pzf(:,kx))^2 + sum(MUI_PZF(:,kx)) + 1 );
    rate_pzf(kx) = log2(1 + sinr_pzf);
end


%Rate computation for PPZF
rate_ppzf=zeros(1,K);
pck_ppzf=sum(PC_ppzf,1).^2; PiCont_ppzf=sum(pck_ppzf,3);

for kxx=1:K
    num = sum(desired_signal_ppzf(:,kxx))^2;
    sinr_ppzf = num/( PiCont_ppzf(kxx) - sum(desired_signal_ppzf(:,kxx))^2 + sum(MUI_PPZF(:,kxx)) + 1 );
    rate_ppzf(kxx) = log2(1 + sinr_ppzf);
end



sumrate_PPZF = sum(rate_ppzf);
sumrate_PZF = sum(rate_pzf);
sumrate_FZF=sum(rate_fzf);
sumrate_MRT = sum(rate);
end