tic
clear
clc
disp(['Program Execution Starts !!!' ]);%datetime('now')

M=100; Nrx=4; %Number of APs and antenna arrays per AP 
Ntx=Nrx;MNrx=M*Nrx;
MNtx=M*Ntx; %Number of RF antenna chains;
Ku=30; Kd=Ku;%Number of UL/DL users

tauD = 5; %downlink pilot length
tauU = 5; %Uplink pilot length
tau = tauD + tauU;
[U,S,V]=svd(randn(tau,tau));
Uu = U(:,[1:tauU]);
Ud = U(:,[1+tauD:tau]);
Tc=200; Td=(Tc-tau)/Tc; %Td=1;

Td1 = 0.5*(1 - tau/Tc);

noise_pow=-92;
pudB=23.01-noise_pow; pddB=30-noise_pow; 
%ppdB=23.01-noise_pow;%Transmit powers
ppdB = (10*log10(100))-noise_pow;
sigma_si_sigma_0=100;
RIdB = sigma_si_sigma_0 +  noise_pow - 30;
si_canc = 10.^(RIdB./10);

%bit=7; %kappa1=(pi*sqrt(3)/2)*2^(-2*bit);
kappa=0.0; kappa1=kappa;  
epsilon=(1-kappa1).*ones(1,Kd);%ADC Resoln
sigmahd=zeros(M,M); rhohd=zeros(Kd,Ku); mW = 200; puhd = (10*log10(mW))-noise_pow; Mw = 2000; pdhd = (10*log10(Mw))-noise_pow;


lsexperiment = 1000;
rateu_theory = zeros(1,lsexperiment);
rated_theory = zeros(1,lsexperiment);

rateu_theoryhd = zeros(1,lsexperiment);
rated_theoryhd = zeros(1,lsexperiment);
rate_MRT = zeros(1,lsexperiment);
alpha = (1-kappa).*ones(1,M);
for ix = 1:lsexperiment
       
         Phii=zeros(tau,Ku);
         for k=1:Ku
          Point=randi([1,tauU]);
          Phii(:,k)=Uu(:,Point);
         % Phii(:,k)=Uu(:,k);
         end

         Phiid=zeros(tau,Kd);
         for kk=1:Kd
          point=randi([1,tauD]);
          Phiid(:,kk)=Ud(:,point);
%           Phiid(:,kk)=Ud(:,kk);
         end
       
        [bu,bd,sigma,rho,~,~] = LSF(M,Ku,Kd,si_canc);
        %[rateu_theory(ix), rated_theory(ix)] = functionTheoretical(M,Nrx,Ntx,alpha,epsilon,Ku,Kd,bu,bd,sigma,rho,pudB,pddB,ppdB,tau,Phii,Phiid); 
        [rateu_theoryhd(ix), rated_theoryhd(ix)] = functionTheoretical(M,Nrx,Ntx,alpha,epsilon,Ku,Kd,bu,bd,sigmahd,rhohd,puhd,pdhd,ppdB,tau,Phii,Phiid);
        [rate_MRT(ix),~] = functionClosedForm(M,Ntx,Kd,pdhd,ppdB,bd,Phiid,tauD); %Sum rate 
        ix
end
    %sumSETheory = Td*(rateu_theory + rated_theory);
    UL = 0.5*Td*(rateu_theoryhd );
    DL = 0.5*Td*( rated_theoryhd);
    DLT = Td1*( rate_MRT);
figure(1)
hold on; h = cdfplot(DLT);
hold on; g =  cdfplot(DL);
set(h,'LineStyle','-','Color','b');
set(g,'LineStyle','--','Color','r');
legend('W','P','interpreter','latex')
xlabel('Sum SE (bps/Hz)','interpreter','latex');
ylabel('Cummulative Distribution','interpreter','latex')

