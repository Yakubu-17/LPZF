tic
clear
clc
disp(['Program Execution Starts !!!' ]);%datetime('now')

M=20:20:100; Nrx=2; %Number of APs and antenna arrays per AP 
Ntx=Nrx;MNrx=M*Nrx;
MNtx=M*Ntx; %Number of RF antenna chains;
Ku=15; Kd=Ku;%Number of UL/DL users

tauD = 5; %downlink pilot length
tauU = 5; %Uplink pilot length
tau = tauD + tauU;
[U,S,V]=svd(randn(tau,tau)); %generate orthogonal pilots
Uu = U(:,[1:tauU]); %Uplink user pilot sequences
Ud = U(:,[1+tauD:tau]); %DL pilot sequences


[Ud1,~,~]=svd(randn(tauD,tauD)); %generate orthogonal pilots


pilot_book=Ud1;%Pilot book matrix

Tc=200; Td=(Tc-tau)/Tc; %Pre-log factor

Td1 = 0.5*(1 - tauD/Tc);

noise_pow=-92; %noise power in dBm
pudB=23.01-noise_pow; pddB=30-noise_pow; 
%ppdB=23.01-noise_pow;%Normalized DL, UL, pilot powers
sigma_si_sigma_0=40; %normalized residual SI (normalized by noise)
RIdB = sigma_si_sigma_0 +  noise_pow - 30; %Convert the normalized noise to absolute dB value
si_canc = 10.^(RIdB./10); %Absolute SI power

rhohd=zeros(Kd,Ku); mW = 200; puhd = (10*log10(mW))-noise_pow; Mw = 2000; pdhd = (10*log10(Mw))-noise_pow;
ppdB = (10*log10(100))-noise_pow;





kappa=0.0; kappa1=kappa; %Resolution of ADCs: refer to paper for values. 0 refers to perfect ADC  
epsilon=(1-kappa1).*ones(1,Kd);%ADC Resoln for the DL users

%Intialize for storage

sumSEdTheoryhd = zeros(1,length(M));
sumSEsim = zeros(1,length(M));
sumSEdMonteCarlohd = zeros(1,length(M));
sumSEcloseform = zeros(1,length(M));

lsexperiment = 100; %large-scale fading experiments.
experiment = 100; %Monte-Carlo simulations
for idx=1:length(M)
    %initialize for averaging over large realizations
    rateu_theoryhd = zeros(1,lsexperiment); %Analytical UL rate per realization
    rated_theoryhd = zeros(1,lsexperiment); %Analytical DL rate per realization
    rateu_MonteCarlohd = zeros(1,lsexperiment); %Simluated UL
    rated_MonteCarlohd = zeros(1,lsexperiment); %Simluated DL
    rate_MonteCarloMRT = zeros(1,lsexperiment);
    rate_MRT = zeros(1,lsexperiment);
    alpha = (1-kappa).*ones(1,M(idx)); %Use the same resolution for all APs
    for ix = 1:lsexperiment
       
        %Assign pilot sequences randomly
         Phii=zeros(tau,Ku); %Initialize for UL pilots
         for k=1:Ku
          Point=randi([1,tauU]);
          Phii(:,k)=Uu(:,Point);
         % Phii(:,k)=Uu(:,k); %Assign unique pilot sequences if available
         end
         Phiid=zeros(tau,Kd); %Initialize for DL pilots
         for kk=1:Kd
            point=randi([1,tauD]);
            Phiid(:,kk)=Ud(:,point);
%           Phiid(:,kk)=Ud(:,kk); %Assign unique pilot sequences if available
         end


         Phiid1=zeros(tauD,Kd); %Initialize for DL pilots
         Pilot_indices = zeros(1,Kd);
        for kk=1:Kd
            point=randi([1,tauD]);
            Pilot_indices(kk)=point;
            Phiid1(:,kk)=Ud1(:,point);
         end


       
        [bu,bd,sigma,rho,~,~] = LSF(M(idx),Ku,Kd,si_canc); %return the large scale fading coefficients (UL, DL, SI and UL-to-DL) channel gains 
        sigmahd=zeros(M(idx),M(idx)); 
       
        [~, rated_MonteCarlohd(ix)] = functionSimulationP(M(idx),Nrx,Ntx,alpha,epsilon,Ku,Kd,bu,bd,sigmahd,rhohd,puhd,pdhd,ppdB,tau,Phii,Phiid,experiment); %Simulated rate per realization
        [rate_MonteCarloMRT(ix),~] = functionSimulation(M(idx),Ntx,Kd,pdhd,ppdB,bd,Phiid1,tauD,Pilot_indices,pilot_book,experiment); %Sum rate 
        [~, rated_theoryhd(ix)] = functionTheoretical(M(idx),Nrx,Ntx,alpha,epsilon,Ku,Kd,bu,bd,sigmahd,rhohd,puhd,pdhd,ppdB,tau,Phii,Phiid);
        [rate_MRT(ix),~] = functionClosedForm(M(idx),Ntx,Kd,pdhd,ppdB,bd,Phiid1,tauD); %Sum rate 
    end
    %Average rates per large scale realization and impose pre-log factor.
    %Store values for per AP number
    sumSEsim(idx) = Td1*mean(rate_MonteCarloMRT); 
    sumSEdTheoryhd(idx) = 0.5*Td*mean(rated_theoryhd);
    sumSEcloseform(idx) =Td1*mean(rate_MRT);
    sumSEdMonteCarlohd(idx) =0.5* Td*mean(rated_MonteCarlohd);
     
   %Display status of the computation
    disp([ 'Status:  ' num2str( (idx/length(M))*100   ) '% '  'UL Theoretical: '  num2str(sumSEsim(idx)) '   Monte Carlo:  ' num2str(sumSEcloseform(idx))...
        '  ||  DL Theoretical: ' num2str(sumSEdTheoryhd(idx))  '  Monte Carlo: '  num2str(sumSEdMonteCarlohd(idx))])    
end
disp('Execution Complete !!!')
figure(4)
%ax1 = subplot(2,1,1);
hold on; plot(M,sumSEsim,'*r',M,sumSEcloseform,'-k',M,sumSEdMonteCarlohd,'sr',M,sumSEdTheoryhd,'^b');
legend('They(MCS)','They(CF)','Us(MCS)','Us(CF)','interpreter','latex')
xlabel('Number of APs $$M$$','interpreter','latex')
ylabel('DL sum SE (bps/Hz)','interpreter','latex')
box on;

% ax2=subplot(2,1,2);
% hold on; plot(M,sumSEdMonteCarlohd,'sr',M,sumSEdTheoryhd,'-k');
% legend('Simulations','Closed-form','interpreter','latex')
% xlabel('Number of APs $$M$$','interpreter','latex')
% ylabel('DL sum SE (bps/Hz)','interpreter','latex')
% box on
toc