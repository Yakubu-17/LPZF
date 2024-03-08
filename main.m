
tic
clear
clc
disp(['Program Execution Starts !!!' ]);

L=10; %Number of APs 
K=7;%Number of UEs in the network
M=5;%Number of antennas per AP
Tc=200;%Length of the coherence block
tau=4;%Number of pilots per coherence block
Kd=K;
si_canc=0;


xi=0.5; %symmetric TDD frame.
Td = xi*(1 - tau/Tc); %Pre log factor.
[U,~,~]=svd(randn(tau,tau)); %generating othogonal pilot sequences.
pilot_book=U;%Pilot book matrix


noise_power = -92; %Noise power in dBm 
pp=100; %Uplink transmit power per UE for pilot(mW)
ppdB=(10*log10(pp))-noise_power; %pilot power in dBm 
p_max = 200; %maximum power(mW)
p_maxdB = (10*log10(p_max))-noise_power; %p_max in dBm.

threshold = 1; %Threshold for grouping APs

%Storing rate for each instance.
sum_MonteCarloMRT = zeros(1,length(L));
sum_MonteCarloFZF = zeros(1,length(L));
sum_MonteCarloPZF = zeros(1,length(L));
sum_MonteCarloPPZF = zeros(1,length(L));
sum_closedForm_MRT = zeros(1,length(L));
sum_closedForm_FZF = zeros(1,length(L));
sum_closedForm_PZF = zeros(1,length(L));
sum_closedForm_PPZF = zeros(1,length(L));



experiment = 70; %Monte-Carlo simulation experiment
lsf_experiment = 50; %Experiment for large scale fading

for id=1:length(L)
rate_MonteCarloMRT = zeros(1,lsf_experiment); %Simluated 
rate_MonteCarloFZF = zeros(1,lsf_experiment); %Simluated
rate_MonteCarloPZF = zeros(1,lsf_experiment); %Simluated
rate_MonteCarloPPZF = zeros(1,lsf_experiment); %Simluated
closedForm_MRT = zeros(1,lsf_experiment);
closedForm_FZF = zeros(1,lsf_experiment);
closedForm_PZF = zeros(1,lsf_experiment);
closedForm_PPZF = zeros(1,lsf_experiment);

    for ic=1:lsf_experiment
        Phii=zeros(tau,K);
        Pilot_indices = zeros(1,K);
        for k=1:K
            Point=randi([1,tau]);
            Pilot_indices(k)=Point;
            Phii(:,k)=U(:,Point);
        end

     [~,beta_matrix,~,~,~,~] = LSF(L(id),K,Kd,si_canc); %Large scale fading
     
     
     [S,W,Z,X,E_S_temp,tau_S_l,delta] = functionUEgrouping(M,beta_matrix,Pilot_indices,tau,threshold); %UE grouping

     [rate_MonteCarloMRT(ic),rate_MonteCarloFZF(ic),rate_MonteCarloPZF(ic),rate_MonteCarloPPZF(ic)] = functionSimulation(L(id),M,K,p_maxdB,ppdB,beta_matrix,Phii,tau,Pilot_indices,pilot_book,S,W,Z,X,E_S_temp,tau_S_l,experiment); %Sum rate 
     [closedForm_MRT(ic),closedForm_FZF(ic),closedForm_PZF(ic),closedForm_PPZF(ic)] = functionClosedForm(L(id),M,K,p_maxdB,ppdB,beta_matrix,Phii,tau,delta,tau_S_l); %Sum rate 
     %[rate_MonteCarloMRT(ic), rate_MonteCarloFZF(ic)] = functionSimulationTest(L(id),M,K,p_maxdB,ppdB,beta_matrix,Phii,tau, Pilot_indices,pilot_book,experiment); %Sum rate 
    end

    %Averaging rate and applying pre log factor.
    sum_MonteCarloMRT(id)=Td*mean(rate_MonteCarloMRT);
    sum_MonteCarloFZF(id)=Td*mean(rate_MonteCarloFZF);
    sum_MonteCarloPZF(id)=Td*mean(rate_MonteCarloPZF);
    sum_MonteCarloPPZF(id)=Td*mean(rate_MonteCarloPPZF);
    sum_closedForm_MRT(id)=Td*mean(closedForm_MRT);
    sum_closedForm_FZF(id)=Td*mean(closedForm_FZF);
    sum_closedForm_PZF(id)=Td*mean(closedForm_PZF);
    sum_closedForm_PPZF(id)=Td*mean(closedForm_PPZF);

disp(['Setup ' num2str(id) ' out of ' num2str(length(L))]);
end
disp('Execution Complete !!!')

% figure(14);
% hold on;plot(L,sum_MonteCarloMRT,'*r',L,sum_closedForm_MRT,'-r','LineWidth', 1.2);
% legend('Simulation','closed-form','interpreter','latex')
% xlabel('Number of APs $$L$$','interpreter','latex')
% ylabel('DL sum SE (bps/Hz)','interpreter','latex')
% title('MRT')
% box on
% grid on
% 
% figure(15);
% hold on;plot(L,sum_MonteCarloFZF,'sb',L,sum_closedForm_FZF,'-b','LineWidth', 1.2);
% legend('Simulation','closed-form','interpreter','latex')
% xlabel('Number of APs $$L$$','interpreter','latex')
% ylabel('DL sum SE (bps/Hz)','interpreter','latex')
% title('FZF')
% box on
% grid on
% 
% figure(16);
% hold on;plot(L,sum_MonteCarloPZF,'-b',L,sum_closedForm_PZF,'-k','LineWidth', 1.2);
% legend('Simulation','closed-form','interpreter','latex')
% xlabel('Number of APs $$L$$','interpreter','latex')
% ylabel('DL sum SE (bps/Hz)','interpreter','latex')
% title('PZF')
% box on
% grid on
% 
% figure(17);
% hold on;plot(L,sum_MonteCarloPPZF,'-b',L,sum_closedForm_PPZF,'-k','LineWidth', 1.2);
% legend('Simulation','closed-form','interpreter','latex')
% xlabel('Number of APs $$L$$','interpreter','latex')
% ylabel('DL sum SE (bps/Hz)','interpreter','latex')
% title('PPZF')
% box on
% grid on
toc
