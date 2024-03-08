
tic
clear
clc
disp(['Program Execution Starts !!!' ]);

L=20; %Number of APs 
K=5;%Number of UEs in the network
M=5:1:10;%Number of antennas per AP
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



%Storing rate for each instance.
sum_MonteCarloMRT = zeros(1,length(M));
sum_MonteCarloFZF = zeros(1,length(M));
sum_closedForm_MRT = zeros(1,length(M));
sum_closedForm_FZF = zeros(1,length(M));

experiment = 50; %Monte-Carlo simulation experiment
lsf_experiment = 30; %Experiment for large scale fading

for id=1:length(M)
rate_MonteCarloMRT = zeros(1,lsf_experiment); %Simluated 
rate_MonteCarloFZF = zeros(1,lsf_experiment); %Simluated 
closedForm_MRT = zeros(1,lsf_experiment);
closedForm_FZF = zeros(1,lsf_experiment);

    for ic=1:lsf_experiment
        Phii=zeros(tau,K);
        Pilot_indices = zeros(1,K);
        for k=1:K
            Point=randi([1,tau]);
            Pilot_indices(k)=Point;
            Phii(:,k)=U(:,Point);
        end

     [~,beta_matrix,~,~,~,~] = LSF(L,K,Kd,si_canc); %Large scale fading
     %
     %beta_matrix = ones(L(id),K);
     [rate_MonteCarloMRT(ic),rate_MonteCarloFZF(ic)] = functionSimulation(L,M(id),K,p_maxdB,ppdB,beta_matrix,Phii,tau, Pilot_indices,pilot_book,experiment); %Sum rate 
     [closedForm_MRT(ic),closedForm_FZF(ic)] = functionClosedForm(L,M(id),K,p_maxdB,ppdB,beta_matrix,Phii,tau); %Sum rate 
     %[rate_MonteCarloMRT(ic), rate_MonteCarloFZF(ic)] = functionSimulationTest(L(id),M,K,p_maxdB,ppdB,beta_matrix,Phii,tau, Pilot_indices,pilot_book,experiment); %Sum rate 
    end

    %Averaging rate and applying pre log factor.
    sum_MonteCarloMRT(id)=Td*mean(rate_MonteCarloMRT);
    sum_MonteCarloFZF(id)=Td*mean(rate_MonteCarloFZF);
    sum_closedForm_MRT(id)=Td*mean(closedForm_MRT);
    sum_closedForm_FZF(id)=Td*mean(closedForm_FZF);

disp(['Setup ' num2str(id) ' out of ' num2str(length(M))]);
end
disp('Execution Complete !!!')
gap = sum_closedForm_FZF-sum_MonteCarloFZF ;
figure(1);
hold on;plot(M,gap,'-r','LineWidth', 1.2);
%legend('Simulation','closed-form','interpreter','latex')
xlabel('Number of antennas per AP $$M$$','interpreter','latex')
ylabel('gap sum SE (bps/Hz)','interpreter','latex')
title('Performance gap')
box on
grid on

figure(2);
hold on;plot(M,sum_MonteCarloFZF,'sr',M,sum_closedForm_FZF,'-k','LineWidth', 1.2);
legend('Simulation','closed-form','interpreter','latex')
xlabel('Number of antennas per AP $$M$$','interpreter','latex')
ylabel('DL sum SE (bps/Hz)','interpreter','latex')
title('FZF')
box on
grid on
toc
