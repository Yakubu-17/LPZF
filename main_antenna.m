
tic
clear
clc
disp(['Program Execution Starts !!!' ]);

L=50; %Number of APs 
K=10;%Number of UEs in the network
M=2:4:14;%Number of antennas per AP
Tc=200;%Length of the coherence block
tau=7;%Number of pilots per coherence block
Kd=K;
si_canc=0;


xi=0.5; %symmetric TDD frame.
Td = xi*(1 - tau/Tc); %Pre log factor.
[U,~,~]=svd(randn(tau,tau)); %generating othogonal pilot sequences.
pilot_book=sqrt(tau)*U;%Pilot book matrix


noise_power = -92; %Noise power in dBm 
pp=100; %Uplink transmit power per UE for pilot(mW)
ppdB=(10*log10(pp))-noise_power; %pilot power in dBm 
p_max = 200; %maximum power(mW)
p_maxdB = (10*log10(p_max))-noise_power; %p_max in dBm.

threshold = 0.95; %Threshold for grouping APs

%Storing rate for each instance.
sum_MonteCarloMRT = zeros(K,length(M));
sum_MonteCarloFZF = zeros(K,length(M));
sum_MonteCarloPZF = zeros(K,length(M));
sum_MonteCarloPPZF = zeros(K,length(M));
sum_closedForm_MRT = zeros(K,length(M));
sum_closedForm_FZF = zeros(K,length(M));
sum_closedForm_PZF = zeros(K,length(M));
sum_closedForm_PPZF = zeros(K,length(M));



experiment = 100; %Monte-Carlo simulation experiment
lsf_experiment = 50; %Experiment for large scale fading

for id=1:length(M)
rate_MonteCarloMRT = zeros(lsf_experiment,K); %Simluated 
rate_MonteCarloFZF = zeros(lsf_experiment,K); %Simluated
rate_MonteCarloPZF = zeros(lsf_experiment,K); %Simluated
rate_MonteCarloPPZF = zeros(lsf_experiment,K); %Simluated
closedForm_MRT = zeros(lsf_experiment,K);
closedForm_FZF = zeros(lsf_experiment,K);
closedForm_PZF = zeros(lsf_experiment,K);
closedForm_PPZF = zeros(lsf_experiment,K);

    for ic=1:lsf_experiment
        Phii=zeros(tau,K);
        Pilot_indices = zeros(1,K);
        for k=1:K
            Point=randi([1,tau]);
            Pilot_indices(k)=Point;
            Phii(:,k)=U(:,Point);
        end

     [~,beta_matrix,~,~,~,~] = LSF(L,K,Kd,si_canc); %Large scale fading
     
     
     [S,W,Z,X,E_S_temp,tau_S_l] = functionUEgrouping(M(id),beta_matrix,Pilot_indices,tau,threshold); %UE grouping

     [rate_sim_MRT,rate_sim_FZF,rate_sim_PZF,rate_sim_PPZF] = functionSimulation(L,M(id),K,p_maxdB,ppdB,beta_matrix,Phii,tau,Pilot_indices,pilot_book,S,W,Z,X,E_S_temp,tau_S_l,experiment); %Sum rate 
     rate_MonteCarloMRT(ic,:) = rate_sim_MRT;
     rate_MonteCarloFZF(ic,:) = rate_sim_FZF;
     rate_MonteCarloPZF(ic,:) = rate_sim_PZF;
     rate_MonteCarloPPZF(ic,:) = rate_sim_PPZF;

     %[closedForm_MRT(ic),closedForm_FZF(ic),closedForm_PZF(ic),closedForm_PPZF(ic)] = functionClosedForm(L,M(id),K,p_maxdB,ppdB,beta_matrix,Phii,tau,S,tau_S_l); %Sum rate 
     %[rate_MonteCarloMRT(ic), rate_MonteCarloFZF(ic)] = functionSimulationTest(L(id),M,K,p_maxdB,ppdB,beta_matrix,Phii,tau, Pilot_indices,pilot_book,experiment); %Sum rate 
    end

    %Averaging rate and applying pre log factor.
    sum_MonteCarloMRT(:,id)=Td*mean(rate_MonteCarloMRT);
    sum_MonteCarloFZF(:,id)=Td*mean(rate_MonteCarloFZF);
    sum_MonteCarloPZF(:,id)=Td*mean(rate_MonteCarloPZF);
    sum_MonteCarloPPZF(:,id)=Td*mean(rate_MonteCarloPPZF);
    sum_closedForm_MRT(:,id)=Td*mean(closedForm_MRT);
    sum_closedForm_FZF(:,id)=Td*mean(closedForm_FZF);
    sum_closedForm_PZF(:,id)=Td*mean(closedForm_PZF);
    sum_closedForm_PPZF(:,id)=Td*mean(closedForm_PPZF);

disp(['Setup ' num2str(id) ' out of ' num2str(length(M))]);
end
disp('Execution Complete !!!')

figure(1);
hold on;
plot(M,mean(sum_MonteCarloMRT),'-ob','LineWidth', 1.2);
plot(M,mean(sum_MonteCarloFZF),'-or','LineWidth', 1.2);
plot(M,mean(sum_MonteCarloPZF),'-ok','LineWidth', 1.2);
plot(M,mean(sum_MonteCarloPPZF),'-o','Color',"#EDB120",'LineWidth', 1.2);
legend('MRT','FZF','PZF','PPZF','interpreter','latex')
xlabel('Number of antennass $$M$$','interpreter','latex')
ylabel('DL sum SE (bps/Hz)','interpreter','latex')

box on
grid on


toc
