
tic
clear
clc
disp('Program Execution Starts !!!' );

L=200; %Number of APs 
K=20;%Number of UEs in the network
M=16;%Number of antennas per AP
Tc=200;%Length of the coherence block
tau=15;%Number of pilots per coherence block
Kd=K;
si_canc=0;


xi=0.5; %symmetric TDD frame.
Td = xi*(1 - tau/Tc); %Pre log factor.
[U,~,~]=svd(randn(tau,tau)); %generating othogonal pilot sequences.
pilot_book=U;%Pilot book matrix


noise_power = -92; %Noise power in dBm 
pp=100; %Uplink transmit power per UE for pilot(mW)
ppdB=10*log10(pp)-noise_power; %pilot power in dBm 
p_max = 200; %maximum power(mW)
p_maxdB = 10*log10(p_max)-noise_power; %p_max in dBm.


threshold = 0.85;



experiment = 20; %Monte-Carlo simulation experiment
lsf_experiment =400; %Experiment for large scale fading
%sum_MonteCarlo = zeros(1,exp);
%for id=1:exp
all_rate_MRT = zeros(lsf_experiment,K); 
all_rate_FZF = zeros(lsf_experiment,K); 
all_rate_PZF = zeros(lsf_experiment,K);
all_rate_PPZF = zeros(lsf_experiment,K); 

    for ic=1:lsf_experiment
        Phii=zeros(tau,K);
        Pilot_indices = zeros(1,K);
        for k=1:K
            Point=randi([1,tau]);
            Pilot_indices(k)=Point;
            Phii(:,k)=U(:,Point);
        end

     [~,beta_matrix,~,~,~,~] = LSF(L,K,Kd,si_canc); %Large scale fading
     [S,W,Z,X,E_S_temp,tau_S_l] = functionUEgrouping(M,beta_matrix,Pilot_indices,tau,threshold); %UE grouping

    % [rate_MonteCarloMRT(ic), rate_MonteCarloFZF(ic) ] = functionSimulation(L,M,K,p_maxdB,ppdB,beta_matrix,Phii,tau, Pilot_indices,pilot_book,experiment); %Sum rate 
     [rate_MRT, rate_FZF,rate_PZF,rate_PPZF ] = functionClosedForm(L,M,K,p_maxdB,ppdB,beta_matrix,Phii,tau,S ,tau_S_l); %Sum rate 

     all_rate_MRT(ic,:) = rate_MRT;
     all_rate_FZF(ic,:) = rate_FZF;
     all_rate_PZF(ic,:) = rate_PZF;
     all_rate_PPZF(ic,:) = rate_PPZF;


   disp(['Setup ' num2str(ic) ' out of ' num2str(lsf_experiment)]);
    end


   %Averaging rate and applying pre log factor.

    MRT_sumrate=Td*(sum(all_rate_MRT,2));
    FZF_sumrate=Td*(sum(all_rate_FZF,2));
    PZF_sumrate=Td*(sum(all_rate_PZF,2));
    PPZF_sumrate=Td*(sum(all_rate_PPZF,2));
    MRT_peruser=Td*(all_rate_MRT(:));
    FZF_peruser=Td*(all_rate_FZF(:));
    PZF_peruser=Td*(all_rate_PZF(:));
    PPZF_peruser=Td*(all_rate_PPZF(:));
    
    
disp('Execution Complete !!!')

figure(1);
hold on;
a = cdfplot(MRT_sumrate);
b = cdfplot(FZF_sumrate);
c = cdfplot(PZF_sumrate);
d = cdfplot(PPZF_sumrate);
set(a,'LineStyle','-','Color','b','LineWidth', 1.2);
set(b,'LineStyle','-','Color','r','LineWidth', 1.2);
set(c,'LineStyle','-','Color','k','LineWidth', 1.2);
set(d,'LineStyle','-','Color',"#EDB120",'LineWidth', 1.2);
legend('MRT','FZF','PZF','PPZF','interpreter','latex')
xlabel('DL Sum SE [bit/s/Hz]','interpreter','latex');
ylabel('CDF','interpreter','latex')



figure(2)
hold on;
a = cdfplot(MRT_peruser);
b = cdfplot(FZF_peruser);
c = cdfplot(PZF_peruser);
d = cdfplot(PPZF_peruser);
set(a,'LineStyle','-','Color','b','LineWidth', 1.2);
set(b,'LineStyle','-','Color','r','LineWidth', 1.2);
set(c,'LineStyle','-','Color','k','LineWidth', 1.2);
set(d,'LineStyle','-','Color',"#EDB120",'LineWidth', 1.2);
legend('MRT','FZF','PZF','PPZF','interpreter','latex')
xlabel('DL SE [bit/s/Hz/user]','interpreter','latex');
ylabel('CDF','interpreter','latex')

smaller_axes = axes('Position', [0.6, 0.2, 0.3, 0.3]);
hold on; 
a = cdfplot(MRT_peruser);
b = cdfplot(FZF_peruser);
c = cdfplot(PZF_peruser);
d = cdfplot(PPZF_peruser);
set(a,'LineStyle','-','Color','b','LineWidth', 1.2);
set(b,'LineStyle','-','Color','r','LineWidth', 1.2);
set(c,'LineStyle','-','Color','k','LineWidth', 1.2);
set(d,'LineStyle','-','Color',"#EDB120",'LineWidth', 1.2);
xlim([1.5, 3]); % Set the desired xmin and xmax
ylim([0.01, 0.3]); % Set the desired ymin and ymax
xlabel('');
ylabel('');
title('');
toc
