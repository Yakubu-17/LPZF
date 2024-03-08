
tic
clear;
clc;
disp('Program Execution Starts !!!' );

L=100; %Number of APs 
K=10;%Number of UEs in the network
M=6;%Number of antennas per AP
Tc=200;%Length of the coherence block
tau=5;%Number of pilots per coherence block
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






experiment = 20; %Monte-Carlo simulation experiment
lsf_experiment =100; %Experiment for large scale fading
all_rate_MRT_sim = zeros(lsf_experiment,K); %Simluated MRT
all_rate_FZF_sim = zeros(lsf_experiment,K); %Simluated FZF
all_rate_MRT = zeros(lsf_experiment,K); %Simluated MRT
all_rate_FZF = zeros(lsf_experiment,K); %Simluated FZF

    for ic=1:lsf_experiment
        Phii=zeros(tau,K);
        Pilot_indices = zeros(1,K);
        for k=1:K
            Point=randi([1,tau]);
            Pilot_indices(k)=Point;
            Phii(:,k)=U(:,Point);
        end

     [~,beta_matrix,~,~,~,~] = LSF(L,K,Kd,si_canc); %Large scale fading
     [rate_MonteCarloMRT, rate_MonteCarloFZF ] = functionSimulation(L,M,K,p_maxdB,ppdB,beta_matrix,Phii,tau, Pilot_indices,pilot_book,experiment); %Sum rate 
     [rate_MRT, rate_FZF ] = functionClosedForm(L,M,K,p_maxdB,ppdB,beta_matrix,Phii,tau); %Sum rate 
      

     all_rate_MRT_sim(ic,:) = rate_MonteCarloMRT;
     all_rate_FZF_sim(ic,:) = rate_MonteCarloFZF;
     all_rate_MRT(ic,:) = rate_MRT;
     all_rate_FZF(ic,:) = rate_FZF;


   disp(['Setup ' num2str(ic) ' out of ' num2str(lsf_experiment)]);
    end
   


    %Averaging rate and applying pre log factor.
    MRT_peruser_sim=Td*(all_rate_MRT_sim(:));
    FZF_peruser_sim=Td*(all_rate_FZF_sim(:));
    MRT_sumrate_sim=Td*(sum(all_rate_MRT_sim,2));
    FZF_sumrate_sim=Td*(sum(all_rate_FZF_sim,2));

    MRT_peruser=Td*(all_rate_MRT(:));
    FZF_peruser=Td*(all_rate_FZF(:));
    MRT_sumrate=Td*(sum(all_rate_MRT,2));
    FZF_sumrate=Td*(sum(all_rate_FZF,2));

disp('Execution Complete !!!')


xd = 1/lsf_experiment:1/lsf_experiment:1;
xp = 1/(lsf_experiment*K):1/(lsf_experiment*K):1;
ps = [1 3 7 10 13 16 20 23 27 30 34 37 40 43 47 50 53 57 61 65 68 72 76 80 83 87 90 93 97 100];


MRT_sumrate_sim_sorted = sort(MRT_sumrate_sim);
FZF_sumrate_sim_sorted = sort(FZF_sumrate_sim);
MRT_peruser_sim_sorted = sort(MRT_peruser_sim);
FZF_peruser_sim_sorted = sort(FZF_peruser_sim);


position1 = round(ps / 100 * numel(MRT_sumrate_sim_sorted));
position2 = round(ps / 100 * numel(FZF_sumrate_sim_sorted));
position3 = round(ps / 100 * numel(MRT_peruser_sim_sorted));
position4 = round(ps / 100 * numel(FZF_peruser_sim_sorted));





MRT_sumrate_sim1 = MRT_sumrate_sim_sorted(position1);
FZF_sumrate_sim1 = FZF_sumrate_sim_sorted(position2);
MRT_peruser_sim1 = MRT_peruser_sim_sorted(position3);
FZF_peruser_sim1 = FZF_peruser_sim_sorted(position4);  





%Interpolate the CDF values at the specified percentiles
interpolated_sumrate_MRT = interp1(MRT_sumrate_sim_sorted, xd, MRT_sumrate_sim1, 'linear', 'extrap');
interpolated_sumrate_FZF = interp1(FZF_sumrate_sim_sorted, xd, FZF_sumrate_sim1, 'linear', 'extrap');
interpolated_peruser_MRT = interp1(MRT_peruser_sim_sorted, xp, MRT_peruser_sim1, 'linear', 'extrap');
interpolated_peruser_FZF = interp1(FZF_peruser_sim_sorted, xp, FZF_peruser_sim1, 'linear', 'extrap');




figure(2);
hold on;
plot(sort(MRT_sumrate), xd, '-r', sort(FZF_sumrate), xd, '-b','LineWidth', 1.2);
plot(MRT_sumrate_sim1, interpolated_sumrate_MRT, '*k','LineWidth', 0.6);
plot(FZF_sumrate_sim1, interpolated_sumrate_FZF, '*k','LineWidth', 0.6);
legend('MRT','FZF','MCS','','interpreter','latex')
xlabel('DL Sum SE [bit/s/Hz]','interpreter','latex');
ylabel('CDF','interpreter','latex')
grid on
box on 
title('');

figure(3)
hold on; 
plot(sort(MRT_peruser), xp, '-r', sort(FZF_peruser), xp, '-b','LineWidth', 1.2);
plot(MRT_peruser_sim1, interpolated_peruser_MRT, '*k','LineWidth', 0.6);
plot(FZF_peruser_sim1, interpolated_peruser_FZF, '*k','LineWidth', 0.6);
legend('MRT','FZF','MCS','','interpreter','latex')
xlabel('DL SE [bit/s/Hz/user]','interpreter','latex');
ylabel('CDF','interpreter','latex')
grid on
box on 
title('');



toc
