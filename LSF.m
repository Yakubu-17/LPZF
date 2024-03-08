%%This funtion returns the large-scale fading coefficients and
%%line-of-sight gains

function [BETAA,BETAAd,BETAARi,BETAA_udi,rku,rkd] = LSF(M,K,Kd,si_canc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rng('default')
%Downlink
%Consider a square are of DxD m^2
%M distributed APs serves K terminals, they all randomly located in the area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inital parameters
%M=6; %number of access points
%Ntx = 2;%number of transmit antennas per AP
%Nrx = 2;%number of receive antennas per AP
%K=2; %number of terminals
%Kd=2;
D=1; %in kilometer
%tau=20;%training length
%[U,S,V]=svd(randn(tau,tau));%U includes tau orthogonal sequences 

% B=20; %Mhz
% Hb = 15; % Base station height in m
% Hm = 1.65; % Mobile height in m
% f = 1900; % Frequency in MHz
% aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
%L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;

% power_f=0.2; %downlink power: 200 mW
% noise_p = 10^((-203.975+10*log10(20*10^6)+9)/10); %noise power
% Pd = power_f/noise_p;%nomalized receive SNR
% Pp=Pd;%pilot power

d0=0.01;%km
d1=0.05;%km
BSI = -81.1846;
L = 140.7;
%%%%%Randomly locations of M APs%%%%
AP=zeros(M,2,9);
AP(:,:,1)=unifrnd(-D/2,D/2,M,2);

%Wrapped around (8 neighbor cells)
D1=zeros(M,2);
D1(:,1)=D1(:,1)+ D*ones(M,1);
AP(:,:,2)=AP(:,:,1)+D1;

D2=zeros(M,2);
D2(:,2)=D2(:,2)+ D*ones(M,1);
AP(:,:,3)=AP(:,:,1)+D2;

D3=zeros(M,2);
D3(:,1)=D3(:,1)- D*ones(M,1);
AP(:,:,4)=AP(:,:,1)+D3;

D4=zeros(M,2);
D4(:,2)=D4(:,2)- D*ones(M,1);
AP(:,:,5)=AP(:,:,1)+D4;

D5=zeros(M,2);
D5(:,1)=D5(:,1)+ D*ones(M,1);
D5(:,2)=D5(:,2)- D*ones(M,1);
AP(:,:,6)=AP(:,:,1)+D5;

D6=zeros(M,2);
D6(:,1)=D6(:,1)- D*ones(M,1);
D6(:,2)=D6(:,2)+ D*ones(M,1);
AP(:,:,7)=AP(:,:,1)+D6;

D7=zeros(M,2);
D7=D7+ D*ones(M,2);
AP(:,:,8)=AP(:,:,1)+D7;

D8=zeros(M,2);
D8=D8- D*ones(M,2);
AP(:,:,9)=AP(:,:,1)+D8;

%Randomly locations of K terminals:
Ter=zeros(K,2,9);
Ter(:,:,1)=unifrnd(-D/2,D/2,K,2);

%Wrapped around (8 neighbor cells)
D1=zeros(K,2);
D1(:,1)=D1(:,1)+ D*ones(K,1);
Ter(:,:,2)=Ter(:,:,1)+D1;

D2=zeros(K,2);
D2(:,2)=D2(:,2)+ D*ones(K,1);
Ter(:,:,3)=Ter(:,:,1)+D2;

D3=zeros(K,2);
D3(:,1)=D3(:,1)- D*ones(K,1);
Ter(:,:,4)=Ter(:,:,1)+D3;

D4=zeros(K,2);
D4(:,2)=D4(:,2)- D*ones(K,1);
Ter(:,:,5)=Ter(:,:,1)+D4;

D5=zeros(K,2);
D5(:,1)=D5(:,1)+ D*ones(K,1);
D5(:,2)=D5(:,2)- D*ones(K,1);
Ter(:,:,6)=Ter(:,:,1)+D5;

D6=zeros(K,2);
D6(:,1)=D6(:,1)- D*ones(K,1);
D6(:,2)=D6(:,2)+ D*ones(K,1);
Ter(:,:,7)=Ter(:,:,1)+D6;

D7=zeros(K,2);
D7=D7+ D*ones(K,2);
Ter(:,:,8)=Ter(:,:,1)+D7;

D8=zeros(K,2);
D8=D8- D*ones(K,2);
Ter(:,:,9)=Ter(:,:,1)+D8;

%Randomly locations of Kd DL terminals:
Terd=zeros(Kd,2,9);
Terd(:,:,1)=unifrnd(-D/2,D/2,Kd,2);

%Wrapped around (8 neighbor cells)
D1=zeros(Kd,2);
D1(:,1)=D1(:,1)+ D*ones(Kd,1);
Terd(:,:,2)=Terd(:,:,1)+D1;

D2=zeros(Kd,2);
D2(:,2)=D2(:,2)+ D*ones(Kd,1);
Terd(:,:,3)=Terd(:,:,1)+D2;

D3=zeros(Kd,2);
D3(:,1)=D3(:,1)- D*ones(Kd,1);
Terd(:,:,4)=Terd(:,:,1)+D3;

D4=zeros(Kd,2);
D4(:,2)=D4(:,2)- D*ones(Kd,1);
Terd(:,:,5)=Terd(:,:,1)+D4;

D5=zeros(Kd,2);
D5(:,1)=D5(:,1)+ D*ones(Kd,1);
D5(:,2)=D5(:,2)- D*ones(Kd,1);
Terd(:,:,6)=Terd(:,:,1)+D5;

D6=zeros(Kd,2);
D6(:,1)=D6(:,1)- D*ones(Kd,1);
D6(:,2)=D6(:,2)+ D*ones(Kd,1);
Terd(:,:,7)=Terd(:,:,1)+D6;

D7=zeros(Kd,2);
D7=D7+ D*ones(Kd,2);
Terd(:,:,8)=Terd(:,:,1)+D7;

D8=zeros(Kd,2);
D8=D8- D*ones(Kd,2);
Terd(:,:,9)=Terd(:,:,1)+D8;


sigma_shd=8; %in dB

%%%%M correlated shadowing cofficients of M APs:
Dist=zeros(M,M);%distance matrix
BETAARi=zeros(M,M); %SI matrix
Z_shdri = sigma_shd * normrnd(0,1,[M,M]);
% Z_shadow_u = lognrnd(0,db2mag(sigma_shd),[M,Kd]);
% Z_shd = db(Z_shadow_u);
for m1=1:M
   for m2=1:M
       Dist(m1,m2) = min([norm(AP(m1,:,1)-AP(m2,:,1)), norm(AP(m1,:,1)-AP(m2,:,2)),...
           norm(AP(m1,:,1)-AP(m2,:,3)),norm(AP(m1,:,1)-AP(m2,:,4)),norm(AP(m1,:,1)-AP(m2,:,5)),...
           norm(AP(m1,:,1)-AP(m2,:,6)),norm(AP(m1,:,1)-AP(m2,:,7)),norm(AP(m1,:,1)-AP(m2,:,8)),...
           norm(AP(m1,:,1)-AP(m2,:,9)) ]); %distance between AP m1 and AP m2
      % Cor(m1,m2)=exp(-log(2)*Dist(m1,m2)/D_cor);
       
    %SI Channel gain
    if Dist(m1,m2)<=d0
         betadBRi=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((Dist(m1,m2)>d0) && (Dist(m1,m2)<=d1))
         betadBRi= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(Dist(m1,m2));
    else
    betadBRi = -L - 35*log10(Dist(m1,m2)) + Z_shdri(m1,m2); %large-scale in dB
    end
    
    BETAARi(m1,m2)=(10^(betadBRi/10)).*si_canc; 
    if m1==m2
       BETAARi(m1,m2) = (10^(BSI/10)).*si_canc;  
    end
    
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create an MxK large-scale coefficients beta_mk
BETAA = zeros(M,K);
dist=zeros(M,K);
Z_shd = sigma_shd * normrnd(0,1,[M,K]);
% Z_shadow_u = lognrnd(0,db2mag(sigma_shd),[M,Kd]);
% Z_shd = db(Z_shadow_u);
rku = zeros(M,K);
for m=1:M  
    for k=1:K
    [dist(m,k),index] = min([norm(AP(m,:,1)-Ter(k,:,1)), norm(AP(m,:,2)-Ter(k,:,1)),...
        norm(AP(m,:,3)-Ter(k,:,1)),norm(AP(m,:,4)-Ter(k,:,1)),norm(AP(m,:,5)-Ter(k,:,1)),...
        norm(AP(m,:,6)-Ter(k,:,1)),norm(AP(m,:,7)-Ter(k,:,1)),norm(AP(m,:,8)-Ter(k,:,1)),...
        norm(AP(m,:,9)-Ter(k,:,1)) ]); %distance between Terminal k and AP m
        if dist(m,k)<d0
            betadB=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
        elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
            betadB= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
        else
         betadB = -L - 35*log10(dist(m,k)) + Z_shd(m,k); %large-scale in dB
        end
    
    BETAA(m,k)=10^(betadB/10);
    rku(m,k) = 10^(1.3-0.003*dist(m,k)*1000);
    end

end

%Create an MxK large-scale coefficients beta_mk
BETAAd=zeros(M,Kd);
distd=zeros(M,Kd);
Z_shdd = sigma_shd * normrnd(0,1,[M,Kd]); 

rkd = zeros(M,Kd);
for md=1:M  
    for kd=1:Kd
    [distd(md,kd),index1] = min([norm(AP(md,:,1)-Terd(kd,:,1)), norm(AP(md,:,2)-Terd(kd,:,1)),...
        norm(AP(md,:,3)-Terd(kd,:,1)),norm(AP(md,:,4)-Terd(kd,:,1)),norm(AP(md,:,5)-Terd(kd,:,1)),...
        norm(AP(md,:,6)-Terd(kd,:,1)),norm(AP(md,:,7)-Terd(kd,:,1)),norm(AP(md,:,8)-Terd(kd,:,1)),...
        norm(AP(md,:,9)-Terd(kd,:,1)) ]); %distance between Terminal k and AP m
    if distd(md,kd)<d0
         betadBd=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((distd(md,kd)>=d0) && (distd(md,kd)<=d1))
         betadBd= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(distd(md,kd));
    else
    betadBd = -L - 35*log10(distd(md,kd)) + Z_shdd(md,kd); %large-scale in dB
    end
    
    BETAAd(md,kd)=10^(betadBd/10); 
    end
    rkd(md,kd)=10^(1.3-0.003*distd(md,kd)*1000);

end

%Create a Kd-by-K UDI channel gains
BETAA_udi=zeros(Kd,K);
Distudi = zeros(Kd,K);
Z_shd_udi = sigma_shd * normrnd(0,1,[Kd,K]); 
for kkd=1:Kd
    for kk=1:K
        Distudi(kkd,kk)=min([norm(Terd(kkd,:,1)-Ter(kk,:,1)), norm(Terd(kkd,:,1)-Ter(kk,:,2)),...
            norm(Terd(kkd,:,1)-Ter(kk,:,3)),norm(Terd(kkd,:,1)-Ter(kk,:,4)),...
            norm(Terd(kkd,:,1)-Ter(kk,:,5)),norm(Terd(kkd,:,1)-Ter(kk,:,6)),...
            norm(Terd(kkd,:,1)-Ter(kk,:,7)),norm(Terd(kkd,:,1)-Ter(kk,:,8)),...
            norm(Terd(kkd,:,1)-Ter(kk,:,9)) ]); %distance between Terminal k1 and Terminal k2
       if Distudi(kkd,kk)<d0
           betadBudi=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
       elseif ((Distudi(kkd,kk)>=d0) && (Distudi(kkd,kk)<=d1))
           betadBudi= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(Distudi(kkd,kk));
       else
           betadBudi = -L - 35*log10(Distudi(kkd,kk)) + Z_shd_udi(kkd,kk); %large-scale in dB
       end
    
    BETAA_udi(kkd,kk)=(10^(betadBudi/10)).*si_canc; 
    end
    
end


end
