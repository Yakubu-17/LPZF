clear;
clc;


L=20; %Number of APs 
K=10;%Number of UEs in the network
M=5;%Number of antennas per AP
Tc=200;%Length of the coherence block
tau=4;%Number of pilots per coherence block
nbrOfRealizations = 50;
p=1e14;
LM = L*M;
si_canc=0;
Kd=K;

[U,~,~]=svd(randn(tau,tau)); %generating othogonal pilot sequences.


pilot_book =U;
Phii=zeros(tau,K);
pilotIndex = zeros(1,K);
        for k=1:K
            Point=randi([1,tau]);
            pilotIndex(k)=Point;
            Phii(:,k)=U(:,Point);
        end




% function [Hprecoder_MRT1, Hprecoder_FZF1, Hhat] = test_CE(p,tau,Phii,H_matrix,B,pilot_book,pilotIndex)
% 
% [rows, columns]=size(H_matrix);
% LM=rows;
% K=columns;







[~,Beta,~,~,~,~] = LSF(L,K,Kd,si_canc); %Large scale fading
B = repelem(Beta,M,1);

%Fi_matrix = sqrt(tau)*eye(tau);

H = ((randn(LM,K)+ 1i*randn(LM,K))/sqrt(2)).*(B).^(1/2);
N = ((randn(LM,tau)+ 1i*randn(LM,tau))/sqrt(2));
Yp1 = sqrt(p)*tau*H*Phii' + sqrt(tau)*N;
Yp = sqrt(p)*H*Phii' + N;
yp = Yp1*Phii;
H_bar = Yp1*pilot_book;
e = eye(tau);
% Initialize the precoding matrix H_precoder of size LM by K
Hprecoder_MRT = zeros(LM,K);
Hprecoder_FZF = zeros(LM,K);
% Loop through each user (k) to form the precoding matrix
for k = 1:K
   
    %H_bar=Yp1*pilot_book_matrix(:,:,k);
    Hprecoder_MRT(:,k) = H_bar*e(:,pilotIndex(k));
    Hprecoder_FZF(:,k) = H_bar*(H_bar'*H_bar)^(-1)*e(:,pilotIndex(k));
end

yl = zeros(LM,tau);
for l=1:LM
    yl_temp =0;
    for k =1:K
       yl_temp = yl_temp + H(l,k)*Phii(:,k)' ;
    end
  yl(l,:) = sqrt(p)*tau*yl_temp +  sqrt(tau)*N(l,:) ;
end

H_bar1 = zeros(LM,tau); 
for l=1:LM
    H_bar1(l,:) =yl(l,:)*pilot_book;
end

Hprecoder_MRT1 = zeros(LM,K);
Hprecoder_FZF1 = zeros(LM,K);

for l =1:LM
    for k = 1:K
        Hprecoder_MRT1(l,k) = H_bar1(l,:)*e(:,pilotIndex(k));
        Hprecoder_FZF1(l,k) = H_bar1(l,:)*(H_bar1(l,:)'*H_bar1(l,:))^(-1)*e(:,pilotIndex(k));
    end
end   



lambda = zeros(LM,K);

for m=1:LM
    for k=1:K
        den = norm( (B(m,:).^(1/2)).*(Phii(:,k)'*Phii))^2;
        lambda(m,k)=(tau*p*den +1); 
    end
end




Hhat = (sqrt(p).*B./(lambda)).*yp;
c_matrix = (sqrt(p).*B./(lambda));
%New Hhat
Hhat1 = zeros(LM,K);
for l=1:LM
    for k = 1:K
        Hhat1(l,k) =  c_matrix(l,k)*H_bar1(l,:)*e(:,pilotIndex(k));
    end
end


[S,W,Z,M,E_S_temp,tau_S_l2] = functionUEgrouping(LM,K,B,pilotIndex,tau);

Hprecoder_PZF = zeros(LM,K);
Hprecoder_PZF_MRT = zeros(LM,K);

  for k = 1:K

    %Strong UE
    [~,Z_k_ind] = find(Z(k,:) == 1);
    Z_k_ind_length = length(Z_k_ind);
    
    for l_Z_k = 1:Z_k_ind_length
            E_S_l_ind = find(E_S_temp(Z_k_ind(l_Z_k),:)~=0);
            tau_S_l = length(E_S_l_ind);
            E_S_l = zeros(tau,tau_S_l);
            eyetau_S_l = eye(tau_S_l );

            for r = 1:tau_S_l
                E_S_l(:,r) = e(:,E_S_l_ind(r));
            end

       epso = eyetau_S_l(:,find(E_S_l_ind == pilotIndex(k)));
       Hprecoder_PZF(Z_k_ind(l_Z_k),k) = H_bar1(Z_k_ind(l_Z_k),:)*E_S_l*pinv(E_S_l'*H_bar1(Z_k_ind(l_Z_k),:)'*H_bar1(Z_k_ind(l_Z_k),:)*E_S_l)*epso;
    end


    %Weak UE
    [~,M_k_ind] = find(M(k,:) == 1);
    M_k_ind_length = length(M_k_ind);

     for  l_M_k = 1:M_k_ind_length

           Hprecoder_PZF_MRT(M_k_ind(l_M_k),k) = H_bar1(M_k_ind(l_M_k),:)*e(:,pilotIndex(k));
     end
   
 end

%end
%Hprecoder_MRT2 = zeros(LM,K);


% for l =1:LM
%     for k = 1:K
%         %Hprecoder_MRT1(l,k) = H_bar1(l,:)*e(:,pilotIndex(k));
%         
%     end
% end   

size(Hhat1);
%end

