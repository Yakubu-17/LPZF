clear
clc

L=20; %Number of APs 
K=10;%Number of UEs in the network
N=8;%Number of antennas per AP
Tc=200;%Length of the coherence block
tau_p=7;%Number of pilots per coherence block
nbrOfRealizations = 50;
p=100000;
[U,~,~]=svd(randn(tau_p,tau_p)); %generating othogonal pilot sequences.


Phii=zeros(tau_p,K);
pilotIndex = zeros(1,K);
        for k=1:K
            Point=randi([1,tau_p]);
            pilotIndex(k)=Point;
            Phii(:,k)=U(:,Point);
        end



%Generate uncorrelated Rayleigh fading channel realizations
H = (randn(L*N,nbrOfRealizations,K)+1i*randn(L*N,nbrOfRealizations,K))*sqrt(0.5);
%H = (randn(L*N,K,nbrOfRealizations)+1i*randn(L*N,K,nbrOfRealizations));

%Fi_matrix
Fi_matrix = sqrt(tau_p)*eye(tau_p);


%% Perform channel estimation
 

%Store identity matrix of size N x N
eyeN = eye(N);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,L,tau_p) + 1i*randn(N,nbrOfRealizations,L,tau_p));

%Prepare to store results
Hhat = zeros(L*N,nbrOfRealizations,K);

H_bar = zeros(L*N,nbrOfRealizations,tau_p);

for l=1:L

    yq  = zeros(N,tau_p,nbrOfRealizations);
    
    for n = 1:nbrOfRealizations
        
        for  k = 1:K
            
            yq(:,:,n) = yq(:,:,n) + sqrt(p)*H((l-1)*N+1:l*N,n,k)*Fi_matrix(:,pilotIndex(k))';
            
        end
        
        H_bar((l-1)*N+1:l*N,n,:) = yq(:,:,n)*Fi_matrix;
        
    end

end