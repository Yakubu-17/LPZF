%% This function returns the Monte-Carlo results for the sum SE in [ref01]
%ref[01] Anokye et. al., "Full-duplex cell-free massive MIMO with
%low-resolution ADCs, IEEE Trans. Veh. Technol., 
function [sumrateumc, sumratedmc] = functionSimulationP(M,Nrx,Ntx,alpha,epsilon,Ku,Kd,bu,bd,sigma,rho,pudB,pddB,ppdB,tau,Phii,Phiid,experiment)
pu=10^(pudB/10); pd=10^(pddB/10); pp=10^(ppdB/10);
lambda=zeros(M,Ku);
nu=zeros(M,Kd);zmk=zeros(M,Kd);eta=zeros(M,Kd);
% Solve for lambda and nu: the variance of the channel estimates
for m=1:M
    for k=1:Ku
        c=0;
        for kk=1:Ku
            c=c+ bu(m,kk)*norm(Phii(:,k)'*Phii(:,kk));
        end
        lambda(m,k)=1 + tau*pp*c;       
    end
    for idx=1:Kd
        d=0;
        for iidx=1:Kd
            d=d+ bd(m,iidx)*norm(Phiid(:,idx)'*Phiid(:,iidx));
        end
        nu(m,idx)=1 + tau*pp*d;
        zmk(m,idx)=pp*tau*bd(m,idx)^2*(1/nu(m,idx))*alpha(m);
    end    
end
%power normalization factor
for n=1:M
   for kk=1:Kd
      eta(n,kk)= 1/(Ntx*sum(zmk(n,:))); 
   end
end

H = (randn(M*Nrx,Ku,experiment) + 1i*randn(M*Nrx,Ku,experiment))*sqrt(0.5); %Uplink channel
Hsi=(randn(M*Nrx,M*Ntx,experiment) + 1i*randn(M*Nrx,M*Ntx,experiment))*sqrt(0.5); %residual SI/IAI channel
G = (randn(M*Ntx,Kd,experiment) + 1i*randn(M*Ntx,Kd,experiment))*sqrt(0.5); %Downlink channel
Gudi = (randn(Kd,Ku,experiment) + 1i*randn(Kd,Ku,experiment))*sqrt(0.5); %UL-to-DL interference channel

Bu = repelem(bu,[Nrx],[1]);
Bd = repelem(bd,[Ntx],[1]);
Sigma = repelem(sigma,[Nrx],[Ntx]);
Eta=repelem(eta,[Ntx],[1]);

%Downlink Terms
Termd1 = zeros(Kd,experiment); %initialize the desired signal power
Termd2 = zeros(Kd,experiment); %initialize for multi-user interference
Termd3 = zeros(Kd,experiment); %initialize for UL-to-DL interference
%Uplink Terms
Term1=zeros(Ku,experiment);%initialize for the desired signal
Term2=zeros(Ku,experiment); %initialize for multi-user interference
Term3=zeros(Ku,experiment); %initialize for noise
Term4=zeros(Ku,experiment); %initialize for residual SI and IAI
Term5=zeros(Ku,experiment); %Initialize for quantization noise


for expr=1:experiment
    Hu = H(:,:,expr).*sqrt(Bu); 
    Hhat = ChannelEstimator(tau,pp,Hu,Bu,Phii,alpha);%estimate the channel
    Q = Hsi(:,:,expr).*sqrt(Sigma);
    
    Gd = G(:,:,expr).*sqrt(Bd);
    Ghat = ChannelEstimator(tau,pp,Gd,Bd,Phiid,alpha);%estimate the channel
    Gsi = Gudi(:,:,expr).*sqrt(rho);
    %Uplink Computations
    for ix=1:Ku
        
        Term1(ix,expr) = (Hhat(:,ix)'*Hu(:,ix));%desired signal
        Term2(ix,expr) = sum( abs(Hhat(:,ix)'*Hu).^2); %MUI + beamforming gain
        Term3(ix,expr) = abs(Hhat(:,ix)'*Hhat(:,ix));%Noise
        %Compute residual SI/IAI
        for kx=1:Kd
           Eta2=diag(Eta(:,kx)); 
           Term4(ix,expr) = Term4(ix,expr) + alpha(1)^2*abs(Hhat(:,ix)'*Q*Ghat(:,kx)*Ghat(:,kx)'*Eta2*Q'*Hhat(:,ix));
        end
        
        %Compute Quantization noise (QN)
        Quantsi = zeros(M*Nrx,M*Nrx); %Initialize for QN due to SI/IAI
        for kx=1:Kd
           Eta2=diag(Eta(:,kx)); 
           Quantsi =  Quantsi + Q*Ghat(:,kx)*Ghat(:,kx)'*Eta2*Q';
        end        
        Term5(ix,expr) = Hhat(:,ix)'*diag( diag(pu*Hu*Hu'  + eye(M*Nrx,M*Nrx) + pd*Quantsi) )*Hhat(:,ix);
    
    end
    
    %Downlink Solution
    for jx=1:Kd
        Eta3 = diag(Eta(:,jx));%power normalization
        Termd1(jx,expr) = Gd(:,jx)'*sqrt(Eta3)*Ghat(:,jx); %desired signal
        
        %Compute MUI and pilot contamination
         for jxx=1:Kd
             Eta4 = diag(Eta(:,jxx));
             Termd2(jx,expr) = Termd2(jx,expr) + abs(Gd(:,jx)'*sqrt(Eta4)*Ghat(:,jxx)).^2;
         end
         %Compute UL-to-DL interference
         for ll=1:Ku
             Termd3(jx,expr) = Termd3(jx,expr) + abs(Gsi(jx,ll)*Gsi(jx,ll)');
         end
    end
    
end

Term1x=real(mean(Term1,2));
Term2x=mean(Term2,2);
Term3x=mean(Term3,2);
Term4x=mean(Term4,2);
Term5x=real(mean(Term5,2));
rateumc=zeros(1,Ku); %initialize for rate
for ixx=1:Ku
   sinrmc = alpha(1)^2*pu*Term1x(ixx)^2/( alpha(1)^2*(Term3x(ixx) + pu*(Term2x(ixx)-Term1x(ixx)^2) + pd*Term4x(ixx))...
                 + (alpha(1)-alpha(1)^2)*(Term5x(ixx))  ); %sinr calc
   rateumc(ixx) = log2(1+sinrmc);%rate
    
end
sumrateumc = sum(rateumc); %solve the sumrate


Termd1x = real(mean(Termd1,2));
Termd2x = mean(Termd2,2);
Termd3x = mean(Termd3,2);
ratedmc = zeros(1,Kd);%initialize for rate
for dxx=1:Kd   
    sinrdmc = epsilon(dxx)*pd*Termd1x(dxx)^2/( pd*(Termd2x(dxx)-Termd1x(dxx)^2)...
        + pu*Termd3x(dxx) + 1 + pd*(1-epsilon(dxx))*Termd1x(dxx)^2 );%sinr
    ratedmc(dxx) = log2(1 + sinrdmc); %solve the rate
end
sumratedmc = sum(ratedmc);%compute the sumrate
end
