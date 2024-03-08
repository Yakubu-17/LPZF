%% This function returns the closed-form solution for the sum SE in [ref01]
%ref[01] Anokye et. al., "Full-duplex cell-free massive MIMO with
%low-resolution ADCs, IEEE Trans. Veh. Technol., (Under major review)

function [sumrateu, sumrated] = functionTheoretical(M,Nrx,Ntx,alpha,epsilon,Ku,Kd,bu,bd,sigma,rho,pudB,pddB,ppdB,tau,Phii,Phiid)

pu=10^(pudB/10); pd=10^(pddB/10); pp=10^(ppdB/10);
lambda=zeros(M,Ku);
nu=zeros(M,Kd);zmk=zeros(M,Kd);eta=zeros(M,Kd);
% Solve for lambda
for m=1:M
    
    for k=1:Ku
        c=0;
        for kk=1:Ku
            c=c+ bu(m,kk)*norm(Phii(:,k)'*Phii(:,kk));
        end
       lambda(m,k)=(1 + tau*pp*c)*alpha(m); 
     %   lambda(m,k)=(1 + tau*pp*c); 
    end
    
    for idx=1:Kd
        d=0;
        for iidx=1:Kd
            d=d+ bd(m,iidx)*norm(Phiid(:,idx)'*Phiid(:,iidx));
        end
        nu(m,idx)=(1 + tau*pp*d)*alpha(m);
        zmk(m,idx)=pp*tau*bd(m,idx)^2*(1/nu(m,idx))*alpha(m);
    end
    
end
%power normalization solution
for n=1:M
   for kkx=1:Kd
      eta(n,kkx)= 1/(Ntx*sum(zmk(n,:))); 
   end
end
%Uplink Initialization
signal=zeros(M,Ku);%desired signal power
MUI=zeros(M,Ku); %includes multi-user interference and beamforming error
noise=zeros(M,Ku);%thermal noise
SI=zeros(M,Ku);%residual SI and inter-AP interference
qns = zeros(M,Ku);%Quantization noise
PiCont=zeros(M,Ku,Ku);%pilot contamination effects

%Downlink Initialization
signald = zeros(M,Kd);%desired signal
MUId = zeros(M,Kd);%Multi-user interference
PCd = zeros(M,Kd,Kd);%pilot contamination
for m=1:M
    for i=1:Ku
       signal(m,i) = pp*tau*bu(m,i)^2*alpha(m)^2*(1/lambda(m,i))*Nrx; 
       MUI(m,i) = alpha(m)^3*pp*tau*bu(m,i)^2*(1/lambda(m,i))*Nrx*sum(bu(m,:));
       noise(m,i) = alpha(m)^3*pp*tau*bu(m,i)^2*(1/lambda(m,i))*Nrx;
       
       pcq=zeros(1,Ku); 
       for jj=1:Ku
          if Phii(:,jj)==Phii(:,i)%any(l==Pset(:,k))
              PiCont(m,i,jj)=(pp*tau*Nrx*bu(m,i)*(1/lambda(m,i))*alpha(m)^2*bu(m,jj));
              pcq(jj) = (pp^2*tau^2*Nrx*bu(m,i)^2*(1/lambda(m,i)^2)*alpha(m)^2*bu(m,jj)^2);
          end 
       end
       
       rsik=0;
       for kk=1:Kd
           rsi=0;
           for mm=1:M
               rsi = rsi + alpha(m)^4*pp^2*tau^2*bu(m,i)^2*(1/lambda(m,i))*bd(mm,kk)^2*(1/nu(mm,kk))*sigma(m,mm)*eta(mm,kk)*Nrx*Ntx;
           end
           rsik=rsik+rsi;
       end
       SI(m,i)=rsik;
       
       qns(m,i) = (pu*pp^2*tau^2*bu(m,i)^4*(1/lambda(m,i)^2)*alpha(m)^2*Nrx + pu*pp*tau*bu(m,i)^2*(1/lambda(m,i))*alpha(m)*Nrx*sum(bu(m,:)) ...
           + pu*(sum(pcq)- pcq(i)) + pp*tau*bu(m,i)^2*(1/lambda(m,i))*alpha(m)*Nrx + pd*(rsik/alpha(m)^2))*(alpha(m)-alpha(m)^2);
    end
    
    
    for kdx=1:Kd
        signald(m,kdx) = sqrt(eta(m,kdx))*bd(m,kdx)^2*alpha(m)*(1/nu(m,kdx))*pp*tau*Ntx;
          for kdxx=1:Kd
               if Phiid(:,kdxx)==Phiid(:,kdx)%any(l==Pset(:,k))
                    PCd(m,kdx,kdxx) = sqrt(eta(m,kdxx))*bd(m,kdx)*bd(m,kdxx)*(1/nu(m,kdxx))*alpha(m)*pp*tau*Ntx;
               end
             MUId(m,kdx) = MUId(m,kdx)+  bd(m,kdx)*alpha(m)*pp*tau*Ntx*eta(m,kdxx) * (1/nu(m,kdxx)) * bd(m,kdxx).^2 ;    
          end
    end       
end
rateu=zeros(1,Ku);%initialize for the per user SE
PiCont2=sum(PiCont,1).^2; PiContk=sum(PiCont2,3);
for jdx=1:Ku
    sinr = pu*sum(signal(:,jdx))^2/( pu*sum(MUI(:,jdx)) + sum(noise(:,jdx))...
        + pu*(PiContk(jdx)-sum(signal(:,jdx))^2) + pd*sum(SI(:,jdx)) + sum(qns(:,jdx)));
    rateu(jdx) = log2(1 + sinr);  
  
end
sumrateu = sum(rateu); %sum all the per user SE to obtain the sum SE

rated=zeros(1,Kd); %initialize for the per user SE
pckd=sum(PCd,1).^2; PiContd=sum(pckd,3);
for index=1:Kd
   %PCkd=sum(pckd(:,index,:),3);   
   sigpower = epsilon(index)*pd*sum(signald(:,index))^2;
   sinrd = sigpower/( pd*sum(MUId(:,index)) + pd*(PiContd(index)-sum(signald(:,index))^2)...
          + 1 + pu*sum(rho(index,:))+ pd*(1-epsilon(index))*(sum(signald(:,index))^2) );
   rated(index) = log2(1 + sinrd);
end
sumrated = sum(rated);

end
