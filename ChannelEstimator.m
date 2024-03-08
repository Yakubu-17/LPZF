%clear
function [Hu_est, Hzf] = ChannelEstimator(tau,pp,H,B,Phii,alpha)
%tau=2;
% [U,S,V]=svd(randn(tau,tau));
% M=20;Nrx=2;Ntx=Nrx;
% MNrx=M*Nrx; MNtx=M*Ntx;
%pp=0.2;
[rows, columns]=size(H);
M=rows;
K=columns;
%Kd=Ku;
%si_canc=0;
% [lsfu,lsfd,lsfusi,lsfudi,~,~]=LSF(M,Ku,Kd,si_canc);
% Bu=repelem(lsfu,[Nrx],[1]);
% Bd=repelem(lsfd,[Ntx],[1]);
%% Pilot Asignment: (random choice)
% Phii=zeros(tau,Ku);
% for k=1:Ku
%     %Point=randi([1,tau]);
%     Point=k;
%     Phii(:,k)=U(:,Point);
% end

%H = ((randn(M,Ku)+ 1i*randn(M,Ku))/sqrt(2)).*sqrt(B);
%taup=tau+tau;
%tau = tau*2;
W = ((randn(M,tau)+ 1i*randn(M,tau))/sqrt(2));
% W = ((randn(M,taup)+ 1i*randn(M,taup))/sqrt(2));
Yp = sqrt(pp)*tau*H*Phii' + sqrt(tau)*W;
yp = Yp*Phii;

lambda = zeros(M,K);%alpha(1)=1;
for m=1:M
    for k=1:K
        %lambda(m,k)=pp*tau*sum(Bu(m,:)*Phii(:,k))
        %num = sqrt(pp)*Bu(m,k);
        den = norm( (B(m,:).^(1/2)).*(Phii(:,k)'*Phii))^2;
        lambda(m,k)=(tau*pp*den +1);
    end
end


% for m=1:M
%     
%     for k=1:K
%         c=0;
%         for kk=1:K
%             c=c+ B(m,kk)*norm(Phii(:,k)'*Phii(:,kk));
%         end
%         lambda(m,k)=1 + tau*pp*c;       
%     end
%     
%     
% end








Hu_est = alpha(1).*(sqrt(pp).*B./(lambda)).*yp;

Hzf = Hu_est*pinv(Hu_est'*Hu_est);

%H_err = Hu-Hu_est;
% H=Hu_est+H_err;
%h=H-Hu;
end





