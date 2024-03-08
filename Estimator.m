function [Hprecoder_MRT, Hprecoder_FZF] = Estimator(pp,tau,Phii,H_matrix,B,pilot_book,Pilot_indices)

[rows, columns]=size(H_matrix);
LM=rows;
K=columns;

% Finding for P_k ....Define the matrix P_k to represent users by users
% P_k = zeros(K,K);
% for k=1:K
%     % Find the indices of UEs with the same pilot sequence as UE k
%     P_k(:,k) = all(Phii(:,k) == Phii,1)';
% end


%Forming pilot-book matrix
%Generate the pilot book matrix by concatenating the pilot sequences row-wise
pilot_book_matrix = reshape(Phii, tau, 1, K);
pilot_book_matrix = repmat(pilot_book_matrix, 1, tau, 1);


H = ((randn(LM,K)+ 1i*randn(LM,K))/sqrt(2)).*sqrt(B);
W = ((randn(LM,tau)+ 1i*randn(LM,tau))/sqrt(2));
%Yp1 = sqrt(pp)*tau*H*Phii' + sqrt(tau)*W;
Yp = sqrt(pp)*H*Phii' + W;
%yp = Yp1*Phii;
H_bar = Yp*pilot_book;




e = eye(tau);
% Initialize the precoding matrix H_precoder of size LM by K
Hprecoder_MRT = zeros(LM,K);
Hprecoder_FZF = zeros(LM,K);
% Loop through each user (k) to form the precoding matrix
for k = 1:K
   
    %H_bar=Yp1*pilot_book_matrix(:,:,k);
    Hprecoder_MRT(:,k) = H_bar*e(:,Pilot_indices(k));
    Hprecoder_FZF(:,k) = H_bar*(H_bar'*H_bar)^(-1)*e(:,Pilot_indices(k));
end


lambda = zeros(LM,K);
for m=1:LM
    for k=1:K
        den = norm( (B(m,:).^(1/2)).*(Phii(:,k)'*Phii))^2;
        lambda(m,k)=(tau*pp*den +1);
    end
end




%Hhat = (sqrt(pp).*B./(lambda)).*yp;




%Working on new H bar

end