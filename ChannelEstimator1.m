%clear
function [Hprecoder_MRT, Hprecoder_FZF,Hprecoder_PZF,Hprecoder_PZF_MRT,Hprecoder_PPZF_MRT] = ChannelEstimator1(tau,pp,H,M,Phii,pilot_book,Pilot_indices,S,Y,Z,X,E_S_temp)

[rows, columns]=size(H);
LM=rows;
K=columns;
L=LM/M;

W = ((randn(LM,tau)+ 1i*randn(LM,tau))/sqrt(2));

Yp=zeros(LM,tau); 
for l=1:L
    Yp((l-1)*M+1:M*l,:) = sqrt(pp)*tau*H((l-1)*M+1:M*l,:)*Phii' + sqrt(tau)*W((l-1)*M+1:M*l,:);
end


e = eye(tau);
H_bar = zeros(LM,tau); 
for l=1:L
    H_bar((l-1)*M+1:M*l,:) =Yp((l-1)*M+1:M*l,:)*pilot_book;
end

Hprecoder_MRT = zeros(LM,K);
Hprecoder_FZF = zeros(LM,K);

for l =1:L
    
    for k = 1:K
        
        Hprecoder_MRT((l-1)*M+1:M*l,k)= H_bar((l-1)*M+1:M*l,:)*e(:,Pilot_indices(k));
        Hprecoder_FZF((l-1)*M+1:M*l,k) = H_bar((l-1)*M+1:M*l,:)/(H_bar((l-1)*M+1:M*l,:)'*H_bar((l-1)*M+1:M*l,:))*e(:,Pilot_indices(k));
       
    end
end 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[S,W,Z,M,E_S_temp,tau_S_l2] = functionUEgrouping(LM,K,B,pilotIndex,tau);
eyeM =eye(M);
Hprecoder_PZF = zeros(LM,K);
Hprecoder_PZF_MRT = zeros(LM,K);
Hprecoder_PPZF_MRT = zeros(LM,K);

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
      
       epso = eyetau_S_l(:,find(E_S_l_ind == Pilot_indices(k)));
       Hprecoder_PZF((Z_k_ind(l_Z_k)-1)*M+1:M*Z_k_ind(l_Z_k),k) = (H_bar((Z_k_ind(l_Z_k)-1)*M+1:M*Z_k_ind(l_Z_k),:)*E_S_l)/(E_S_l'*H_bar((Z_k_ind(l_Z_k)-1)*M+1:M*Z_k_ind(l_Z_k),:)'*H_bar((Z_k_ind(l_Z_k)-1)*M+1:M*Z_k_ind(l_Z_k),:)*E_S_l)*epso;
    end


    %Weak UE
    [~,M_k_ind] = find(X(k,:) == 1);
    M_k_ind_length = length(M_k_ind);

    
    
     for  l_M_k = 1:M_k_ind_length
          E_S_l_ind = find(E_S_temp(M_k_ind(l_M_k),:)~=0);

          if sum(E_S_l_ind)~=0
                tau_S_l = length(E_S_l_ind);
                E_S_l = zeros(tau,tau_S_l);
                %eyetau_S_l = eye(tau_S_l );
                
                for r = 1:tau_S_l
                    E_S_l(:,r) = e(:,E_S_l_ind(r));
                end
                B_matrix_temp = H_bar( (M_k_ind(l_M_k)-1)*M+1:M*M_k_ind(l_M_k),:)*E_S_l;
                B_matrix = eyeM - B_matrix_temp/(B_matrix_temp'*B_matrix_temp)*B_matrix_temp';
                
          else
                %B_matrix = eye(1);
                B_matrix = eyeM;
          end

           
           Hprecoder_PZF_MRT( (M_k_ind(l_M_k)-1)*M+1:M*M_k_ind(l_M_k),k) = H_bar( (M_k_ind(l_M_k)-1)*M+1:M*M_k_ind(l_M_k),:)*e(:,Pilot_indices(k));

           Hprecoder_PPZF_MRT( (M_k_ind(l_M_k)-1)*M+1:M*M_k_ind(l_M_k),k) =  B_matrix*H_bar( (M_k_ind(l_M_k)-1)*M+1:M*M_k_ind(l_M_k),:)*e(:,Pilot_indices(k));
    end
   
 end









end





