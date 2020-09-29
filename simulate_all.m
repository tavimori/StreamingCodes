function final_output = simulate_all(alpha,beta,sim_length, eps)   

t=11;b=8; a=4;
% t=10;b=6; a=3;

del = b-a;
k = t+1-a;
n = t+1+del;
M=4; %Fritchman channel parameter
n_mds = 12; a_mds = 6; 

G = construction_A(t,b,a); %load('G_fk_11_8_4')
% disable fong
% G_fk = construction_fong_khisti(t,b,a);
%For t=b, Martinian trott is just repetition
G_c = construction_C(t,b,a);

% MDS code
m_mds = 4; %Field size for MDS code
P = cauchygen(a_mds,a_mds,m_mds);
I_mat = gf(eye(a_mds),m_mds);
G_mds = [I_mat P];

error_array = zeros(1,6);
counter = 0;
ep = zeros(1,t+k);  % Erasure pattern. Need t values ahead, (k-1) values back.

% at every state, we check if kth packet can be decoded
% the source symbols in kth packets are only corresponding to 
%     1. (k-1) previous coded packet 
%     2. kth packets itself
%     3. t lator coded packet
% thus we made the memory size to be (k-1)+1+t = k+t
ep_full = zeros(n,t+k);  % every encoded packet has length n

ch_state = 0;

% eps_full indicates the erasure distribution
eps_full = zeros(1,n+1);
% with eps prob the whole batch is lost
eps_full(end) = eps;
eps_full(1) = 1- eps;


while(counter<sim_length)  
%     fprintf("%d\n", counter);
     counter = counter+1;
%UNCOMMENT THESE 2 LINES FOR GE CHANNEL     
%     [ep_next, ch_state_out] = GE_sim(eps,alpha,beta,ch_state);
    [ep_next_full, ch_state_out] = batchedGE_sim(eps_full,alpha,beta,ch_state);
    th_uncoded =  (alpha+eps*beta)/(alpha+beta);
%     [ep_next, ch_state_out] = Fritchman_sim(eps,alpha,beta,M,ch_state);
%     [th_uncoded,~] = fritchman_uncoded_error(alpha,beta,M,eps);
    
    ch_state = ch_state_out;
    
    % move first column to the end, like a sliding window
%     ep = circshift(ep,[0,-1]); 
    ep_full = circshift(ep_full,[0,-1]);
    % load new error pattern to the end (overwrite the one out of the window)
%     ep(end) = ep_next;
    ep_full(:,end) = ep_next_full;
    
    % ep: [Y_(t-(k+tau)) .... Y_t]
    
    
    
    
    
    
    %Decode an erased codeword
    % different schemes are tested:
    % [1,2,3] construct A,C, and Fong are streaming codes, so correction of sliding window error is guaranteed
    % if not, we check the linear structure of the receiving matrix
    
    % mt: test repetition condition only
    % mds: test mds receiving is satisfied
    
    % at current situation, we have full information on solving packet k
%     if(ep(k)==1)    %if an erasure happens
    if(ep_full(1:k,k)~=zeros(k,1))  % if erasure happeds
%         mt_success = (ep(k+t)==0); %Due to repetition coding
        
%         if (sw_erasure(ep,t,b,a))
%             sw_flag = 1;
%         else
%             sw_flag = 0;
%         end
        
        % bypass sw check
        if (0)
%         if (sw_erasure(ep,t,b,a))  %If the erasure pattern is a valid sliding window one -- can always decode
            a_success = 1;
            fo_kh_success =1;
            c_success =1;
        else
            a_success = decode_codeword(G,ep_full,n,k,t);
            % disable fong
            % fo_kh_success = decode_codeword(G_fk,ep,n,k,t);
%             fo_kh_success = 0;
            c_success = decode_codeword(G_c,ep_full,n,k,t);
            
%             if ((a_success == 0) || (c_success == 0)) && (sw_flag == 1)
%                 fprintf("error: the sw case is actually unsolvable!!!")
%             end
        end
        
%         if (sum(ep)<=a_mds)
%             mds_success =1;
%         else
%             mds_success = decode_codeword_mds(G_mds,ep,n_mds,k,a_mds,t);
%         end
        fo_kh_success = 1;
        mt_success = 1;
        mds_success =1;
        % the error count, add one to correcponding scheme if decoding is failure
        error_array = error_array+ ones(1,6) - [0, a_success, fo_kh_success,mt_success, c_success,mds_success];
    end
    
%         if (mod(counter,sim_length*0.1)==0)
%             fprintf('-eps=%f-%d%%--',eps,100*counter/sim_length)            
%         end
    
    
end
error_rates = error_array/counter;
final_output = [eps alpha beta sim_length error_rates th_uncoded];
end


function success = decode_codeword(G,ep_full,n,k,t)
success = 1;
for i = 1:k
    % check if ith source symbol can be decoded
    
    % ul is length of G?
    ul = min(t+i,n);
%     ep_i = ep(k-i+1:k-i+ul);
    ep_i = diag(ep_full(:,k-i+1:k-i+ul));
    G_punctured = G(:,find(~ep_i));
    G_punctured_aug = [G_punctured G(:,i)];
    success = success & (rank(G_punctured)==rank(G_punctured_aug));
    if(success==0)
        break;
    end
    
end
end

function success = decode_codeword_mds(G_mds,ep,n_mds,k,a_mds,t)
success = 1;
for i = 1:a_mds
    ul = min(t+i,n_mds);
    ep_i = ep(k-i+1:k-i+ul);
    G_punctured = G_mds(:,find(~ep_i));
    G_punctured_aug = [G_punctured G_mds(:,i)];
    success = success & (rank(G_punctured)==rank(G_punctured_aug));
    
    if(success==0)
        break;
    end
end
end
