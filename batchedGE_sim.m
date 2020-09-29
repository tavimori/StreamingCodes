function [erasure, out_state] = batchedGE_sim(eps,alpha,beta,in_state)
%batchedGE_sim generalization of 
%   alpha - prob(good->burst)
%   beta - prob(burst->good)
%   eps - eps(i)=prob(# random erasures = i)

% erasure = 0;
out_state = in_state; 
rn = rand(1,2); 

batch_size = length(eps) - 1; % eps has M+1 entries


if (in_state == 0) % Random erasure state
    % sample number of erasures
    
    num_erasure = randsample(0:batch_size,1,true,eps);
%     fprintf("re: %d\n", num_erasure);
    erasure_ind = randsample(batch_size, num_erasure, false);
    erasure = zeros(batch_size, 1);
    erasure(erasure_ind) = 1;
    if (rn(2)<alpha)
        out_state =1; 
    end
else                %Burst erasure state
%     fprintf("be\n");
    erasure = ones(batch_size, 1);
    if(rn(2)<beta)
        out_state = 0;
    end
end


end
