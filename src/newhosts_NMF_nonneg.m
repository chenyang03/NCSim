function [new_w, new_h] =  newhosts_NMF(w,h, newd_out, newd_in, exist_lam, weight_vec)
% [new_w, new_h] =  newhosts_NMF(w,h, newd_out, newd_in, exist_lam)
% w: Nxd, h: dxN. The position vectors of all landmarks
% N: number of landmarks
% M: number of new hosts
% newd_out: Mxn_landmark out distance
% newd_in: n_landmarkxM in distance
% exist_lam: indexes of measured landmarks
% return: new_w: Mxd, new_h: dxM, such that new_w * h = newd
    [N,d] = size(w);
    [M,n_landmark] = size(newd_out);
    % nnls, mldivide
    w = w(exist_lam,:);
    h = h(:,exist_lam);
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%      x = lsqnonneg(C,d)
%      x = lsqnonneg(C,d) returns the vector x that minimizes norm(C*x-d) subject to x >= 0. C and d must be real.
    total=0;total_count=0;
    for i=1:M
       
%         t = lsqnonneg(h', newd_out(i,exist_lam)');
%         new_w(i,:)=t';
%         new_h(:,i)=lsqnonneg(w, newd_in(exist_lam,i));
        t = single(lsqnonneg(h', newd_out(i,exist_lam)'));
        new_w(i,:)=t';
        new_h(:,i)=single(lsqnonneg(w, newd_in(exist_lam,i)));
        
    end
    %fprintf('\n');