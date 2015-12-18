function rerr = relative_error(estimate, real)
%    rerr=abs(estimate-real)./(min(estimate,real)+1);
    %rerr=abs(log((estimate+0.01)./(real+0.01))./log(2));
    
    estimate = max(estimate, 0); % avoid the negative distance caused by DMF/IDES
    
    rerr=abs(estimate-real)./(min(real,estimate)+0.1);
    %rerr=abs(estimate-real)./(real+0.1);
    %clear estimate;
    real = real(:);
    tmp = find(real<=0);
    rerr = rerr(:);
    %clear real;    
    rerr(tmp) = [];
    
%mask  = (real>0) | (abs(estimate-real)>10);
    
    
%     mask  = (real>0);  
%     rerr=abs(estimate-real)./(real+0.1);
%    
%     rerr= rerr.*mask; 
%     rerr = rerr(:); 
    
    %% square -> line %%
    %a=estimate(:);b=real(:);

%    a(1:10)'
%    b(1:10)'
    

