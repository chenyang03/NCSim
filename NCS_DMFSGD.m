% Downloaded from http://www.run.montefiore.ulg.ac.be/~liao/DMFSGD.html

function [U,V] = NCS_DMFSGD(X, K, params)
%function [U,V] = NCS_DMFSGD(X, W, params)
%% decentralized matrix factorization by stochastic gradient descent 
%% X: data matrix to be factorized
%% W: Weight matrix with 0 representing missing data

%% parameter initialization
if ~exist('params')
    params.loss = 2; %1: L_2; 2: L_1.
    params.minibatch = 1;
    params.linesearch = 1;

    params.lr = 1e-1;%learning rate
    params.lambda = 1e0;%regularization coeffient
    params.dimension = 8;%number of ranks, called dimension in consistent with vivaldi
 
    params.maxIters = 50;
    params.forceNonNegative = 1;%sign of non-negative factorization
end

lr = params.lr;
lambda = params.lambda;
d = params.dimension;

if isfield(params,'loss')
    loss = params.loss;
else
    loss = 2;
end
if loss~=1 && loss~=2
    loss = 2;
end

minibatch = params.minibatch;
linesearch = params.linesearch;
if linesearch
    maxSearches = 10;
    scalingFactor = 1/2;
    torl = [loss==1,loss==2]*[1e2;1e1];
end
nStartSlow = 1;

maxNorm = 1e+2;
maxv = 1e+4;

if isfield(params,'maxIters')
    maxIters = params.maxIters;
else
    maxIters = 50;
end

if isfield(params,'forceNonNegative')
    NonNegative = params.forceNonNegative;
else
    NonNegative = 1;
end

n = size(X,1);
W = zeros(n, n);

for i=1:n
    randseq = randperm(n);
    for j=1:K
        W(i, randseq(j)) = 1;
    end
end

epoch = length(find(W>0));

%% coordinate initilization
U = rand(n,d);
V = rand(n,d);
Xhat = U*V';

% Uall = zeros(n,d,maxIters);
% Vall = zeros(n,d,maxIters);
% Uall(:,:,1) = U;
% Vall(:,:,1) = V;
% 
% err = zeros(maxIters,1);
% mae = zeros(maxIters,1);
% 
% err(1) = stress(X,Xhat); 
% mae(1) = medianAbsoluteError(X,Xhat); 

%% initialize data structure for each node
s = cputime;
for i = 1 : n %initialize data structure for the neighbor set of each node
    node{i}.neighbor_out.id = find(W(i,:));%ids of nodes in the neighbor set
    node{i}.neighbor_out.k = length(node{i}.neighbor_out.id);%number of neighbors
    node{i}.neighbor_out.x = -1*ones(node{i}.neighbor_out.k,1);%distance measurements of neighbors
    node{i}.neighbor_out.V = zeros(node{i}.neighbor_out.k,d);%coordinates of neighbors
    node{i}.neighbor_out.time = s*ones(node{i}.neighbor_out.k,1);

    node{i}.neighbor_in.id = find(W(:,i));%ids of nodes in the neighbor set
    node{i}.neighbor_in.k = length(node{i}.neighbor_in.id);%number of neighbors
    node{i}.neighbor_in.x = -1*ones(node{i}.neighbor_in.k,1);%distance measurements of neighbors
    node{i}.neighbor_in.U = zeros(node{i}.neighbor_in.k,d);%coordinates of neighbors
    node{i}.neighbor_in.time = s*ones(node{i}.neighbor_in.k,1);
end

%% simulation
k = 1;
ll = epoch*maxIters;
for l = 1 : ll
    skl = cputime;

    % select node i to update
    i = ceil(n*rand(1,1));% pick a node randomly

    % select one neighbor as a reference
    j_id = ceil(node{i}.neighbor_out.k*rand(1,1));
    j = node{i}.neighbor_out.id(j_id);
    i_id = find(node{j}.neighbor_in.id==i);

    ui = U(i,:);
    vj = V(j,:);

    xij = X(i,j);
    if xij>maxv, xij = maxv; end%filter large values
    
    node{i}.neighbor_out.x(j_id) = xij;
    node{i}.neighbor_out.V(j_id,:) = vj;
    node{i}.neighbor_out.time(j_id) = skl;
    node{j}.neighbor_in.x(i_id) = xij;
    node{j}.neighbor_in.U(i_id,:) = ui;
    node{j}.neighbor_in.time(i_id) = skl;

    err_out = node{i}.neighbor_out.x-node{i}.neighbor_out.V*ui';
    err_in = node{j}.neighbor_in.x-node{j}.neighbor_in.U*vj';
    
    if minibatch==1
        jids = 1:node{i}.neighbor_out.k;
        iids = 1:node{j}.neighbor_in.k;
    else
        jids = j_id;
        iids = i_id;
    end

    if loss==1%L_2 loss
        g_ui = - err_out(jids)' * node{i}.neighbor_out.V(jids,:);
        g_vj = - err_in(iids)' * node{j}.neighbor_in.U(iids,:);
    else%L_1 loss
        g_ui = - sign(err_out(jids))' * node{i}.neighbor_out.V(jids,:);
        g_vj = - sign(err_in(iids))' * node{j}.neighbor_in.U(iids,:);
    end

    if k<=nStartSlow && linesearch>0
        lr0 = 1e-5;
    else
        lr0 = lr;
    end

    ui_r = (1-lr0*lambda)*ui-lr0*g_ui;
    vj_r = (1-lr0*lambda)*vj-lr0*g_vj;

    if linesearch && k>nStartSlow
        if loss==1
            e2out = err_out'*err_out;
        else
            e2out = sum(abs(err_out));
        end
        lr0 = lr;
        for r = 1 : maxSearches
            e_outr = node{i}.neighbor_out.x-node{i}.neighbor_out.V*ui_r';
            if loss==1
                e2outr = e_outr'*e_outr;
            else
                e2outr = sum(abs(e_outr));
            end

            if (e2outr<e2out+torl)
                break;
            end

            lr0 = lr0*scalingFactor;
            ui_r = (1-lr0*lambda)*ui-lr0*g_ui;
        end

        if loss==1
            e2in = err_in'*err_in;
        else
            e2in = sum(abs(err_in));
        end
        lr0 = lr;
        for r = 1 : maxSearches
            e_inr = node{i}.neighbor_in.x-node{i}.neighbor_in.U*vj_r';
            if loss==1
                e2inr = e_inr'*e_inr;
            else
                e2inr = sum(abs(e_inr));
            end

            if (e2inr<e2in+torl)
                break;
            end

            lr0 = lr0*scalingFactor;
            vj_r = (1-lr0*lambda)*vj-lr0*g_vj;
        end
    end

    ui = ui_r;
    vj = vj_r;

    nn_u = sqrt(ui*ui');
    if nn_u>maxNorm
        ui = ui*maxNorm/nn_u;
    end
    nn_v = sqrt(vj*vj');
    if nn_v>maxNorm;
        vj = vj*maxNorm/nn_v;
    end

    if NonNegative
        id_neg = find(ui<0);
        if ~isempty(id_neg), ui(id_neg) = 0; end
        id_neg = find(vj<0);
        if ~isempty(id_neg), vj(id_neg) = 0; end
    end

    U(i,:) = ui;
    V(j,:) = vj;

    if mod(l,epoch)~=0, continue; end
    k = l/epoch+1;

    Xhat = U*V';
    Xhat(find(Xhat<0)) = 1;
    
    err(k) = stress(X,Xhat);
    mae(k) = medianAbsoluteError(X,Xhat);
    Uall(:,:,k) = U;
    Vall(:,:,k) = V;

    %disp(sprintf('the current stress after the %dth iteration is %f',k,err(k)))
end

V=V';

% %%
% err = err(1:k);
% mae = mae(1:k); 
% 
% disp('-------------decentralized matrix factorization-------------------');
% disp('----------------------------------------------------------------');
% disp('|         error        | number of iterations | running time(seconds) |');
% disp('----------------------------------------------------------------');
% disp(sprintf('|       %0.5g          |          %d          |          %f          |',err(k),k,cputime-s));
% disp('----------------------------------------------------------------');

return
