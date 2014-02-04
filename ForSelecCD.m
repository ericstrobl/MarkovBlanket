function [Ranked,KCDM] = ForSelecCD(x,TarIndx,kernel_type,stopNum)
% Inputs:
% (1) x = data matrix, where rows are instances and columns are features
% (2) TarIndx = column index of the target
% (3) kernel_type = 'lin' for linear kernel, 'rbf' for rbf kernel
% (4) stopNum = number of variables to return
%
% Outputs:
% (1) Ranked = ranking of features in descending order (most to least likely
%     in Markov blanket)
% (2) KCDM = kernel-based conditional dependence measure when each feature
%     is eliminated
%
% Citation: Strobl EV & Visweswaran S, Markov Blanket Ranking using
% Kernel-based Measures of Conditional Dependence, NIPS Workshop on
% Causality, 2013.

[r,c] = size(x);
y = x(:,TarIndx);
x(:,TarIndx) = [];
xindices = 1:c;
xindices(:,TarIndx) = [];
x=zscore(x);
y=zscore(y);

doty = y*y';
Q=eye(r)-1/r;
Ky = KernelType(doty,kernel_type);
Gy = Q*Ky*Q;

toTest = 1:c-2;
KCDM = zeros(1,c-2);
indxDelAcc = zeros(1,c-2);
dotx = zeros(r,r);
for t1=1:stopNum, 
    KCDMt = zeros(1,length(toTest));
    for t=toTest,
        dotT = dotx + x(:,t)*x(:,t)';
        Kx = KernelType(dotT,kernel_type);
        Gx = Q*Kx*Q + r*0.01*eye(r);
        KCDMt(find(t==toTest)) = trace(Gy/Gx);
    end
    KCDMtmin = min(KCDMt);
    KCDM(t1) = KCDMtmin;
    indxDel = find(KCDMt==KCDMtmin);
    dotx = dotx + x(:,indxDel(1))*x(:,indxDel(1))';
    indxDelAcc(t1) = toTest(indxDel(1));
    toTest(indxDel(1)) = [];
    disp(['Selecting feature: ', num2str(xindices(indxDelAcc(t1)))])
end
Ranked = xindices(indxDelAcc(1:stopNum));
KCDM = KCDM(1:stopNum);
end

function sig = DetermineSig(dot)
Dis = PairWiseDistance(dot);
xDis = triu(Dis);
xDis = xDis(:);
xDis(xDis==0)=[];
sig = median(xDis);
if isnan(sig),
    sig = 1;
end
end

function Dis = PairWiseDistance(dot)
n=size(dot,1);
Md = diag(dot);
Dis =repmat(Md',n,1)+repmat(Md,1,n)-2*(dot);
Dis = sqrt(Dis);
end

function K = rbf(Dot,sig)
n=size(Dot,1);
K=Dot/sig^2;
d=diag(K);
K=K-ones(n,1)*d'/2;
K=K-d*ones(1,n)/2;
K=exp(K);
end

function K = KernelType(dot,kernel_type)
if strcmp(kernel_type,'rbf')
    sig = DetermineSig(dot);
    K = rbf(dot,sig);
elseif strcmp(kernel_type,'lin')
    K = dot;
end
end

