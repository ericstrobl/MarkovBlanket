function [Ranked,KCDM] = ForSelecCD(x,TarIndx,kernel_type,stopNum)
% Inputs:
% (1) x = data matrix, where rows are instances and columns are features
% (2) TarIndx = column index of the target
% (3) kernel_type = 'lin' for linear kernel, 'rbf' for rbf kernel
% (4) stopNum = specified number of variables to return
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
if strcmp(kernel_type,'rbf')
    sigy = DetermineSig(doty);
	Kyt = rbf(doty,sigy);
elseif strcmp(kernel_type,'lin')
    Kyt = doty;
end
Q=eye(r)-1/r;
Gy = Q*Kyt*Q;

toTest = 1:c-2;
KCDM = zeros(1,c-2);
indxDelAcc = zeros(1,c-2);
dotx = x*x';
for t1=1:stopNum, 
    KCDMt = zeros(1,length(toTest));
    for t=toTest,
        dotT = dotx - x(:,t)*x(:,t)';
        if strcmp(kernel_type,'rbf')
            sigx = DetermineSig(dotT);
            Kx = rbf(dotT,sigx);
        elseif strcmp(kernel_type,'lin')
            Kx = dotT;
        end
        Gx = Q*Kx*Q + 0.01*eye(r);
        KCDMt(find(t==toTest)) = trace(Gy/Gx);
    end
    KCDMtmax = max(KCDMt);
    KCDM(t1) = KCDMtmax;
    indxDel = find(KCDMt==KCDMtmax);
    indxDelAcc(t1) = toTest(indxDel(1));
    toTest(indxDel(1)) = [];
    disp(['Eliminating feature: ', num2str(xindices(indxDelAcc(t1)))])
end
Ranked = xindices(indxDelAcc(1:stopNum));
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

