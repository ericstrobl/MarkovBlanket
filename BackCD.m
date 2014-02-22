function [Ranked,KCDM] = BackCD(x,TarIndx,kernel_type,reg)
% Markov blanket discovery by backward elimination
% 
% Inputs:
% (1) x = data matrix, where rows are instances and columns are features
% (2) TarIndx = column index of the target
% (3) kernel_type = 'lin' for linear kernel, 'rbf' for rbf kernel
%     (default='rbf')
% (4) reg = regularization value (default=1E-6)
% 
% Outputs:
% (1) Ranked = ranking of features in ascending order (least to most likely
%     in Markov blanket)
% (2) KCDM = kernel-based conditional dependence measure when each feature
%     is eliminated
%
% Citation: Strobl EV & Visweswaran S, Markov Blanket Ranking using
% Kernel-based Measures of Conditional Dependence, NIPS Workshop on
% Causality, 2013.

SetDefaultValue(3,'kernel_type','rbf');
SetDefaultValue(4,'reg',1E-4);

[r,c] = size(x);
x = copulaTransform(x);
y = x(:,TarIndx);
x(:,TarIndx) = [];
xindices = 1:c;
xindices(:,TarIndx) = [];

doty = y*y';        
Q=eye(r)-1/r;
Ky = KernelType(doty,kernel_type);
Ky = Q*(Ky)*Q;

toTest = 1:c-2;
KCDM = zeros(1,c-2);
indxDelAcc = zeros(1,c-2);
dotx = x*x';
for t1=1:c-2, 
    KCDMt = zeros(1,length(toTest));
    for t=toTest,
        dotT = dotx - x(:,t)*x(:,t)';
        Kx = KernelType(dotT,kernel_type);
        Kx = Q*Kx*Q + r*reg*eye(r);
        KCDMt(find(t==toTest)) = trace(Ky/Kx);
    end
    KCDMtmin = min(KCDMt);
    KCDM(t1) = KCDMtmin;
    if isnan(KCDMtmin);
        indxDelAcc(t1) = toTest;
        disp(['Eliminating feature: ', num2str(xindices(indxDelAcc(t1)))])
        break
    end
    indxDel = find(KCDMt==KCDMtmin);
    indxDelAcc(t1) = toTest(indxDel(1));
    toTest(indxDel(1)) = [];
    dotx = dotx - x(:,indxDelAcc(t1))*x(:,indxDelAcc(t1))';
    disp(['Eliminating feature: ', num2str(xindices(indxDelAcc(t1)))])
end
Ranked = xindices(indxDelAcc);
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

function SetDefaultValue(position, argName, defaultValue)
% Author: Richie Cotton
if evalin('caller', 'nargin') < position || ...
      isempty(evalin('caller', argName))
   assignin('caller', argName, defaultValue);
end
end

function x = copulaTransform(x)
[r,c] = size(x);
disp('Copula Transform...')
for t=1:c,
[f,x1] = ecdf(x(:,t));
for t1=1:r,
   x(t1,t)=max(f(x1<=x(t1,t)));
end
end
end
