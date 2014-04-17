function [Ranked,KCDM] = BackCD(x,TarIndx,task_type,kernel_type,reg)
% Markov blanket discovery by backward elimination
% 
% Inputs:
% (1) x = data matrix, where rows are instances and columns are features
% (2) TarIndx = column index of the target
% (3) task_type = 'class' for classification, 'reg' for regression. Uses
%     Kronecker delta kernel for classification.
% (4) kernel_type = 'lin' for linear kernel, 'rbf' for rbf kernel
%     (default='rbf')
% (5) reg = regularization value (default=1E-4)
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

SetDefaultValue(4,'kernel_type','rbf');
SetDefaultValue(5,'reg',1E-4);

[r,c] = size(x);
x = zscore(x);
y = x(:,TarIndx);
x(:,TarIndx) = [];
xindices = 1:c;
xindices(:,TarIndx) = [];

if strcmp(task_type, 'class')
    Ky = kronDel(y);
elseif strcmp(task_type, 'reg')
    doty = y*y';
    Ky = KernelType(doty,kernel_type);
end  
Q=eye(r)-1/r;
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
