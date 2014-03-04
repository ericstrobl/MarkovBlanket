function [Ranked,KDM] = BAHSIC(x,TarIndx,task_type,kernel_type)
% Markov blanket discovery by backward elimination
% 
% Inputs:
% (1) x = data matrix, where rows are instances and columns are features
% (2) TarIndx = column index of the target
% (3) task_type = 'class' for classification, 'reg' for regression. Uses
%     Kronecker delta kernel for classification.
% (4) kernel_type = 'lin' for linear kernel, 'rbf' for rbf kernel
%     (default='rbf')
% 
% Outputs:
% (1) Ranked = ranking of features in ascending order (least to most likely
%     in Markov blanket)
% (2) KDM = kernel-based dependence measure when each feature is eliminated
%
% Citation: Song L, Smola A, Gretton A, Bedo J, Borgwardt K. Feature
% Selection via Dependence Maximization. JMLR, 2007.
% 
% Coded by Eric V Strobl, January 2014

SetDefaultValue(4,'kernel_type','rbf');

[r,c] = size(x);
x = copulaTransform(x);
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

toTest = 1:c-2;
KDM = zeros(1,c-2);
indxDelAcc = zeros(1,c-2);
dotx = x*x';
for t1=1:c-2, 
    KDMt = zeros(1,length(toTest));
    for t=toTest,
        dotT = dotx - x(:,t)*x(:,t)';
        Kx = KernelType(dotT,kernel_type);
        Kx = Q*Kx*Q + r*reg*eye(r);
        KDMt(find(t==toTest)) = trace(Ky*Kx);
    end
    KDMtmin = min(KDMt);
    KDM(t1) = KDMtmin;
    if isnan(KDMtmin);
        indxDelAcc(t1) = toTest;
        disp(['Eliminating feature: ', num2str(xindices(indxDelAcc(t1)))])
        break
    end
    indxDel = find(KDMt==KDMtmin);
    indxDelAcc(t1) = toTest(indxDel(1));
    toTest(indxDel(1)) = [];
    dotx = dotx - x(:,indxDelAcc(t1))*x(:,indxDelAcc(t1))';
    disp(['Eliminating feature: ', num2str(xindices(indxDelAcc(t1)))])
end
Ranked = xindices(indxDelAcc);
end
