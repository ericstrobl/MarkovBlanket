function [Ranked,KDM] = FOHSIC(x,TarIndx,stopNum,task_type,kernel_type)
% Forward selection by HSIC
% 
% Inputs:
% (1) x = data matrix, where rows are instances and columns are features
% (2) TarIndx = column index of the target
% (3) stopNum = number of variables to return
% (4) task_type = 'class' for classification, 'reg' for regression. Uses
%     Kronecker delta kernel for classification.
% (5) kernel_type = 'lin' for linear kernel, 'rbf' for rbf kernel
%     (default='rbf')
%
% Outputs:
% (1) Ranked = ranking of features in descending order (most to least likely
%     in Markov blanket)
% (2) KCDM = kernel-based dependence measure when each feature is eliminated
%
% Citation: Song L, Smola A, Gretton A, Bedo J, Borgwardt K. Feature
% Selection via Dependence Maximization. JMLR, 2007.
%
% Coded by Eric V Strobl, January 2014

SetDefaultValue(5,'kernel_type','rbf');

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

toTest = 1:c-1;
KDM = zeros(1,c-1);
indxDelAcc = zeros(1,c-1);
dotx = zeros(r,r);
for t1=1:stopNum, 
    KDMt = zeros(1,length(toTest));
    for t=toTest,
        dotT = dotx + x(:,t)*x(:,t)';
        Kx = KernelType(dotT,kernel_type);
        Kx = Q*Kx*Q;
        KDMt(find(t==toTest)) = trace(Ky*Kx);
    end
    KDMtmax = max(KDMt);
    KDM(t1) = KDMtmax;
    indxDel = find(KDMt==KDMtmax);
    dotx = dotx + x(:,indxDel(1))*x(:,indxDel(1))';
    indxDelAcc(t1) = toTest(indxDel(1));
    toTest(indxDel(1)) = [];
    disp(['Selecting feature: ', num2str(xindices(indxDelAcc(t1)))])
end
Ranked = xindices(indxDelAcc(1:stopNum));
KDM = KDM(1:stopNum);
end
