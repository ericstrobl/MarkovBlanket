function K = KernelType(dot,kernel_type)
if strcmp(kernel_type,'rbf')
    sig = DetermineSig(dot);
    K = rbf(dot,sig);
elseif strcmp(kernel_type,'lin')
    K = dot;
end
end
