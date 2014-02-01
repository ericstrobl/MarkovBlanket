% examples

[data] = CreateArtificialMB();

[Ranked,KCDM] = BackElimCD(data,4,'lin');

[Ranked,KCDM] = BackElimCD(data,4,'rbf');
