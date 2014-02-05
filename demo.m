% examples

[data] = CreateArtificialMB();

[Ranked,KCDM] = BackCD(data,4,'rbf');

[Ranked,KCDM] = ForCD(data,4,'lin',10);
