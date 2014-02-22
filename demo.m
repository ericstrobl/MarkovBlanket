% examples

[data] = CreateArtificialMB();

[Ranked,KCDM] = BackCD(data,4);

[Ranked,KCDM] = ForCD(data,4,5);
