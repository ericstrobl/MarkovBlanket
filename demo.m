% examples

[data] = CreateArtificialMB();

[Ranked,KCDM] = BackCD(data,4);

[Ranked,KCDM] = ForCD(data,4);

[Ranked,KCDM] = BackCDm(data,4);

[Ranked,KCDM] = ForCDm(data,4,10);
