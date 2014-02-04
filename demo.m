% examples

[data] = CreateArtificialMB();

[Ranked,KCDM] = BackElimCD(data,4,'lin');

[Ranked,KCDM] = ForSelecCD(data,4,'lin',15);
