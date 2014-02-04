% examples

[data] = CreateArtificialMB();

[Ranked,KCDM] = BackElimCD(data,4,'lin','y');

[Ranked,KCDM] = ForSelecCD(data,4,'lin','y',15);
