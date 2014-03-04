% examples

[data] = CreateArtificialMB();

[Ranked,KCDM] = BackCD(data,4,'reg');

[Ranked,KCDM] = ForCD(data,4,5,'reg');

[Ranked,KDM] = BAHSIC(data,4,'reg');

[Ranked,KDM] = FOHSIC(data,4,5,'reg');
