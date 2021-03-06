# viscofric2d_evo
MATLAB code to simulate post- and inter-seismic evolution of fault velocity and strain rate in a viscoelastic channel, considering linear, power-law and linear Burger's rheologies. The simulations are consistent over short-term as well as long-term. 

The long-term strain-rates and stress are taken from analytical solutions by Moore and Parsons (2015) - https://doi.org/10.1093/gji/ggv143, while the short-term evolution is solved using a boundary integral scheme. 

We use RK-4th order integration to solve the governing system of ODEs (IVP).

This code relies on 

a modification/version of mesh2d - https://github.com/dengwirda/mesh2d provided by https://people.sc.fsu.edu/~jburkardt/classes/dis_2014/mesh2d/mesh2d.html. This is incorporated into the package so no need to download separately.

This is optional - Scientific colourmaps from Fabio Crameri - https://www.fabiocrameri.ch/colourmaps/

Greens functions calculations from Lambert&Barbot (2016) https://doi.org/10.1002/2016GL070345 and Barbot (2018) https://doi.org/10.1785/0120180058. Also included in this package.
