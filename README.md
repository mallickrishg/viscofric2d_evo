# viscofric2d_evo
MATLAB code to simulate post- and inter-seismic evolution of a two-dimensional strike-slip plate boundary that are mechanically consistent over the short-term (10-100 years) and long-term (10 ka - 1 Ma).<br> 
The code solves the equilibrium equations to compute fault velocity and strain rate in a viscoelastic channel. The routines can handle arbitrary spatial heterogeneity for viscous parameters corresponding to Newtonian, power-law and linear Burger's rheologies while the elastic elements correspond to a homogeneous half-space.

AUTHOR: Rishav Mallick (Seismological Laboratory, Caltech)

The long-term strain-rates and stress are taken from analytical solutions by Moore and Parsons (2015) - https://doi.org/10.1093/gji/ggv143, while the short-term evolution is solved using a boundary integral scheme. 

We use RK-4th order integration to solve the governing system of ODEs (IVP).

This code relies on 

a modification/version of mesh2d - https://github.com/dengwirda/mesh2d provided by https://people.sc.fsu.edu/~jburkardt/classes/dis_2014/mesh2d/mesh2d.html. This is incorporated into the package so no need to download separately.

This is optional - Scientific colourmaps from Fabio Crameri - https://www.fabiocrameri.ch/colourmaps/

Greens functions calculations from Lambert&Barbot (2016) https://doi.org/10.1002/2016GL070345 and Barbot (2018) https://doi.org/10.1785/0120180058. Also included in this package.

To use this work, please cite:<br>
Mallick, R., Lambert, V., & Meade, B. (2022). On the choice and implications of rheologies that maintain kinematic and dynamic consistency over the entire earthquake cycle. Journal of Geophysical Research: Solid Earth, 127, e2022JB024683. https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022JB024683
