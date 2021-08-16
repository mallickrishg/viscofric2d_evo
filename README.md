# viscofric2d_evo
MATLAB code to simulate post- and inter-seismic evolution of fault velocity and strain rate in a viscoelastic channel. The simulations are consistent over short-term as well as long-term. The long-term strain-rates and stress are taken from analytical solutions by Moore and Parsons (2015), while the short-term evolution is solved using a boundary integral scheme. We use RK-4th order integration to solve the governing ODE (IVP).

This code relies on mesh2d - https://github.com/dengwirda/mesh2d
and scientific colourmaps from Fabio Crameri - https://www.fabiocrameri.ch/colourmaps/

