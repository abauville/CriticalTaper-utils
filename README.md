# CritTaper-utils
A suite of functions and scripts to compute and visualize results of the critical taper theory using the analytical solution of Lehner (1986).The nomenclature of angles follows Lehner (1986). The solution has been benchmarked against the graphical method of	Dahlen (1984).

Please cite this repository and the original paper:
1. Dahlen, F. A. (1984). Noncohesive critical Coulomb wedges: An exact solution. 
Journal of Geophysical Research: Solid Earth, 89(B12):10125–10133.
2. Lehner, F. K. (1986). Comments on "Noncohesive critical Coulomb wedges: an exact
solution" by F. A Dahlen. 	Journal of Geophysical Research, 91(B1):793-796
3. Arthur Bauville. (2019, June 21). CriticalTaper-utils (Version v1.0.1). 
Zenodo. http://doi.org/10.5281/zenodo.3251524

# Content

`basic_script.ipynb` - notebook to compute the critical taper solution and visualize the $\alpha$ (surface angle) vs $\beta$ (basal angle) plot.
`Critical_taper_utils.py` - definition of class Taper().

     import numpy as np
     import matplotlib as plt
     deg = np.pi/180.0
     my_taper = Taper(phi=30.0*deg, phi_b=10.0,
                      Lambda=0.0, Lambda_b=0.0,
                      rho=2500.0, rho_w=1000.0)
     my_taper.computeAlphaVsBeta()
     plt.plot(my_taper.beta_all, my_taper.alpha_all)
     
      
