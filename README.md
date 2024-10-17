# Description

This repository contains an implementation of the three-dimensional foot-placement-estimator (3DFPE). The 3DFPE was derived by Millard et al. (2012) as an extension to the two-dimensional (2D) FPE derived by Wight et al. (2007). There have been a number of application papers to compare the 3DFPE to sit-to-stand movements (Sloot et al. 2020), to analyze the balance of a humanoid robot (Aller et al. 2021), and to analyze the balance of children with cerebral-palsy (Bruijn et al. 2013). The motivation to derive the 3DFPE came after seeing how accurate the 2DFPE was at predicting human foot placement in the sagittal plane during walking (Millard et al. 2009).


All of the code and files in this repository are covered by the license mentioned in the SPDX file header which makes it possible to audit the licenses in this code base using the ```reuse lint``` command from https://api.reuse.software/. A full copy of the license can be found in the LICENSES folder. To keep the reuse tool happy even this file has a license:

 SPDX-FileCopyrightText: 2024 Matthew Millard <millard.matthew@gmail.com>

 SPDX-License-Identifier: MIT

## Modelling Papers

 - D.Wight DL, E.Kubica, D.Wang (2007). Introduction of the foot placement estimator: A dynamic measure of balance for bipedal robotics. ASME Journal of Computtational and Nonlinear Dynamics, 3(1): 011009 (10 pages). https://doi.org/10.1115/1.2815334

 - M.Millard, D.Wight, J.McPhee, E.Kubica, and D.Wang (2009). Human foot placement and balance in the sagittal plane. ASME Journal of Biomechanical Engineering 131(12). https://doi.org/10.1115/1.4000193

 - M.Millard, J.McPhee, and E.Kubica (2012). Foot Placement and Balance in 3D. ASME Journal of Computational and Nonlinear Dynamics 7(2). https://doi.org/10.1115/1.4005462


## Application Papers

 - F.Aller, M.Harant, S.Sontag, M.Millard, & K.Mombaur (2021). I3SA: The increased step size stability assessment benchmark and its application to the humanoid robot REEM-C. IEEE International Conference on Intelligent Robots and Systems, Prague, Czech Republic, September 27-October 1. https://doi.org/10.1109/IROS51168.2021.9636429

 - S.Bruijn, M.Millard, L. Van Gestel, P.Meyns, I.Jonkers, and K.Desloovere (2013). Gait stability in children with cerebral palsy. Journal of Research in Developmental Disabilities, 34(5). https://doi.org/10.1016/j.ridd.2013.02.011

 - L.H.Sloot*, M.Millard*, C.Werner, & K.Mombaur (2020). Slow but steady: similar sit-to-stand balance at seat-off in older versus younger adults. Frontiers in Sports and Active Living 2(144). https://doi.org/10.3389/fspor.2020.548174 ( * equal first authors)

# Quick start
1. Start Matlab in the src directory
2. Run main_testFpe3d.m. You should see the following output:

```
computing the fpe ...
done: 
    5       bisection iterations
    6       Newton titerations
    1.083578e-13    final solution error
testing the fpe against a pre-computed solution
success: non-derivative quantities passed the numerical check
success: derivative quantities passed the numerical check
```
3. Read the documentation of calc3DFootPlacementEstimatorInfo.m so that you know what inputs to give to the function.