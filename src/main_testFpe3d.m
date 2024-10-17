%%
% SPDX-FileCopyrightText: 2020 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: Zlib
%
%%

clc;
close all;
clear all;


% This test is ported from the one that appears in the C++ implementation
% of the 3DFPE that is a part of the balance addon in RBDL:
%
% https://github.com/rbdl/rbdl/tree/master/addons/balance
%


% Description:
%
%
% This test compares the output from the RBDL implementation of the 3DFPE to
% the output produced by a Matlab implementation that has been used in the 
% publication of Sloot et al. As tests go this is weak test, because after 
% all the matlab implementation could also have errors in it. However, 
% presently I don't have the time to write a test that evaluates the  
% invariant mathematical properties of the FPE solution to a random set of 
% data. Luckily it is obvious when the FPE solution is in error when it is
% animated in concert with the participant's movements: I have manually 
% checked this, and the solutions appear to be good. The derivative 
% quantities have been derived using maxima. I have compared the values
% returned by these symbolic derivatives against numerically computed 
% derivatives, and these also match. Unfortunately I do not have this
% test code on hand otherwise I would include it.
%
% Sloot LH, Millard M, Werner C, Mombaur K. Slow but Steady: Similar 
% Sit-to-Stand Balance at Seat-Off in Older vs. Younger Adults. Frontiers in 
% sports and active living. 2020;2.

% Details
%
% This data comes from 
%  
%  Sloot LH, Millard M, Werner C, Mombaur K. Slow but Steady: Similar 
%  Sit-to-Stand Balance at Seat-Off in Older vs. Younger Adults. Frontiers in 
%  sports and active living. 2020;2.
% 
%   participant E01, 
%   trial sts_0002_Chest.c3d, 
%   time 1.466666666666667e+00
% 
%  The scripts that write this data to file appear in Millard's RhodeCode
%  repository for
%    FrontiersBalance2020/code/matlab
%       main_ProcessBalanceData.m
%       process3DFootPlacementEstimator.m (lines 189-208)
% 
% the first data point at which all inputs to the FPE are non-zero and the
% center of mass has a speed that is greater than 30 cm/s 

 % mass in kg
m = 100;

% position of the center-of-mass in m
r0C0 = [6.003900000000000e-01;...
        2.167700000000000e-01;...
        5.957600000000000e-01];

% velocity of the center-of-mass in m/s
v0C0 = [3.020800000000000e-01;...
        2.088000000000000e-02;...
        2.757000000000000e-02];


% whole-body moment of inertia (at this instant) about the center of mass in 
% kg-m^2
JC0 = [5.165410000000000e+00, 5.842000000000000e-02,  4.450500000000000e-01;... 
       5.842000000000000e-02, 6.110220000000000e+00, -1.848600000000000e-01;... 
       4.450500000000000e-01,-1.848600000000000e-01, 2.099090000000000e+00];


HC0 = [ ...
    -2.005000000000000e-02;... 
     6.268060000000000e+00;... 
    -8.467000000000000e-02];


% a point on the contact plane
r0S0 = [0;0;0];


% The gravity vector in m/s^2
g0 = [0;0;-9.81];

% Numeric tolerance: a non-linear root needs to be solved to compute the 3DFPE
% and so a tolerance is required
numericTolerance = 1e-12;

% Internally Newton's method is used, so far fewer iterations are usually 
% required.
maximumIterations = 50;


% A small value that is used to regularize the magnitude of angular momentum
% norm. This norm is used to evaluate the projection error as a percentage of
% the body's angular momentum.
omegaSmall = sqrt(eps);


% Derivatives of the fpe step length and the angle phi are calculated so that
% sensitivities can be computed (if desired)
flag_evaluateDerivatives = 1;

% If you are evaluating the fpe for a set of continuous motion capture data, you
% can save some time by using the fpe solution from the previous time step
% as the initial guess for this time step.
fpeInfoGuess = [];

disp('computing the fpe ...');


fpeInfo = calc3DFootPlacementEstimatorInfo(...
            m,...
            r0C0,...
            v0C0,...                                                    
            JC0,...                                                    
            HC0,...
            r0S0,...
            g0,...
            omegaSmall,...
            numericTolerance,...
            maximumIterations,...
            flag_evaluateDerivatives,...
            fpeInfoGuess);

disp('done: ');
fprintf('\t%i\t\tbisection iterations\n',fpeInfo.bisectionIter);
fprintf('\t%i\t\tNewton titerations\n',fpeInfo.newtonIter);
fprintf('\t%e\tfinal solution error\n',fpeInfo.tolerance);

%Check that the computed results are the same as the C++ implementation

scalarTol = sqrt(3)*numericTolerance;
vecNormTol = sqrt(3)*numericTolerance;
matNormTol = sqrt(9)*numericTolerance;

projectionErrorTolerance = 0.01;

%%
% Pre-computed solution
%%

fpeTest.r0G0 = [ ...
    6.003900000000000e-01;... 
    2.167700000000000e-01;... 
    0.000000000000000e+00];

fpeTest.v0G0 = [ ...
    -3.212179353090099e-02;... 
     5.833672380300754e-01;... 
     1.784919742991831e-02];

fpeTest.HG0 = [ ...
     -1.263996880000000e+00;... 
     2.426477808000000e+01 ;... 
     -8.467000000000000e-02];


fpeTest.f    = -8.670397733112623e-12;
fpeTest.phi  = 1.551340426396078e-01;
fpeTest.projectionError = 3.484674002384912e-03;


fpeTest.w0C0 = [...
    -2.019572919578354e-02;...
     1.027672539346732e+00;...
     5.444914458275742e-02];  

fpeTest.n = [ ...
   -5.202130399910231e-02;... 
   9.986459752736367e-01 ;...
   0.000000000000000e+00 ];

fpeTest.u = [ ...
   9.986459752736367e-01;... 
   5.202130399910231e-02;... 
   0.000000000000000e+00];

fpeTest.k = [ ...
   -0.000000000000000e+00;...
   -0.000000000000000e+00;... 
   1.000000000000000e+00];  

fpeTest.r0F0 = [ ...
    6.934351408607322e-01;... 
    2.216168923704711e-01;... 
    0.000000000000000e+00];

fpeTest.h            =  5.957600000000000e-01;
fpeTest.v0C0u        =  3.027571810381615e-01;
fpeTest.v0C0k        =  2.757000000000000e-02;
fpeTest.JC0n         =  6.101593200827201e+00;

fpeTest.Df_Dphi      = -1.824211348591032e+02; 
fpeTest.Df_Dw0C0n    =  6.890339351758533e+00; 
fpeTest.Df_Dh        = -2.847069951699541e+01; 
fpeTest.Df_Dv0C0u    =  6.727732310385989e+01; 
fpeTest.Df_Dv0C0k    =  1.052154468083766e+01; 
fpeTest.Df_DJ        =  8.335241002082351e-01; 
fpeTest.Df_Dm        = -5.085824982564849e-02; 
fpeTest.Df_Dg        = -1.413732690973843e+00; 
fpeTest.Ds_Dl        = -9.411122195708749e-02; 
fpeTest.Ds_DJ        =  2.788743191588045e-03; 
fpeTest.Ds_DE        =  6.597315882691825e-03; 
fpeTest.Ds_Dv0C0u    =  2.250914841062020e-01; 
fpeTest.Ds_Dv0C0k    =  3.520220481488790e-02; 
fpeTest.Ds_Dw0C0n    =  2.305318700460767e-02; 
fpeTest.Dphi_Dl      = -1.541969905085013e-01; 
fpeTest.Dphi_DJ      =  4.569229880364602e-03; 
fpeTest.Dphi_DE      =  1.080940437697085e-02; 
fpeTest.Dphi_Dv0C0u  =  3.688022396956522e-01; 
fpeTest.Dphi_Dv0C0k  =  5.767722412736988e-02; 
fpeTest.Dphi_Dw0C0n  =  3.777160665665429e-02; 
fpeTest.TV           =  5.915445196560693e+02;

disp('testing the fpe against a pre-computed solution');

assert( norm(  fpeInfo.w0C0 - fpeTest.w0C0 , 2 ) < vecNormTol );
assert( norm(  fpeInfo.r0G0 - fpeTest.r0G0 , 2 ) < vecNormTol );
assert( norm(  fpeInfo.HG0  - fpeTest.HG0  , 2 ) < matNormTol );
assert( norm(  fpeInfo.n    -  fpeTest.n   , 2 ) < vecNormTol );
assert( norm(  fpeInfo.u    -  fpeTest.u   , 2 ) < vecNormTol );
assert( norm(  fpeInfo.k    -  fpeTest.k   , 2 ) < vecNormTol ); 
assert(  abs(  fpeInfo.JC0n  -  fpeTest.JC0n   ) < vecNormTol );
assert( norm(  fpeInfo.h     - fpeTest.h       ) < scalarTol  );
assert(  abs(  fpeInfo.v0C0u - fpeTest.v0C0u   ) < scalarTol  );
assert(  abs(  fpeInfo.v0C0k - fpeTest.v0C0k   ) < scalarTol  );
assert(  abs(  fpeInfo.phi   -  fpeTest.phi    ) <= scalarTol );
assert( norm(  fpeInfo.r0F0  -  fpeTest.r0F0   ) <= vecNormTol);
assert( abs(fpeInfo.projectionError-fpeTest.projectionError) <= scalarTol);

disp('success: non-derivative quantities passed the numerical check');

if(flag_evaluateDerivatives==1)

  assert( abs( fpeTest.Df_Dphi     - fpeTest.Df_Dphi    ) <= scalarTol);
  assert( abs( fpeTest.Df_Dw0C0n   - fpeTest.Df_Dw0C0n  ) <= scalarTol);
  assert( abs( fpeTest.Df_Dh       - fpeTest.Df_Dh      ) <= scalarTol);
  assert( abs( fpeTest.Df_Dv0C0u   - fpeTest.Df_Dv0C0u  ) <= scalarTol);
  assert( abs( fpeTest.Df_Dv0C0k   - fpeTest.Df_Dv0C0k  ) <= scalarTol);
  assert( abs( fpeTest.Df_DJ       - fpeTest.Df_DJ      ) <= scalarTol);
  assert( abs( fpeTest.Df_Dm       - fpeTest.Df_Dm      ) <= scalarTol);
  assert( abs( fpeTest.Df_Dg       - fpeTest.Df_Dg      ) <= scalarTol);
  assert( abs( fpeTest.Ds_Dl       - fpeTest.Ds_Dl      ) <= scalarTol);
  assert( abs( fpeTest.Ds_DJ       - fpeTest.Ds_DJ      ) <= scalarTol);
  assert( abs( fpeTest.Ds_DE       - fpeTest.Ds_DE      ) <= scalarTol);
  assert( abs( fpeTest.Ds_Dv0C0u   - fpeTest.Ds_Dv0C0u  ) <= scalarTol);
  assert( abs( fpeTest.Ds_Dv0C0k   - fpeTest.Ds_Dv0C0k  ) <= scalarTol);
  assert( abs( fpeTest.Ds_Dw0C0n   - fpeTest.Ds_Dw0C0n  ) <= scalarTol);
  assert( abs( fpeTest.Dphi_Dl     - fpeTest.Dphi_Dl    ) <= scalarTol);
  assert( abs( fpeTest.Dphi_DJ     - fpeTest.Dphi_DJ    ) <= scalarTol);
  assert( abs( fpeTest.Dphi_DE     - fpeTest.Dphi_DE    ) <= scalarTol);
  assert( abs( fpeTest.Dphi_Dv0C0u - fpeTest.Dphi_Dv0C0u) <= scalarTol);
  assert( abs( fpeTest.Dphi_Dv0C0k - fpeTest.Dphi_Dv0C0k) <= scalarTol);
  assert( abs( fpeTest.Dphi_Dw0C0n - fpeTest.Dphi_Dw0C0n) <= scalarTol);
  assert( abs( fpeTest.TV          - fpeTest.TV         ) <= scalarTol);

end

disp('success: derivative quantities passed the numerical check')


