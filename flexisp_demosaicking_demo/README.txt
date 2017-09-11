Instructions
-----------------------

Dependency
-----------------------
BM3D: http://www.cs.tut.fi/~foi/GCF-BM3D/
LASIP Image Restoration Toolbox: http://www.cs.tut.fi/~lasip/2D/
Gabriel Peyre Toolbox: div.m, getoptions.m grad.m https://www.ceremade.dauphine.fr/~peyre/codes/
Matlab Toolbox: dst.m, dct.m, kaiser.m
Misc: BM3D also depends on check_order.m if you run into any error: https://searchcode.com/codesearch/view/13666439/


BM3D Configuration
-----------------------
CBM3D comes with lots of parameters. Here is the parameter we have used for better result:

transform_2D_HT_name     = 'dct';
transform_2D_Wiener_name = 'dst';

%%%% Hard-thresholding (HT) parameters:
N1                  = 4;
Nstep               = 1;
N2                  = 8;
Ns                  = 49;
tau_match           = 3000*3;
lambda_thr2D        = 0;
lambda_thr3D        = 2.7;
beta                = 2.0;

%%%% Wiener filtering parameters:
N1_wiener           = 4;
Nstep_wiener        = 1;
N2_wiener           = 8*2;
Ns_wiener           = 39;
tau_match_wiener    = 400*3;
beta_wiener         = 2.0;

%%%% Block-matching parameters:
stepFS              = 1;  %% step that forces to switch to full-search BM, "1" implies always full-search
smallLN             = 'not used in np'; %% if stepFS > 1, then this specifies the size of the small local search neighb.
stepFSW             = 1;
smallLNW            = 'not used in np';
thrToIncStep        = 1;  %% used in the HT filtering to increase the sliding step in uniform regions

We also modified CBM3D so that it takes two color spaces, one for hard thresholding and one for Wiener filter. For hard thresholding, we use yCbCr and opp for Wiener.
You need to modify CBM3D.m in order to reproduce the same result. For convenience, we have put a patch file to CBM3D.m in 3rdparty/BM3D that performs the above changes.


Execution
-----------------------
Please execute:
apps/demosaic/demo_demosaic.m
 