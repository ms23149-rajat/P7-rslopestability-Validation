/*
##############################################################################
#
# MODULE:       r.slope.stability.main.c
#
# AUTHOR:       Martin Mergili
# CONTRIBUTORS: Massimiliano Alvioli, Ivan Marchesini, Markus Neteler
#
# PURPOSE:      The slope stability model - Main script
#
# COPYRIGHT:    (c) 2009 - 2021 by the author
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2009 - 2021 by the BOKU University, Vienna
#               (c) 2015 - 2021 by the University of Vienna
#               (c) 2000 - 2021 by the GRASS Development Team
#
# VERSION:      20210625 (25 June 2021)
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
##############################################################################
*/


// Libraries

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grass/gmath.h>
#include <grass/gis.h>
#include <grass/gprojects.h>
#include <grass/glocale.h>
#include <grass/raster.h>
#include <fcntl.h>
#include <time.h>
#include <unistd.h>
#include <grass/segment.h>

#include <string.h>
#include <unistd.h>
#include <sys/stat.h>


// Operational variables

char *mapset1, **strtodhlp=NULL, path[1200], wkdir[1000], sparam[200], slayers[102], xchar[500], *yint, *TEMPDIR, *MAPSET;


//Constant and user-defined parameters

int NUM_MCSIM, TRUNCATED, LAYERS, SEEPAGE, PSEUDOSTATIC, NEWMARK, NEWMARKREF, INFMOD, ELLIPS, CLASSES, CLASSINF, REVISED, SUMMARY, MAXLAYERS, USESEG, PF, PFDEPTH,
    ADDITIONAL, PFSAMPLING, pf_nsteps_csoil, pf_nsteps_phi, pf_nsteps_depth, pdftype_csoil, pdftype_phi, pdftype_depth, xint;

float PF_STDEVS, UNDEFMIN, UNDEFPLUS, UNDEFLYR, CHANGEDM, PI, LENGTHMAX, LENGTHMIN, WIDTHMAX, WIDTHMIN, DEPTHMIN, DEPTHMAX, ZBMAX, ZBMIN,
    CENTERX, CENTERY, GRAVITY, GAMMA_WATER, ELLAMAX, ELLAMIN, ELLBMAX, ELLBMIN, ELLCMAX, ELLCMIN, RAT_GEOM_MIN, RAT_GEOM_MAX, ELLDENS,
    ELL_ALPHADEF, ELL_BETADEF, MEMORYSEG, pf_depth_mean, pf_depth_stdev, pf_depth_min, pf_depth_max, seiscoef, eqdepth, magnitude, pars1, pars2, regpar1, regpar2, newmarkcoef[5], fact_revised;

SEGMENT seg_soilclass, seg_elev, seg_pga, seg_distance, seg_newmarkref, seg_gwdepth, seg_depthl, seg_dzdx, seg_dzdy, seg_ctrl_lyr, seg_num_lyr;


// Control and input variables

int i, ipar, j, jmax, j2, j3, j4, k, l, m, mx, n, nx, o, p, me, ne, x, y, ctrl_continue, ctr_valid, pathtest, time_start, time_end,
    time_stepstart, time_stepend, nlinesl, lbtm, lbtm2, nlayersmax, kmax, ctrl_valtopo, linit, lloc, lnew, lglob, num_cellval,
    lbtm_max, ctrl_lyr, num_lyr, nsegs, nsegrc, nsoilclass, nlayer[102], nlayersl, nly, pscs, pf_ctr_max, pf_testctr, pf_nsteps, pparias[5],
    iarias, iariasmin, iariasmax, jarias, narias; //jctrl;

float csz, bd_west, bd_north, time_elapsed, cell_change, my_depth, my_csoil, my_phi, sigma_depth, sigma_csoil, sigma_phi, mode,
    pf_depth_int;


// Input files

FILE *fparam1, *fparam2, *fparam3, *flayersf, *flayers;


// Timer

time_t timer;


// Variable lists

struct GModule *module;
struct Cell_head cellhd;

/*int pf_comp ( const int * f, const int * s ) { // sorting arrays
    if ( *f > *s ) return 1;
    if ( *f < *s ) return -1;
    return 0;
}*/

int pf_comp (const void * f, const void * s) {
   return ( *(int*)f - *(int*)s );
}

int **G_alloc_imatrix(int mrows, int mcols) { // allocating integer arrays
    int **mm;
    int mi;
    mm = (int **)G_calloc(mrows, sizeof(int *));
    mm[0] = (int *)G_calloc(mrows * mcols, sizeof(int));
    for (mi = 1; mi < mrows; mi++)
    mm[mi] = mm[mi - 1] + mcols;
    return mm;
}

float **G_alloc_dmatrix(int mrows, int mcols) { // allocating float arrays
    float **mm;
    int mi;
    mm = (float **)G_calloc(mrows, sizeof(float *));
    mm[0] = (float *)G_calloc(mrows * mcols, sizeof(float));
    for (mi = 1; mi < mrows; mi++)
    mm[mi] = mm[mi - 1] + mcols;
    return mm;
}

float ***alloc_dmatrix3 (int mrows, int mcols, int mdepths) { // function for allocating 3D float arrays
    float ***mm;
    int mi, mj;
    mm = (float***)calloc(mrows, sizeof(float**));

    for (mi=0; mi<mrows; mi++) {
        mm[mi] = (float**)calloc(mcols, sizeof(float*));

        for (mj=0; mj<mcols; mj++) {
            mm[mi][mj] = (float*)calloc(mdepths, sizeof(float));
        }
    }
    return mm;
}

void free_dmatrix3 (float ***mm,int mi,int mj)
{
    int i,j;

    for(i=0;i<mi;i++)
    {
        for(j=0;j<mj;j++)
        {
            free(mm[i][j]);
        }
        free(mm[i]);
    }
    free(mm);
}

int fiparam ( FILE *gfparam ) { // reading integer parameters

    int gparam;
    char gsparam[100], **gstrtodhlp=NULL;
    if( fgets(gsparam, 100, gfparam) == NULL ) return -1;
    gparam=strtod(gsparam, gstrtodhlp);
    return gparam;
}

float fdparam ( FILE *gfparam ) { // reading float parameters

    float gparam;
    char gsparam[100], **gstrtodhlp=NULL;
    if( fgets(gsparam, 100, gfparam) == NULL ) return -1;
    gparam=strtod(gsparam, gstrtodhlp);
    return gparam;
}

char *fcparam ( FILE *gfparam ) { // reading character parameters

    char *gsparam = (char*) malloc ( sizeof(char) * 1000 );
    if( fgets(gsparam, 1000, gfparam) == NULL ) return FALSE;
    gsparam[strnlen(gsparam, 1010) - 1] = '\0';
    return gsparam;
}

SEGMENT finrasti ( char *gm, float gundef ) { // input of integer GRASS raster maps

    int gf, gx, gy, gin;
    CELL *gc;
    SEGMENT gseg_in;

    Segment_open (&gseg_in, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(int), nsegs);

    gc = Rast_allocate_c_buf();
    gf = Rast_open_old ( gm, MAPSET );

    for ( gx = 0; gx < m; gx++ ) {

        Rast_get_c_row ( gf, gc, gx );

        for ( gy = 0; gy < n; gy++ ) {

            if ( !Rast_is_c_null_value ( gc + gy ) ) gin = ( int ) gc[gy];
            else gin = (int)gundef;
            if ( gin < 1 ) gin = (int)gundef;
            Segment_put(&gseg_in, &gin, gx, gy);
	}
    }
    Segment_flush(&gseg_in);
    G_free ( gc );
    Rast_close ( gf );

    return gseg_in;
}

SEGMENT finrastd ( char *gm, float gundef ) { // input of float GRASS raster maps

    int gf, gx, gy;
    float gin;
    DCELL *gc;
    SEGMENT gseg_in;

    Segment_open (&gseg_in, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs);

    gc = Rast_allocate_d_buf();
    gf = Rast_open_old ( gm, MAPSET );

    for ( gx = 0; gx < m; gx++ ) {

        Rast_get_d_row ( gf, gc, gx );

        for ( gy = 0; gy < n; gy++ ) {

            if ( !Rast_is_d_null_value ( gc + gy ) ) gin = ( float ) gc[gy];
            else gin = gundef;
            Segment_put(&gseg_in, &gin, gx, gy);
	}
    }

    Segment_flush(&gseg_in);
    G_free ( gc );
    Rast_close ( gf );

    return gseg_in;
}

void foutrasti ( char *gm, SEGMENT gseg_v ) { // output of integer GRASS raster maps

    int gf, gx, gy;
    int gv;
    CELL *gc;

    if ( gm != NULL ) {
        gc = Rast_allocate_c_buf();
	gf = Rast_open_c_new( gm );

        for ( gx = 0; gx < m; gx++ ) {

            if ( gm != NULL ) {
	        for ( gy = 0; gy < n; gy++ ) {
                    Segment_get(&gseg_v, &gv, gx, gy);
                    gc[gy] = (CELL) gv;
                }
	        Rast_put_c_row( gf, gc );
            }
        }
        G_free( gc );
        Rast_close( gf );
    }

    return;
}

void foutrastd ( char *gm, SEGMENT gseg_v, int gdigits ) { // output of float GRASS raster maps

    int gf, gx, gy;
    float gv;
    DCELL *gc;

    if ( gm != NULL ) {
        gc = Rast_allocate_d_buf();
        gf = Rast_open_fp_new( gm );

        for ( gx = 0; gx < m; gx++ ) {

            if ( gm != NULL ) {
	        for ( gy = 0; gy < n; gy++ ) {
                    Segment_get(&gseg_v, &gv, gx, gy);
                    gv = round( pow( 10, gdigits ) * gv ) / pow( 10, gdigits );
                    gc[gy] = (DCELL) gv;
                }
	         Rast_put_d_row( gf, gc );
            }
        }
        G_free( gc );
        Rast_close( gf );
    }

    return;
}

float *fnewmark ( float fos, float beta, float magnitude, float distance, float eqdepth, float pars1, float pars2, float *newmarkcoef, float newmarkref, 
    int iariasmin, int iariasmax, int *pparias, int narias ) { // Newmark displacement

    int jarias, iarias, iall; 
    float *gnewmark, pnewmark_raw, arias[5], parias0[5], parias[5], ariastest, newmarktest, newmarkmean, critacc;
    gnewmark = (float*) calloc( 3, sizeof(float));

    if ( fos > 1 ) {

        critacc = ( fos - 1 ) * 9.81 * sin( beta );

        pnewmark_raw = 0;
        newmarkmean = 0;
        iall = 0;
        for ( jarias = 0; jarias <5; jarias++ ) parias0[jarias] = 0;
        for ( iarias = iariasmin; iarias < iariasmax; iarias++ ) {
        
            if ( pparias[0] == 1 ) {
                arias[0] = 0.98 * magnitude - 1.35 * log10( pow( pow( distance, 2) + pow( eqdepth, 2 ), 0.5 )) - 4.9;
                parias[0] = 1 / (float)( iariasmax - iariasmin );
            }
            if ( pparias[1] == 1 ) {
                arias[1] = magnitude - 2 * log10( pow( pow( distance, 2) + pow( eqdepth, 2 ), 0.5 )) - 4.1 - 3 * 0.44 + 6 * 0.44 / (float)( iariasmax - iariasmin ) * ( iarias - 0.5 );
                ariastest = ( -3 * 0.44 + 6 * 0.44 / (float)( iariasmax - iariasmin ) * iarias ) / 0.44;
                parias[1] = 1 / ( exp( -358 * ariastest / 23 + 111 * atan( 37 * ariastest / 294 )) + 1 ) - parias0[1];
                parias0[1] = 1 / ( exp( -358 * ariastest / 23 + 111 * atan( 37 * ariastest / 294 )) + 1 );
            }
            if ( pparias[2] == 1 ) {
                arias[2] = 0.601 * magnitude - log10( distance ) - 0.011 * distance - 2.659 
                    - 3 * 0.43 + 6 * 0.43 / (float)( iariasmax - iariasmin ) * ( iarias - 0.5 );
                ariastest = ( -3 * 0.43 + 6 * 0.43 / (float)( iariasmax - iariasmin ) * iarias ) / 0.43;
                parias[2] = 1 / ( exp( -358 * ariastest / 23 + 111 * atan( 37 * ariastest / 294 )) + 1 ) - parias0[2];
                parias0[2] = 1 / ( exp( -358 * ariastest / 23 + 111 * atan( 37 * ariastest / 294 )) + 1 );
            }
            if ( pparias[3] == 1 ) {
                arias[3] = -4.066 + 0.911 * magnitude - 1.818 * log10( pow( pow( distance, 2 ) + 28.09, 0.5 )) + 0.139 * pars1 + 0.244 * pars2 
                    - 3 * 0.397 + 6 * 0.397 / (float)( iariasmax - iariasmin ) * ( iarias - 0.5 );
                ariastest = ( - 3 * 0.397 + 6 * 0.397 / (float)( iariasmax - iariasmin ) * iarias ) / 0.397;
                parias[3] = 1 / ( exp( -358 * ariastest / 23 + 111 * atan( 37 * ariastest / 294 )) + 1 ) - parias0[3];
                parias0[3] = 1 / ( exp( -358 * ariastest / 23 + 111 * atan( 37 * ariastest / 294 )) + 1 );
            }
            if ( pparias[4] == 1 ) {
                arias[4] = 0.810 * magnitude - log10( pow( pow( distance, 2) + pow( eqdepth, 2 ), 0.5 )) 
                - 0.02 * pow( pow( distance, 2) + pow( eqdepth, 2 ), 0.5 ) - 3.88;
                parias[4] = 1 / (float)( iariasmax - iariasmin );
            }
            
            for ( jarias = 0; jarias <5; jarias++ ) {
        
                if ( pparias[jarias] == 1 ) {
        
                    gnewmark[0] = newmarkcoef[0] * arias[jarias] + newmarkcoef[1] * critacc + newmarkcoef[2] * arias[jarias] * critacc + newmarkcoef[3];
                    newmarktest = ( gnewmark[0] - log10( newmarkref )) / newmarkcoef[4];
                    pnewmark_raw += parias[jarias] / narias / ( exp( -358 * newmarktest / 23 + 111 * atan( 37 * newmarktest / 294 )) + 1 );
                    if ( gnewmark[0] > log10( newmarkref )) gnewmark[0] = log10( newmarkref );
                    newmarkmean += pow( 10, gnewmark[0] );
                    iall += 1;
                }
            }
        }
        
        gnewmark[0] = pow( 10, gnewmark[0] );
        gnewmark[1] = pnewmark_raw;
        gnewmark[2] = newmarkmean / iall;

    } else { gnewmark[0] = newmarkref; gnewmark[1] = 1; gnewmark[2] = newmarkref; }

    return gnewmark;
}


// Calling main function

int main ( int argc, char *argv[] ) {


// Defining constants

    PI = 3.1415926536;
    GRAVITY = 9.81;
    GAMMA_WATER = 9810;
    UNDEFPLUS = 9999;
    UNDEFMIN = -9999;
    UNDEFLYR = 0;
    CHANGEDM = 1000;

// Defining variables

    struct Option *input_opt1, *input_opt2, *input_opt3;

// Initializing GIS functions

    G_gisinit( argv[0] );

    module = G_define_module();
    G_add_keyword(_("raster"));
    G_add_keyword(_("factor of safety"));
    G_add_keyword(_("limit equilibrium model"));
    G_add_keyword(_("sliding surface"));
    G_add_keyword(_("slope failure probability"));
    G_add_keyword(_("slope stability"));
    module->description =
	_("The slope stability model");

    input_opt1 = G_define_standard_option(G_OPT_F_BIN_INPUT);
    input_opt1->key = "input1";
    input_opt1->description = _("Parameter file 1");
    input_opt1->required = NO;
    input_opt1->answer = "TEMP/rtemp0/param1.txt";

    input_opt2 = G_define_standard_option(G_OPT_F_BIN_INPUT);
    input_opt2->key = "input2";
    input_opt2->description = _("Parameter file 2");
    input_opt2->required = NO;
    input_opt2->answer = "TEMP/rtemp0/param2.txt";

    input_opt3 = G_define_standard_option(G_OPT_F_BIN_INPUT);
    input_opt3->key = "input3";
    input_opt3->description = _("Parameter file 3");
    input_opt3->required = NO;
    input_opt3->answer = "TEMP/rtemp0/param3.txt";

    if (G_parser(argc, argv))
        exit(EXIT_FAILURE);

    G_get_set_window( &cellhd );


// Setting spatial extent and cell size

    csz=cellhd.ew_res;
    bd_west=cellhd.west;
    bd_north=cellhd.north;
    n=cellhd.cols;
    m=cellhd.rows;


// Processing arguments

    yint = (char*) malloc( sizeof(char) * 100 );   

    char* buffer;

    if ( argc == 0 ) {
        
        exit(0);
        
    } else if ( argc > 1 ) {
        
        sprintf( wkdir, "%s", argv[2] );
        xint = 0;
            
    } else {
        
        buffer=getcwd(NULL, 0);
        sprintf(wkdir, "%s/TEMP/rtemp", buffer );
        xint=0;
        yint=getenv("XINT");
        if ( yint != NULL ) xint=atoi(yint);
        free(buffer);
    }

    G_init_tempfile(); // initializing generation of temporary files
    timer=time(NULL);
    srand ( timer ); // seed for random numbers


// Opening parameter files

    fparam1=fopen(input_opt1->answer, "r"); // opening parameter file 1
    fparam2=fopen(input_opt2->answer, "r"); // opening parameter file 2
    fparam3=fopen(input_opt3->answer, "r"); // opening parameter file 3


// Reading parameter file 1

    if ( fparam1 != NULL ) { if( fgets(sparam, 12, fparam1) == NULL ) return -1; } // reading header of parameter files
    if ( fparam2 != NULL ) { if( fgets(sparam, 12, fparam2) == NULL ) return -1; }
    if ( fparam3 != NULL ) { if( fgets(sparam, 12, fparam3) == NULL ) return -1; }

    TEMPDIR = (char*) malloc( sizeof(char) * 1000 );
    MAPSET = (char*) malloc( sizeof(char) * 1000 );

    if ( fparam1 != NULL ) {

    TEMPDIR = fcparam ( fparam1 ); // temporary directory
    MAPSET = fcparam ( fparam1 ); // mapset
    CLASSES = fiparam ( fparam1 ); // application of class-based model
    INFMOD = fiparam ( fparam1 ); // application of infinite slope stability model (0=no, 1=yes)
    LAYERS = fiparam ( fparam1 ); // application of layer-based model
    REVISED = fiparam ( fparam1 ); // application of revised Hovland model
    TRUNCATED = fiparam ( fparam1 ); // truncation of ellipsoid
    ADDITIONAL = fiparam ( fparam1 ); // writing additional output
    SEEPAGE = fiparam ( fparam1 ); // seepage direction
    PSEUDOSTATIC = fiparam ( fparam1 ); // pseudo-static seismic slope stability analysis
    NEWMARK = fiparam ( fparam1 ); // Newmark seismic slope stability analysis
    NEWMARKREF = fiparam ( fparam1 ); // application of critical displacement map

    if ( PSEUDOSTATIC == 1 ) seiscoef = fdparam ( fparam1 ); // seismic coefficient
    if ( NEWMARK == 1 ) {
        eqdepth = fdparam ( fparam1 ); // earthquake depth
        magnitude = fdparam ( fparam1 ); // moment magnitude
        pars1 = fdparam ( fparam1 ); // parameter S1
        pars2 = fdparam ( fparam1 ); // parameter S2
        
        for ( i = 0; i < 5; i++ ) {
            newmarkcoef[i] = fdparam ( fparam1 ); //coefficients for Newmark displacement
        }     
        
        narias = 0;
        for ( i = 0; i < 5; i++ ) {
            pparias[i] = fiparam ( fparam1 ); //equation for arias intensity
            narias += pparias[i];
        }
    }

    if ( CLASSES == 1 || LAYERS == 1 ) INFMOD = 0; // excluding infinite slope stability mode for ellipsoid-based modes
    if ( CLASSES == 1 ) LAYERS = 0; // excluding layer-based mode for soil class-based mode
    if ( LAYERS == 1 ) INFMOD = 0; // excluding infinite slope stability model for layer-based mode
    if ( LAYERS == 0 && SEEPAGE == 2 ) SEEPAGE = 1; // excluding layer-parallel seepage for soil class-based modes
    if ( CLASSES == 1 || INFMOD == 1 ) CLASSINF = 1; else CLASSINF = 0; // control for application of soil class-based modes
    if ( CLASSES == 1 || LAYERS == 1 ) ELLIPS = 1; else ELLIPS = 0; // control for application of ellipsoids

    SUMMARY = fiparam ( fparam1 ); // control for summary
    nsoilclass = fiparam ( fparam1 ); // number of soil classes
    nlayersl = fiparam ( fparam1 ); // total number of layers
    if ( LAYERS == 1 ) MAXLAYERS = fiparam ( fparam1 ); // maximum number of layers for one ellipsoid
    USESEG = fiparam ( fparam1 ); // whether to use segment library
    PF = fiparam ( fparam1 ); // whether to compute slope failure probability
    PFDEPTH = PF; // control for testing various depths


// Sampling of parameters (for slope failure probability)

    if ( PF == 1 ) {

        PFSAMPLING = fiparam ( fparam1 ); // sampling strategy for computing slope failure probability
        if ( INFMOD == 1 ) PFSAMPLING = 2; // forcing regular sampling for infinite slope stability model
        pf_nsteps_csoil = fiparam ( fparam1 ); // number of steps of c for computing slope failure probability 
        pf_nsteps_phi = fiparam ( fparam1 ); // number of steps of phi for computing slope failure probability
        pf_nsteps_depth = fiparam ( fparam1 ); // number of steps of depth for computing slope failure probability
        
        if ( pf_nsteps_depth < 2 || LAYERS == 1 || ( TRUNCATED == 0 && INFMOD == 0 ) ) {
            pf_nsteps_depth = 1; // setting minimum number of tested depths to 1
            PFDEPTH = 0; // setting control for testing various depths to negative, if necessary
        }


    // Pdf types: 0=rectangular, 1=gaussian, 2=lognormal, 3=exponential, 4=relation c-phi

        pdftype_csoil = fiparam ( fparam1 ); // type of probability density function for c
        pdftype_phi = fiparam ( fparam1 ); // type of probability density function for phi
        pdftype_depth = fiparam ( fparam1 ); // type of probability density function for depth
      
        if ( pdftype_csoil == 4 ) {
        
            regpar1 = fiparam( fparam1 ); // slope for regression
            regpar2 = fiparam( fparam1 ); // intercept for regression
        }
    }


// Defining further variables related to sampling of parameters

    else {
        pf_nsteps_phi = 0;
        pf_nsteps_csoil = 0;
        pf_nsteps_depth = 0;
    }
    pf_nsteps = pf_nsteps_phi * pf_nsteps_csoil * pf_nsteps_depth;

    if ( PF == 0 || PFDEPTH == 1 ) {

        if ( PF == 0 ) pf_ctr_max = 1;
        else pf_ctr_max = pf_nsteps_csoil * pf_nsteps_phi;
    }
    else pf_ctr_max = pf_nsteps;

    if ( PF == 1 && PFSAMPLING == 0 ) {

        pf_nsteps_phi = pf_nsteps;
        pf_nsteps_csoil = pf_nsteps;

        if ( PFDEPTH == 1 ) {
            pf_nsteps_depth = pf_nsteps;
            pf_ctr_max = 1;
        }
        else {
            pf_ctr_max = pf_nsteps;    
        }
    }


// Sampling of truncated depth (for slope failure probability)

    if ( PFDEPTH == 1 ) {

        pf_depth_mean = fdparam ( fparam1 ); // arithmetic mean of depth of layer 1
        pf_depth_stdev = fdparam ( fparam1 ); // standard deviation of depth of layer 1
        pf_depth_min = fdparam ( fparam1 ); // minimum of depth of layer 1
        pf_depth_max = fdparam ( fparam1 ); // maximum of depth of layer 1
        pf_depth_int = pf_depth_max - pf_depth_min; // interval between maximum and minimum of depth of layer 1
    }

    nsegrc = fiparam ( fparam1 ); // segment size
    nsegs = fiparam ( fparam1 ); // number of segments to be held in memory


// Reading parameter file 2

    if ( CLASSINF == 1 ) { // for soil class or infinite slope stability mode:

        nly = 0; // initializing maximum number of layers
        for ( i = 1; i <= nsoilclass; i++ ) {

            if ( fparam2 != NULL ) { if( fgets(sparam, 8, fparam2) == NULL ) return -1; } // number of layers
            nlayer[i]=strtod(sparam, strtodhlp); // converting to number
            if ( nlayer[i] > nly ) nly = nlayer[i]; // updating maximum number of layers
        }
        nlayersmax = nly + 2; // maximum number of layers per ellipsoid
    }
    else nlayersmax = 102; // maximum number of layers per ellipsoid


// Soil parameter arrays

    float **pgamma0 = 0, **pcsoil0 = 0, **pcsoils0 = 0, **pcsoilmin0 = 0, **pcsoilmax0 = 0, **pphi0 = 0, **pphis0 = 0, **pphimin0 = 0, **pphimax0 = 0, **ptheta0 = 0;
    float **pgamma = 0, **pcsoil = 0, **pcsoils = 0, **pcsoilmin = 0, **pcsoilmax = 0, **pphi = 0, **pphis = 0, **pphimin = 0, **pphimax = 0, **ptheta = 0, **param_depth0 = 0;

    pgamma0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
    pcsoil0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
    pphi0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
    ptheta0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
    pgamma = G_alloc_dmatrix(nsoilclass+1, nlayersmax);
    pcsoil = G_alloc_dmatrix(nsoilclass+1, nlayersmax);
    pphi = G_alloc_dmatrix(nsoilclass+1, nlayersmax);
    ptheta = G_alloc_dmatrix(nsoilclass+1, nlayersmax);
    param_depth0 = G_alloc_dmatrix(nsoilclass+1, nlayersmax);

    if ( PF == 1) { // for slope failure probability:

        pcsoils0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
        pcsoilmin0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
        pcsoilmax0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
        pphis0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
        pphimin0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
        pphimax0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
        pcsoils = G_alloc_dmatrix(nsoilclass+1, nlayersmax);
        pcsoilmin = G_alloc_dmatrix(nsoilclass+1, nlayersmax);
        pcsoilmax = G_alloc_dmatrix(nsoilclass+1, nlayersmax);
        pphis = G_alloc_dmatrix(nsoilclass+1, nlayersmax);
        pphimin = G_alloc_dmatrix(nsoilclass+1, nlayersmax);
        pphimax = G_alloc_dmatrix(nsoilclass+1, nlayersmax);
    }


// Reading parameter file 3

    if ( CLASSINF == 1 ) {  // for soil class or infinite slope stability mode:

        for ( i = 1; i <= nsoilclass; i++ ) { // loop over all soil classes:
            for ( j = 1; j <= nlayer[i]; j++ ) { // loop over all layers:

                pgamma[i][j] = fdparam ( fparam3 ); // soil weight
                pcsoil[i][j] = fdparam ( fparam3 ); // soil cohesion
                pphi[i][j] = fdparam ( fparam3 ) * PI / 180; // soil angle of internal friction
                ptheta[i][j] = fdparam ( fparam3 ) / 100; // soil water content

                if ( PF == 1 ) { // for slope failure probability:

                    pcsoils[i][j] = fdparam ( fparam3 ); // soil cohesion std. dev. 
                    pcsoilmin[i][j] = fdparam ( fparam3 ); // soil cohesion minimum 
                    pcsoilmax[i][j] = fdparam ( fparam3 ); // soil cohesion maximum 
                    pphis[i][j] = fdparam ( fparam3 ) * PI / 180; // soil angle of internal friction std. dev. 
                    pphimin[i][j] = fdparam ( fparam3 ) * PI / 180; // soil angle of internal friction minimum
                    pphimax[i][j] = fdparam ( fparam3 ) * PI / 180; // soil angle of internal friction maximum
                }

                if ( SEEPAGE == 0 ) ptheta[i][j]=0; // setting soil water content to zero if dry soil is assumed

                param_depth0[i][j] = fdparam ( fparam3 ); // soil depth
            }
        }
    }

    else { // for layer mode:

        for ( i = 1; i <= nsoilclass; i++ ) { // loop over all soil classes:

            pgamma0[i][0] = fdparam ( fparam2 ); // soil weight
            pcsoil0[i][0] = fdparam ( fparam2 ); // soil cohesion
            pphi0[i][0] = fdparam ( fparam2 ) * PI / 180; // soil angle of internal friction
            ptheta0[i][0] = fdparam ( fparam2 ) / 100; // soil water content

            if ( PF == 1 ) { // for slope failure probability:

                pcsoils0[i][0] = fdparam ( fparam2 ); // soil cohesion std. dev.
                pcsoilmin0[i][0] = fdparam ( fparam2 ); // soil cohesion minimum
                pcsoilmax0[i][0] = fdparam ( fparam2 ); // soil cohesion maximum 
                pphis0[i][0] = fdparam ( fparam2 ) * PI / 180; // soil angle of internal friction std. dev. 
                pphimin0[i][0] = fdparam ( fparam2 ) * PI / 180; // soil angle of internal friction minimum 
                pphimax0[i][0] = fdparam ( fparam2 ) * PI / 180; // soil angle of internal friction maximum 
            }

            if ( SEEPAGE == 0 ) ptheta0[i][0]=0;  // setting soil water content to zero if dry soil is assumed

            pgamma[i][0] = pgamma0[i][0]; // defining parameters for soil classes (to be applied where no layers are defined)
            pcsoil[i][0] = pcsoil0[i][0];
            pphi[i][0] = pphi0[i][0];
            ptheta[i][0] = ptheta0[i][0];

            if ( PF == 1 ) { // for slope failure probability:

                pcsoils[i][0] = pcsoils0[i][0];
                pcsoilmin[i][0] = pcsoilmin0[i][0];
                pcsoilmax[i][0] = pcsoilmax0[i][0];
                pphis[i][0] = pphis0[i][0];
                pphimin[i][0] = pphimin0[i][0];
                pphimax[i][0] = pphimax0[i][0];
             }
        }

        for ( j = 1; j <= nlayersl; j++ ) { // loop over all layers:

            pgamma0[0][j] = fdparam ( fparam3 ); // soil weight
            pcsoil0[0][j] = fdparam ( fparam3 ); // soil cohesion
            pphi0[0][j] = fdparam ( fparam3 ) * PI / 180; // soil angle of internal friction
            ptheta0[0][j] = fdparam ( fparam3 ) / 100; // soil water content

            if ( PF == 1 ) { // for slope failure probability:

                pcsoils0[0][j] = fdparam ( fparam3 );  // soil cohesion std. dev.
                pcsoilmin0[0][j] = fdparam ( fparam3 ); // soil cohesion minimum
                pcsoilmax0[0][j] = fdparam ( fparam3 ); // soil cohesion maximum
                pphis0[0][j] = fdparam ( fparam3 ) * PI / 180; // soil angle of internal friction std. dev.
                pphimin0[0][j] = fdparam ( fparam3 ) * PI / 180; // soil angle of internal friction minimum
                pphimax0[0][j] = fdparam ( fparam3 ) * PI / 180; // soil angle of internal friction maximum
            }

            if ( SEEPAGE == 0 ) ptheta0[0][j]=0; // setting soil water content to zero if dry soil is assumed
        }
    }

    if ( ELLIPS == 1 ) { // if ellipsoid mode is applied:

        LENGTHMAX = fdparam ( fparam3 ); // maximum length of landslide scarp (slope-parallel)
        LENGTHMIN = fdparam ( fparam3 ); // minimum length of landslide scarp (slope-parallel)
        WIDTHMAX = fdparam ( fparam3 ); // maximum width of landslide scarp
        WIDTHMIN = fdparam ( fparam3 ); // minimum width of landslide scarp
        DEPTHMAX = fdparam ( fparam3 ); // maximum depth of landslide scarp (vertical)
        DEPTHMIN = fdparam ( fparam3 ); // minimum depth of landslide scarp (vertical)
        ZBMAX = fdparam ( fparam3 ); // maximum offset of center of ellipsoid from terrain (relative to c, perpendicular to the slope)
        ZBMIN = fdparam ( fparam3 ); // minimum offset of center of ellipsoid from terrain (relative to c, perpendicular to the slope)
        RAT_GEOM_MAX = fdparam ( fparam3 ); // maximum length to depth ratio
        RAT_GEOM_MIN= fdparam ( fparam3 ); // minimum length to depth ratio
        CENTERX = fdparam ( fparam3 ); // x coordinate of center of landslide scarp (on ground)
        CENTERY = fdparam ( fparam3 ); // y coordinate of center of landslide scarp (on ground)
        ELL_ALPHADEF = fdparam ( fparam3 ); // aspect of ellipsoid
        ELL_BETADEF = fdparam ( fparam3 ); // slope of ellipsoid
        ELLDENS = fdparam ( fparam3 ); // ellipsoid density
    }


// Closing parameter files

    fclose (fparam1);
    fclose (fparam2);
    fclose (fparam3);


// Maximum extent of ellipsoid

    if ( ELLIPS == 1 ) {

        if ( LENGTHMAX > WIDTHMAX ) {

            me = (int)(( LENGTHMAX / pow( 1 - ZBMAX, 0.5 )) / csz ) + 8;
            ne = (int)(( LENGTHMAX / pow( 1 - ZBMAX, 0.5 )) / csz ) + 8;
        }

        else {

            me = (int)(( WIDTHMAX / pow( 1 - ZBMAX, 0.5 )) / csz ) + 8;
            ne = (int)(( WIDTHMAX / pow( 1 - ZBMAX, 0.5 )) / csz ) + 8;
        }
    }


// Segment files and variables for input parameters and preprocessing

    float param_elev, param_gwdepth, depthl, *gnewmark;

    char *msoilclass, *melev, *mpga, *mdistance, *mnewmarkref, *mgwdepth, *mdepth0, mdepth[nly+1][20], *mdepthl0, mdepthl[20], *mdzdx, *mdzdy, sxarr[20], syarr[20],
        slayerarr[20], sdeptharr[20];

    SEGMENT seg_depth[nly+1];


// Preparing segment files for input rasters and preprocessing

    if ( LAYERS == 1 ) { // for layer mode:

        Segment_open (&seg_ctrl_lyr, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(int), nsegs); // control for layers
        Segment_open (&seg_num_lyr, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(int), nsegs); // layer number
    }


// Defining names of input raster maps

    msoilclass="xsoilclass";
    melev="xelev";
    mgwdepth="xgwdepth";

    if ( PSEUDOSTATIC == 1 ) mpga="xpga";
    if ( NEWMARK == 1 ) mdistance="xdistance";
    if ( NEWMARKREF == 1 ) mnewmarkref="xnewmarkref";

    if ( CLASSINF == 1 && param_depth0[1][1] < 0 ) {
        mdepth0="xdepth";
        for ( j = 1; j <= nly; j++ ) sprintf(mdepth[j], "%s%i", mdepth0, j);
    }

    if ( LAYERS == 1 ) mdepthl0="xdepth";

    if ( LAYERS == 1 && SEEPAGE == 2 ) {
        mdzdx="xdzdxl";
        mdzdy="xdzdyl";
    }


// Opening input raster maps

    seg_soilclass = finrasti( msoilclass, UNDEFMIN );
    seg_elev = finrastd( melev, UNDEFMIN );
    seg_gwdepth = finrastd( mgwdepth, UNDEFMIN );

    if ( PSEUDOSTATIC == 1 ) seg_pga = finrastd( mpga, UNDEFMIN );
    if ( NEWMARK == 1 ) seg_distance = finrastd( mdistance, UNDEFMIN );
    if ( NEWMARKREF == 1 ) seg_newmarkref = finrastd( mnewmarkref, UNDEFMIN );

    Segment_open (&seg_depth[0], G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs);
    if ( CLASSINF == 1 ) {

        if ( param_depth0[1][1] >= 0 ) {
            for ( j = 0; j <= nly; j++ ) Segment_open (&seg_depth[j], G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs);
        }
        for ( j = 1; j <= nly; j++ ) {

            if ( param_depth0[1][1] < 0 ) seg_depth[j] = finrastd( mdepth[j], UNDEFMIN );
            else {
                for ( x = 0; x < m; x++ ) {
                    for ( y = 0; y < n; y++ ) {

                        Segment_get(&seg_soilclass, &pscs, x, y);
                        if ( pscs > 0 ) i = pscs;
                        else i = 0;
                        Segment_put(&seg_depth[j], &param_depth0[i][j], x, y);
                    }
                }
                Segment_flush(&seg_depth[j]);
            }
        }
    }

    if ( LAYERS == 1 && SEEPAGE == 2 ) {

        seg_dzdx = finrastd( mdzdx, UNDEFMIN );
        seg_dzdy = finrastd( mdzdy, UNDEFMIN );
    }


// Preparing table with layers

    if ( LAYERS == 1 ) {


    // Preparing layer files and variables

        sprintf(path, "%s%i/layers.txt", wkdir, xint); // preparing layer file
        flayers=fopen(path, "w"); // opening layer file

        sprintf(path, "%s%i/layersf.txt", wkdir, xint); // preparing surface layer file
        flayersf=fopen(path, "w"); // opening surface layer file

        fprintf( flayersf, "X\tY\tL\tD\n" ); // writing header for surface layers file

        nlinesl = 0; // initializing number of lines of layer file

        for ( x = 0; x < m; x++ ) { // loop over all cells:
            for ( y = 0; y < n; y++ ) {

                ctrl_lyr = 1; // initializing control for layers
                num_lyr = 0; // initializing number of layers

                Segment_put(&seg_ctrl_lyr, &ctrl_lyr, x, y);
                Segment_put(&seg_num_lyr, &num_lyr, x, y);
            }
        }
        Segment_flush(&seg_ctrl_lyr);
        Segment_flush(&seg_num_lyr);

        for ( i = 1; i <= nlayersl; i++ ) { // loop over all layers:


        // Opening raster map with depth of layer bottom

            sprintf(mdepthl, "%s%i", mdepthl0, i);
            seg_depthl = finrastd( mdepthl, UNDEFLYR );


        // For each raster cell, printing information to layers file

            for ( x = 0; x < m; x++ ) {
                for ( y = 0; y < n; y++ ) {

                    Segment_get(&seg_depthl, &depthl, x, y);

                    if ( depthl != UNDEFLYR ) {

                        fprintf ( flayers, "%i %i %i %.3f\n", x, y, i, depthl );
                        Segment_get(&seg_ctrl_lyr, &ctrl_lyr, x, y);

                        if ( ctrl_lyr == 1 ) {

                            fprintf( flayersf, "%.0f\t%.0f\t%.0f\t%.1f\n",
                                bd_west + csz * (float)y, bd_north - csz * (float)x, (float)i, depthl );
                            ctrl_lyr = 0;
                            Segment_put(&seg_ctrl_lyr, &ctrl_lyr, x, y);
                        }

                        Segment_get(&seg_num_lyr, &num_lyr, x, y);
                        num_lyr += 1;
                        Segment_put(&seg_num_lyr, &num_lyr, x, y);
                        nlinesl += 1;
                    }
                }
            }

            Segment_flush(&seg_ctrl_lyr);
            Segment_flush(&seg_num_lyr);
            Segment_release(&seg_depthl);

            printf("Preprocessing layer %i of %i ...\r", i, nlayersl);
            fflush(stdout);
        }

        kmax = 0;
        for ( x = 0; x < m; x++ ) {
            for ( y = 0; y < n; y++ ) {

                Segment_get(&seg_num_lyr, &num_lyr, x, y);
                if ( num_lyr > kmax ) kmax = num_lyr;
                Segment_put(&seg_num_lyr, &num_lyr, x, y);
            }
        }
        Segment_flush(&seg_num_lyr);


    // Closing segment and layer files

        Segment_release(&seg_ctrl_lyr);
        Segment_release(&seg_num_lyr);
        fclose (flayersf);
        fclose (flayers);

        if ( MAXLAYERS == 0 ) nly = 8 * kmax + 1; // maximum number of layers in local reference system
        else nly = MAXLAYERS + 1; 
    }

    if ( LAYERS == 0 ) nlinesl = 1; // setting number of lines of layer array to 1 for other modes
    int **ldata; // array for layer data
    ldata = G_alloc_imatrix(nlinesl, 4);


// Writing content of layer file to array

    if ( LAYERS == 1 ) {

        sprintf(path, "%s%i/layers.txt", wkdir, xint); // searching for layers file
        flayers=fopen(path, "r"); // opening layers file

        for ( k = 0; k < nlinesl; k++ ) { // starting loop over all lines of layers file

            if( fgets(slayers, 100, flayers) == NULL ) return -1; // reading line from layers file
            sscanf(slayers, "%s %s %s %s", sxarr, syarr, slayerarr, sdeptharr); // decomposing string

            ldata[k][0]=strtod(sxarr, strtodhlp); // converting to number and writing to array
            ldata[k][1]=strtod(syarr, strtodhlp);
            ldata[k][2]=strtod(slayerarr, strtodhlp);
            ldata[k][3]=(int)(0.5 + strtod(sdeptharr, strtodhlp)*CHANGEDM);
        }
        fclose (flayers);

        printf("Preprocessing layers ... completed.         \n");
    }

    if ( PFDEPTH == 1 ) jmax = pf_nsteps_depth;
    else jmax = nly;


// Various working and output variables

// ATTENTION! jctrl!

    int **epd[jmax+2], **epe = 0, **epga = 0, **dgw = 0, **eparam_soilclass = 0, **ezslip = 0, epds, epd1s, epd2s, epd_nhs[9], epd_nhbs[9], ldata_trans[jmax+1],
        epes, epgas, epe_nhs[9], dgw_nhs[9], ezslips, ezslip_nhs[9], dgws, dslip_test[jmax+1], enh[9], ctr_mcsim, randmax, ell_centerx, ell_centery,
        xp, xm, yp, ym, ex, ey, epsc, *epsc_nh, interval_ctr, elldens_interval = 0, area_eff[jmax+1],
        validity[jmax+1], num_elldens_interval, ctr_topo, ctrl_content, ell_centerxd, ell_centeryd, lyr_crit, xnh[9], ynh[9], q; // jctrl;

    float **r_cum, **t_cum, **fos_surf, **pnewmark_surf = 0, **dnewmark_surf = 0, ienh[9];

        r_cum = G_alloc_dmatrix(jmax+1, pf_ctr_max+1);
        t_cum = G_alloc_dmatrix(jmax+1, pf_ctr_max+1);
        fos_surf = G_alloc_dmatrix(jmax+1, pf_ctr_max+1);
        if ( NEWMARK == 1 ) {
            dnewmark_surf = G_alloc_dmatrix(jmax+1, pf_ctr_max+1);
            pnewmark_surf = G_alloc_dmatrix(jmax+1, pf_ctr_max+1);
        }

    if ( ELLIPS == 1 && USESEG == 0 ) {

        epe = G_alloc_imatrix(me, ne);
        if ( PSEUDOSTATIC == 1 ) epga = G_alloc_imatrix(me, ne);
        dgw = G_alloc_imatrix(me, ne);
        ezslip = G_alloc_imatrix(me, ne);
        eparam_soilclass = G_alloc_imatrix(me, ne);

        for ( j=0; j <= jmax+1; j++) {
            epd[j] = G_alloc_imatrix(me, ne);
        }
    }

    float ell_a, ell_b, ell_c, ell_max, alpha, beta, zell_0, *efnh, csoil = 0, phi = 0,
        xloc1_ell, xloc2_ell, yloc1_ell, yloc2_ell, zloc2_ell, beta_c, dz_dx, dz_dy, beta_xz,
        beta_yz, beta_m = 0, area, pfail, pfaill[jmax+1], fos_mean, fos_meanl[jmax+1], fos_stdev, fos_stdevl[jmax+1], pfail0, fos_stdev0,
        fosmin = 0, fosdeep = 0, fos_min, fos_deep, dnewmark = 0, d_newmark, dnewmark_stdev, dnewmark_mean, dnewmark_meanl[jmax+1], dnewmark_stdev0, dnewmark_stdevl[jmax+1], pnewmark, pnewmarkl[jmax+1], pnewmark0, 
        depthfail_min, depthfail_deep, depthfail_crit, amin, bmin, zbmin, 
	xdiff, ydiff, delta_xloc2, qeq_a, qeq_b, qeq_c, weight_sub, u_appr, ell_zb, alpha_c, fit, 
        slope, gslope, beta_w, alpha_w, depth_control, ise[9], isdn[9], isd[100], isg[9], ignh[9], isd11, ell_zb_rat,
        f_s, f_sh, f_sv, f_sch, f_smh, f_sc, f_sm, beta_fsc, beta_fsm, n_s, t_s, dz1, dz2, depth_layer, arr_depthr, pccomp,
        r_cumcell = 0, t_cumcell = 0, depthslipmin = 0, ell_centerx_metric, ell_centery_metric,
        ell_length, ell_width, ell_area, ell_depth, fos_critabs, ell_centerx_crit = 0, yell0_crit = 0, ella_crit = 0, ellb_crit = 0, ellc_crit = 0,
        zb_crit = 0, fit_crit = 0, alpha_crit = 0, beta_crit = 0, ell_width_crit = 0, ell_length_crit = 0, ell_area_crit = 0, ell_depth_crit = 0, corr_topo,
        s2_dxtest, s2_dytest, s2_d, s2_dcum[501], s2_data[501][nly+1][8], pga, distance, seisforce, newmarkref;

    int s2_x0, s2_y0, s2_x[501], s2_y[501], s2_ctrl = 0, s2_ctrl2, s2_dprex, s2_dprey, jcrit, jcrit2 = UNDEFMIN, output_add,
        cumfos_critmax = 0, cumfos_crit, cumsurf;

    efnh = (float*) malloc( sizeof(float) * 9 );
    epsc_nh = (int*) malloc( sizeof(int) * 9 );


// Variables for sampling of parameters (slope failure probability)

    int pf_ctr, pf_ctr0, pf_ctr_phi = 0, pf_ctr_csoil = 0, pf_depth[pf_nsteps_depth+1], pf_it,
        **pf_phi0[pf_nsteps_phi+1], **pf_csoil0[pf_nsteps_csoil+1], **pf_phi[pf_nsteps_phi+1], **pf_csoil[pf_nsteps_csoil+1];

    float pf_pdf0, pf_depth0, pf_pdf_test, pf_phi_test, pf_csoil_test, pf_depth_test, pf_pdf_phi, pf_pdf_csoil, pf_pdf_depth, 
        pf_cdf_max_depth = 0, pf_cdf_max_phi, pf_cdf_max_csoil, pf_phi_mean, pf_phi_stdev, pf_phi_min, pf_phi_max, pf_phi_int, pf_csoil_mean,
        pf_csoil_stdev, pf_csoil_min, pf_csoil_max, pf_csoil_int, pf_pdf_peak_depth = 0, pf_pdf_peak_phi00, pf_pdf_peak_csoil00,
        pf_cdf_depth, pf_cdf_phi, pf_cdf_csoil, pf_phi00, pf_csoil00, pf_cdf_pre, pf_weight_depth[pf_nsteps_depth+1], pf_weight_fphi,
        pf_weight_fcsoil, pf_weight_geo[pf_ctr_max+1], **pf_pdf_peak_phi0 = 0, **pf_pdf_peak_csoil0 = 0, **pf_pdf_peak_phi = 0, **pf_pdf_peak_csoil = 0;

    if ( PF == 1) {

        pf_pdf_peak_phi0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
        pf_pdf_peak_phi = G_alloc_dmatrix(nsoilclass+1, nlayersmax);
        pf_pdf_peak_csoil0 = G_alloc_dmatrix(nsoilclass+1, nlayersl+1);
        pf_pdf_peak_csoil = G_alloc_dmatrix(nsoilclass+1, nlayersmax);

        for ( pf_ctr = 0; pf_ctr <= pf_nsteps_phi; pf_ctr++ ) {
            pf_phi0[pf_ctr] = G_alloc_imatrix(nsoilclass+1, nlayersl+1);
            pf_phi[pf_ctr] = G_alloc_imatrix(nsoilclass+1, nlayersmax);
        }

        for ( pf_ctr = 0; pf_ctr <= pf_nsteps_csoil; pf_ctr++ ) {
            pf_csoil0[pf_ctr] = G_alloc_imatrix(nsoilclass+1, nlayersl+1);
            pf_csoil[pf_ctr] = G_alloc_imatrix(nsoilclass+1, nlayersmax);
        }
    }


// Variables for names of output raster maps

    char *mfos_min, *mcumfos_crit, *mcumsurf, *mpfail, *mfos_stdev, *mfos_deep, *mfos_pseudo, *mdnewmark, *mdepthfail_min,
        *mdepthfail_deep, *mdepthfail_crit, *mlyr_crit, *mamin, *mbmin, *mzbmin;

// Data files

    FILE *fsummary1, *fprofile, *fsummary2, *fchange;


// Segment variables for processing and output

    SEGMENT seg_epd[jmax+2], seg_epe, seg_epga, seg_dgw, seg_ezslip, seg_esoilclass, seg_fos_min, seg_cumfos_crit, seg_cumsurf, seg_pfail, seg_fos_stdev,
        seg_fos_deep, seg_fos_pseudo, seg_d_newmark, seg_pnewmark, seg_d_newmark_stdev,
        seg_depthfail_min, seg_depthfail_deep, seg_depthfail_crit, seg_lyr_crit, seg_amin, seg_bmin, seg_zbmin;


// Preparing segment files for processing and output

    if ( ELLIPS == 1 && USESEG == 1 ) { 

        for ( j = 0; j <= jmax+1; j++ ) Segment_open (&seg_epd[j], G_tempfile(), me, ne, nsegrc, nsegrc, sizeof(int), nsegs);
            // depth of layer bottom
        Segment_open (&seg_epe, G_tempfile(), me, ne, nsegrc, nsegrc, sizeof(int), nsegs); // elevation
        if ( PSEUDOSTATIC == 1 ) Segment_open (&seg_epga, G_tempfile(), me, ne, nsegrc, nsegrc, sizeof(float), nsegs); // peak ground acceleration
        Segment_open (&seg_dgw, G_tempfile(), me, ne, nsegrc, nsegrc, sizeof(int), nsegs); // depth of groundwater table
        Segment_open (&seg_ezslip, G_tempfile(), me, ne, nsegrc, nsegrc, sizeof(int), nsegs); // slip surface elevation
        Segment_open (&seg_esoilclass, G_tempfile(), me, ne, nsegrc, nsegrc, sizeof(int), nsegs); // soil class 
    }

    Segment_open (&seg_fos_min, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs); // minimum factor of safety
    if ( PSEUDOSTATIC == 1 ) Segment_open (&seg_fos_pseudo, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs); // factor of safety (pseudostatic seismic slope stability analysis)
    if ( NEWMARK == 1 ) {
        Segment_open (&seg_d_newmark, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs); // Newmark displacement
        Segment_open (&seg_pnewmark, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs);
        Segment_open (&seg_d_newmark_stdev, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs);
    }
    Segment_open (&seg_pfail, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs); // slope failure probability
    Segment_open (&seg_fos_stdev, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs); // standard deviation of factor of safety
    Segment_open (&seg_depthfail_deep, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs); // depth of deepest slip surface with factor of safety lower than 1
    Segment_open (&seg_fos_deep, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs); // factor of safety for deepest slip surface with factor of safety lower than 1

    if ( ADDITIONAL == 1 ) {

        Segment_open (&seg_cumfos_crit, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(int), nsegs);
            // number of slip surfaces with factor of safety lower than 1
        Segment_open (&seg_cumsurf, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(int), nsegs); // total number of slip surfaces
        Segment_open (&seg_depthfail_min, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs);
            // depth of slip surface with minimum factor of safety
        Segment_open (&seg_depthfail_crit, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs);
            // depth of most critical slip surface
        Segment_open (&seg_lyr_crit, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(int), nsegs);
            // id of slip surface with lowest factor of safety ( layer mode)
        Segment_open (&seg_amin, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs);
            // a axis of slip surface with minimum factor of safety
        Segment_open (&seg_bmin, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs);
            // b axis of slip surface with minimum factor of safety
        Segment_open (&seg_zbmin, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs);
            // relative zb of slip surface with minimum factor of safety
    }

    if ( ELLDENS == 0 ) Segment_open (&seg_depthfail_min, G_tempfile(), m, n, nsegrc, nsegrc, sizeof(float), nsegs);
        // depth of slip surface with minimum factor of safety


// Defining names of output raster maps

    mfos_min = "r_fos";
    if ( PSEUDOSTATIC == 1 ) mfos_pseudo = "r_fos_pseudo";
    if ( NEWMARK == 1 ) mdnewmark = "r_dnewmark";
    mpfail = "r_pfail";
    mfos_stdev = "r_fos_stdev";
    mdepthfail_deep = "r_ddepth";
    if ( ADDITIONAL == 1 ) {
        mfos_deep = "r_dfos";
        mcumfos_crit = "r_cumfos_crit";
        mcumsurf = "r_cumsurf";
        mdepthfail_min = "r_depth";
        mdepthfail_crit = "r_depthfail_crit";
        mlyr_crit = "r_layer";
        mamin = "r_lcrit";
        mbmin = "r_wcrit";
        mzbmin = "r_zbcrit";
    }
    if ( ELLDENS == 0 ) mdepthfail_min = "r_depth"; 


// Opening and preparing summary and change files

    time_start = clock(); // start time

    if ( SUMMARY == 1 ) { // if mode includes summary:

        sprintf(path, "%s%i/summary1%i.txt", wkdir, xint, xint); // creating first summary file
        fsummary1=fopen(path, "w"); // opening first summary file

        fprintf ( fsummary1, "id\tell\tx0\ty0\ta\tb\tc\tzb\tfit\talpha\tbeta\tL\tW\tA\tD\tfos" ); // writing header to first summary file
        if ( TRUNCATED == 1 ) for ( j = 1; j <= nly; j++ ) fprintf ( fsummary1, "\tfos%i", j );
        fprintf ( fsummary1, "\ttime\n" );
    }
    if ( ELLDENS == 0 ) { // if only one ellipsoid is to be tested:
        sprintf(path, "%s%i/summary2%i.txt", wkdir, xint, xint); // creating profile file
        fprofile=fopen(path, "w"); // opening profile file
    }
    else { // if more than one ellipsoid is to be tested:
        sprintf(path, "%s%i/evolution%d.txt", wkdir, xint, xint); // creating evolution file
        fchange=fopen(path, "w"); // writing header to change file
    }
    sprintf(path, "%s%i/summary3%i.txt", wkdir, xint, xint); // creating second summary file
    fsummary2=fopen(path, "w"); // opening second summary file


// Initializing and preparing working and output variables

    for ( x = 0; x < m; x++ ) { // loop over all pixels:
        for ( y = 0; y < n; y++ ) {

            fos_min = UNDEFPLUS;
            if ( NEWMARK == 1 ) d_newmark = UNDEFMIN;
            Segment_put(&seg_fos_min, &fos_min, x, y);
            if ( NEWMARK == 1 ) Segment_put(&seg_d_newmark, &d_newmark, x, y);

            if ( PF == 1 ) { // for slope failure probability:
                pfail0=0;
                fos_stdev0=UNDEFMIN;
                Segment_put(&seg_pfail, &pfail0, x, y);
                Segment_put(&seg_fos_stdev, &fos_stdev0, x, y);

                if ( NEWMARK == 1 ) {

                    pnewmark0=0;
                    dnewmark_stdev0=UNDEFMIN;
                    Segment_put(&seg_pnewmark, &pnewmark0, x, y);
                    Segment_put(&seg_d_newmark_stdev, &dnewmark_stdev0, x, y);
                }
            }

            depthfail_deep = UNDEFMIN;
            fos_deep = UNDEFPLUS;
            Segment_put(&seg_depthfail_deep, &depthfail_deep, x, y);
            Segment_put(&seg_fos_deep, &fos_deep, x, y);

            if ( ADDITIONAL == 1 ) { // for additional output:

                cumfos_crit = 0;
                cumsurf = 0;
                depthfail_min = UNDEFMIN;
                depthfail_crit = UNDEFMIN;
                amin = UNDEFMIN;
                bmin = UNDEFMIN;
                zbmin = UNDEFMIN;

                Segment_put(&seg_cumfos_crit, &cumfos_crit, x, y);
                Segment_put(&seg_cumsurf, &cumsurf, x, y);
                Segment_put(&seg_depthfail_min, &depthfail_min, x, y);
                Segment_put(&seg_depthfail_crit, &depthfail_crit, x, y);
                Segment_put(&seg_amin, &amin, x, y);
                Segment_put(&seg_bmin, &bmin, x, y);
                Segment_put(&seg_zbmin, &zbmin, x, y);

                if ( LAYERS == 1 ) { // for layer mode:
                    lyr_crit = UNDEFMIN;
                    Segment_put(&seg_lyr_crit, &lyr_crit, x, y);
                }
            }

            if ( ELLDENS == 0 ) {
                depthfail_min = UNDEFMIN;
                Segment_put(&seg_depthfail_min, &depthfail_min, x, y);
            }

            if ( INFMOD == 1 ) { // for infinite slope stability mode
                isd[0] = 0;
                Segment_put(&seg_depth[0], &isd[0], x, y);
            }
        }
    }
    Segment_flush(&seg_fos_min);
    if ( NEWMARK == 1 ) Segment_flush(&seg_d_newmark);

    if ( PF == 1 ) {
        Segment_flush(&seg_pfail);
        Segment_flush(&seg_fos_stdev);
        if ( NEWMARK == 1 ) Segment_flush(&seg_pnewmark);
        if ( NEWMARK == 1 ) Segment_flush(&seg_d_newmark_stdev);
    }

    Segment_flush(&seg_depthfail_deep);
    Segment_flush(&seg_fos_deep);

    if ( ADDITIONAL == 1 ) {

        Segment_flush(&seg_cumfos_crit);
        Segment_flush(&seg_cumsurf);
        Segment_flush(&seg_depthfail_min);
        Segment_flush(&seg_depthfail_crit);
        Segment_flush(&seg_amin);
        Segment_flush(&seg_bmin);
        Segment_flush(&seg_zbmin);
        if ( LAYERS == 1 ) Segment_flush(&seg_lyr_crit);
    }
    if ( ELLDENS == 0 ) Segment_flush(&seg_depthfail_min); 
    if ( INFMOD == 1 ) Segment_flush(&seg_depth[0]);


// Defining sample of parameters needed for slope failure probability

    if ( PF == 1 ) { // for slope failure probability:

        if ( PFDEPTH == 1 ) { // for sampling of depth values:

            pf_depth[0] = 0;
            my_depth = log ( pow ( pf_depth_mean, 2 ) / pow ( pow ( pf_depth_mean, 2 ) + pow ( pf_depth_stdev, 2 ) , 0.5 ) );
            sigma_depth = pow ( log( 1 + pow ( pf_depth_stdev, 2 ) / pow ( pf_depth_mean, 2 ) ) , 0.5 );


        // Peak of probability density for truncated depth

            if ( PFSAMPLING <= 1 ) { // for random sampling:

                if ( pdftype_depth == 0 ) pf_pdf_peak_depth = 1; // rectangular distribution
                else if ( pdftype_depth == 1 ) pf_pdf_peak_depth = 1 / ( pf_depth_stdev * pow ( 2 * PI, 0.5 )); // normal distribution
                else if ( pdftype_depth == 2 ) { // log-normal distribution
                    mode = exp ( my_depth - pow ( sigma_depth, 2 ));
                    pf_pdf_peak_depth = 1 / ( mode * sigma_depth * pow( 2 * PI, 0.5 ))
                        * exp ( -( pow( log( mode ) - my_depth, 2 ) / ( 2 * pow( sigma_depth , 2 ))));
                }
                else pf_pdf_peak_depth = 1 / pf_depth_mean; // exponential distribution
            }


        // Sampling of truncated depth values

            else if ( PFSAMPLING == 2 ) { // for regular sampling:

                for ( pf_it = 0; pf_it <= 1; pf_it++ ) {

                    pf_cdf_depth = 0;
                    pf_cdf_pre = 0;
                    pf_depth_test = pf_depth_min;
                    pf_ctr0 = 0;
                    pf_ctr = 0;
                    while ( pf_depth_test <= pf_depth_max || ( pf_it == 1 && pf_ctr < pf_nsteps_depth ) ) {

                        if ( pdftype_depth == 0 ) pf_pdf_depth = 1 / pf_depth_int; // rectangular distribution
                        else if ( pdftype_depth == 1 ) pf_pdf_depth = 1 / ( pf_depth_stdev * pow ( 2 * PI, 0.5 ) )
                            * exp ( -pow( pf_depth_test - pf_depth_mean, 2 ) / ( 2 * pow ( pf_depth_stdev, 2 ))); // normal distribution
                        else if ( pdftype_depth == 2 ) // log-normal distribution
                            pf_pdf_depth = 1 / ( pf_depth_test * sigma_depth * pow( 2 * PI, 0.5 ))
                                * exp ( -( pow( log( pf_depth_test ) - my_depth, 2 ) / ( 2 * pow( sigma_depth , 2 ))));
                        else pf_pdf_depth = ( 1 / pf_depth_mean * exp ( -( 1 / pf_depth_mean ) * pf_depth_test ));
                            // exponential distribution
                        pf_cdf_depth += pf_pdf_depth * pf_depth_int / 1000;

                        if ( pf_it == 1 && pf_cdf_depth >= pf_cdf_max_depth * (float)pf_ctr / ( (float)pf_nsteps_depth - 1 ) ) {

                            pf_ctr += 1;
                            pf_depth[pf_ctr] = (int)( CHANGEDM * pf_depth_test + 0.5 );
                        }

                        if ( pf_it == 1
                            && pf_cdf_depth >= pf_cdf_max_depth * ( ( (float)pf_ctr0 ) - 0.5 ) / ( (float)pf_nsteps_depth - 1 ) ) {

                            pf_ctr0 += 1;
                            pf_weight_depth[pf_ctr] = ( pf_cdf_depth - pf_cdf_pre ) / pf_cdf_max_depth;
                            pf_cdf_pre = pf_cdf_depth;
                        }

                        pf_depth_test += pf_depth_int / 1000;
                    }
                    if ( pf_it == 0 ) pf_cdf_max_depth = pf_cdf_depth; // maximum cumulative density
                }
                pf_weight_depth[pf_nsteps_depth] = ( pf_cdf_depth - pf_cdf_pre ) / pf_cdf_max_depth;
            }
        }

        for ( i=0; i<=nsoilclass; i++ ) { // loop over all soil classes and layers
            for ( j=0; j<=nlayersl; j++ ) {

                if ( ( CLASSINF == 1 && i > 0 && j > 0 ) || ( LAYERS == 1 && ( i > 0 || j > 0 ) && ( i == 0 || j == 0 ) ) ) {
                    // for all valid combinations:

                    if ( CLASSINF == 1 ) { // for soil class modes:
 
                        pf_phi_mean = pphi[i][j]; // statistics of internal friction angle (phi)
                        pf_phi_stdev = pphis[i][j];
                        pf_phi_min = pphimin[i][j];
                        pf_phi_max = pphimax[i][j];

                        pf_csoil_mean = pcsoil[i][j]; // statistics of cohesion (csoil)
                        pf_csoil_stdev = pcsoils[i][j];
                        pf_csoil_min = pcsoilmin[i][j];
                        pf_csoil_max = pcsoilmax[i][j];
                    }

                    else if ( LAYERS == 1 ) { // for layer mode:
 
                        pf_phi_mean = pphi0[i][j]; // statistics of internal friction angle (phi)
                        pf_phi_stdev = pphis0[i][j];
                        pf_phi_min = pphimin0[i][j];
                        pf_phi_max = pphimax0[i][j];

                        pf_csoil_mean = pcsoil0[i][j]; // statistics of cohesion (csoil)
                        pf_csoil_stdev = pcsoils0[i][j];
                        pf_csoil_min = pcsoilmin0[i][j];
                        pf_csoil_max = pcsoilmax0[i][j];
                    }

                    pf_phi_int = pf_phi_max - pf_phi_min;
                    pf_csoil_int = pf_csoil_max - pf_csoil_min;

                    my_phi = log ( pow ( pf_phi_mean, 2 ) / pow ( pow ( pf_phi_mean, 2 ) + pow ( pf_phi_stdev, 2 ) , 0.5 ) );
                    sigma_phi = pow ( log( 1 + pow ( pf_phi_stdev, 2 ) / pow ( pf_phi_mean, 2 ) ) , 0.5 );
                    my_csoil = log ( pow ( pf_csoil_mean, 2 ) / pow ( pow ( pf_csoil_mean, 2 ) + pow ( pf_csoil_stdev, 2 ) , 0.5 ) );
                    sigma_csoil = pow ( log( 1 + pow ( pf_csoil_stdev, 2 ) / pow ( pf_csoil_mean, 2 ) ) , 0.5 );

                    if ( PFSAMPLING <= 1 ) { // preparation for random sampling:


                    // Peak of probability density for internal friction angle (phi)

                        if ( pdftype_phi == 0 ) pf_pdf_peak_phi00 = 1;  // rectangular distribution
                        else if ( pdftype_phi == 1 ) pf_pdf_peak_phi00 = 1 / ( pf_phi_stdev * pow ( 2 * PI, 0.5 )); // normal distribution
                        else if ( pdftype_phi == 2 ) { // log-normal distribution
                            mode = exp ( my_phi - pow ( sigma_phi, 2 ));
                            pf_pdf_peak_phi00 = 1 / ( mode * sigma_phi * pow( 2 * PI, 0.5 ))
                                * exp ( -( pow( log( mode ) - my_phi, 2 ) / ( 2 * pow( sigma_phi , 2 ))));
                        }
                        else pf_pdf_peak_phi00 = 1 / pf_phi_mean; // exponential distribution

                        if ( CLASSINF == 1 ) pf_pdf_peak_phi[i][j] = pf_pdf_peak_phi00;
                        else pf_pdf_peak_phi0[i][j] = pf_pdf_peak_phi00;


                    // Peak of probability density for cohesion (csoil)

                        if ( pdftype_csoil == 0 ) pf_pdf_peak_csoil00 = 1;  // rectangular distribution
                        else if ( pdftype_csoil == 1 || pdftype_csoil == 4 )
                            pf_pdf_peak_csoil00 = 1 / ( pf_csoil_stdev * pow ( 2 * PI, 0.5 ));
                            // normal distribution
                        else if ( pdftype_csoil == 2 ) { // log-normal distribution
                            mode = exp ( my_csoil - pow ( sigma_csoil, 2 ));
                            pf_pdf_peak_csoil00 = 1 / ( mode * sigma_csoil * pow( 2 * PI, 0.5 ))
                                * exp ( -( pow( log( mode ) - my_csoil, 2 ) / ( 2 * pow( sigma_csoil , 2 ))));
                        }
                        else pf_pdf_peak_csoil00 = 1 / pf_csoil_mean; // exponential distribution

                        if ( CLASSINF == 1 ) pf_pdf_peak_csoil[i][j] = pf_pdf_peak_csoil00;
                        else pf_pdf_peak_csoil0[i][j] = pf_pdf_peak_csoil00;
                    }

                    else if ( PFSAMPLING == 2 ) { // for regular sampling:


                    // Sampling of internal friction angle (phi) values

                        for ( pf_it = 0; pf_it <= 1; pf_it++ ) {

                            pf_cdf_phi = 0;
                            pf_phi_test = pf_phi_min;
                            pf_ctr = 0;
                            while ( pf_phi_test <= pf_phi_max || ( pf_it == 1 && pf_ctr < pf_nsteps_phi ) ) {

                                if ( pdftype_phi == 0 ) pf_pdf_phi = 1 / pf_phi_int; // rectangular distribution
                                else if ( pdftype_phi == 1 ) pf_pdf_phi = 1 / ( pf_phi_stdev * pow ( 2 * PI, 0.5 ) )
                                    * exp ( -pow( pf_phi_test - pf_phi_mean, 2 ) / ( 2 * pow ( pf_phi_stdev, 2 ))); // normal distribution
                                else if ( pdftype_phi == 2 ) // log-normal distribution
                                    pf_pdf_phi = 1 / ( pf_phi_test * sigma_phi * pow( 2 * PI, 0.5 ))
                                        * exp ( -( pow( log( pf_phi_test ) - my_phi, 2 ) / ( 2 * pow( sigma_phi , 2 ))));
                                else pf_pdf_phi = ( 1 / pf_phi_mean * exp ( -( 1 / pf_phi_mean ) * pf_phi_test ));
                                    // exponential distribution
                                pf_cdf_phi += pf_pdf_phi * pf_phi_int / 1000;

                                if ( pf_it == 1 && pf_cdf_phi >= pf_cdf_max_phi * (float)pf_ctr / ( (float)pf_nsteps_phi - 1 ) ) {

                                    pf_ctr += 1;
                                    if ( CLASSINF == 1 || j == 0 ) pf_phi[pf_ctr][i][j] = (int)( 10000 * pf_phi_test + 0.5 );
                                    else pf_phi0[pf_ctr][i][j] = (int)( 10000 * pf_phi_test + 0.5 );
                                }
                                pf_phi_test += pf_phi_int / 1000;
                            }
                            if ( pf_it == 0 ) pf_cdf_max_phi = pf_cdf_phi; // maximum cumulative density
                        }


                    // Sampling of cohesion (csoil) values

                        for ( pf_it = 0; pf_it <= 1; pf_it++ ) {

                            pf_cdf_csoil = 0;
                            pf_csoil_test = pf_csoil_min;
                            pf_ctr = 0;
                            while ( pf_csoil_test <= pf_csoil_max || ( pf_it == 1 && pf_ctr < pf_nsteps_csoil ) ) {

                                if ( pdftype_csoil == 0 ) pf_pdf_csoil = 1 / pf_csoil_int; // rectangular distribution
                                else if ( pdftype_csoil == 1  || pdftype_csoil == 4 )
                                    pf_pdf_csoil = 1 / ( pf_csoil_stdev * pow ( 2 * PI, 0.5 ) )
                                    * exp ( -pow( pf_csoil_test - pf_csoil_mean, 2 ) / ( 2 * pow ( pf_csoil_stdev, 2 )));
                                    // normal distribution
                                else if ( pdftype_csoil == 2 ) // log-normal distribution
                                    pf_pdf_csoil = 1 / ( pf_csoil_test * sigma_csoil * pow( 2 * PI, 0.5 ))
                                        * exp ( -( pow( log( pf_csoil_test ) - my_csoil, 2 ) / ( 2 * pow( sigma_csoil , 2 ))));
                                else pf_pdf_csoil = ( 1 / pf_csoil_mean * exp ( -( 1 / pf_csoil_mean ) * pf_csoil_test ));
                                    // exponential distribution
                                pf_cdf_csoil += pf_pdf_csoil * pf_csoil_int / 1000;

                                if ( pf_it == 1 && pf_cdf_csoil >= pf_cdf_max_csoil * (float)pf_ctr / ( (float)pf_nsteps_csoil - 1 ) ) {

                                    pf_ctr += 1;
                                    if ( CLASSINF == 1 || j == 0 ) pf_csoil[pf_ctr][i][j] = (int)( pf_csoil_test + 0.5 );
                                    else pf_csoil0[pf_ctr][i][j] = (int)( pf_csoil_test + 0.5 );
                                }
                                pf_csoil_test += pf_csoil_int / 1000;
                            }
                            if ( pf_it == 0 ) pf_cdf_max_csoil = pf_cdf_csoil; // maximum cumulative density
                        }
                    }
                }
            }
        }


    // Weighting combinations of sampled geotechnical parameters

        if ( PFSAMPLING == 2 ) { // for regular sampling:

            pf_ctr_phi = 1; // resetting counter for sampling of phi
            pf_ctr_csoil = 1; // resetting counter for sampling of csoil

            for ( pf_ctr = 1; pf_ctr <= pf_ctr_max; pf_ctr++ ) {
                // starting loop over samples of phi and csoil:

                if ( pf_ctr_phi == 1 || pf_ctr_phi == pf_nsteps_phi ) pf_weight_fphi = 0.5;
                else pf_weight_fphi = 1; // relative weight of phi

                if ( pf_ctr_csoil == 1 || pf_ctr_csoil == pf_nsteps_csoil ) pf_weight_fcsoil = 0.5;
                else pf_weight_fcsoil = 1; // relative weight of csoil

                pf_weight_geo[pf_ctr] = pf_weight_fphi * pf_weight_fcsoil
                    / ( ( (float)pf_nsteps_phi - 1 ) * ( (float)pf_nsteps_csoil - 1 ) );
                    // combined weight of sampled parameters

                pf_ctr_phi += 1; // updating counters of phi and csoil
                if ( pf_ctr_phi > pf_nsteps_phi ) {
                    pf_ctr_phi = 1;
                    pf_ctr_csoil += 1;
                }
            }
        }
    }

    if ( ELLIPS == 1 ) { // if ellipsoid-based mode applies:


    // Initializing data for profile file

        if ( ELLDENS == 0 ) { // if only one ellipsoid is to be tested:

            for ( s2_ctrl2 = 0; s2_ctrl2 <= 500; s2_ctrl2++ ) {
                for ( j=0; j<= nly; j++ ) {

                    s2_data[s2_ctrl2][j][0] = 0;
                    s2_data[s2_ctrl2][j][1] = 0;
                    s2_data[s2_ctrl2][j][2] = 0;
                    s2_data[s2_ctrl2][j][3] = 0;
                    s2_data[s2_ctrl2][j][4] = 0;
                    s2_data[s2_ctrl2][j][5] = 0;
                    s2_data[s2_ctrl2][j][6] = 0;
                    s2_data[s2_ctrl2][j][7] = 0;
                }
            }
        }

        fos_critabs = UNDEFPLUS; // initializing absolute minimum factor of safety
        num_cellval = 0; // initializing number of cells with valid soil class
        lbtm_max = 0; // initializing maximum number of layers in local reference system
        ctr_mcsim = 0; // initializing counter for number of tested ellipsoids
        ctr_valid = 0; // initializing counter for number of valid ellipsoids
        num_elldens_interval=0; // initializing counter for ellipsoid density
        jcrit = UNDEFMIN; // initializing control for critical layer
        if ( ELLDENS == 0 ) depthslipmin = UNDEFPLUS; // lowest slip surface for profile
        
        if ( PF == 1 || PFDEPTH == 1 ) { iariasmin = 1; iariasmax = 101; }
        else { iariasmin = 1; iariasmax = 2; }


    // Topographic correction and number of valid pixels

        corr_topo = 0; // initializing topographic correction factor
        ctr_topo = 0; // initializing counter for number of pixels

        for ( x = 1; x < m-1; x++ ) { // loop over all pixels:
            for ( y = 1; y < n-1; y++ ) {

                ctrl_valtopo = 1; // resetting validity of pixel

                Segment_get(&seg_elev, &efnh[0], x, y); // elevation of pixel neighbourhood
                Segment_get(&seg_elev, &efnh[1], x, y-1);
                Segment_get(&seg_elev, &efnh[2], x+1, y-1);
                Segment_get(&seg_elev, &efnh[3], x+1, y);
                Segment_get(&seg_elev, &efnh[4], x+1, y+1);
                Segment_get(&seg_elev, &efnh[5], x, y+1);
                Segment_get(&seg_elev, &efnh[6], x-1, y+1);
                Segment_get(&seg_elev, &efnh[7], x-1, y);
                Segment_get(&seg_elev, &efnh[8], x-1, y-1);

                for (i=0; i<9; i++) if (efnh[i] == UNDEFMIN ) ctrl_valtopo = 0; // setting validity of cell to negative, if necessary

                if ( ctrl_valtopo == 1 ) { // if pixel is valid:

                    dz_dx = (( efnh[4] + 2 * efnh[3] + efnh[2] ) - ( efnh[6] + 2 * efnh[7] + efnh[8] )) / ( 8 * csz );
                        // slope in x direction
                    dz_dy = (( efnh[8] + 2 * efnh[1] + efnh[2] ) - ( efnh[6] + 2 * efnh[5] + efnh[4] )) / ( 8 * csz );
                        // slope in y direction
                    corr_topo += cos ( atan ( pow ( pow ( dz_dx, 2 ) + pow ( dz_dy, 2 ) , 0.5 ) ) );
                        // updating topographic correction factor
                    ctr_topo += 1; // updating counter for number of cells
                }

                if ( efnh[0] > UNDEFMIN ) ctrl_content = 1; // setting control for valid content to one if at least one cell is valid

                Segment_get(&seg_soilclass, &pscs, x, y); // reading soil class
                if ( pscs > 0 ) num_cellval += 1; // if soil class is valid, increasing number of valid raster cells
            }
        }
        corr_topo = corr_topo / (float)ctr_topo; // final value for topographic correction factor


    // Number of ellipsoids and related parameters

        if ( ELLDENS > 0 ) { // if more than 1 ellipsoid is to be tested:

            fprintf (fchange, "%i\n", num_cellval); // writing number of valid cells to evolution file
            interval_ctr = 0; // counter for number of simulations after full ellipsoid density
            elldens_interval = (int)(( 16 * m * n * pow ( csz, 2 ))
                / ( PI * corr_topo * ( LENGTHMAX + LENGTHMIN ) * ( WIDTHMAX + WIDTHMIN )) + 0.5 );
                // number of ellipsoids to be tested for ellipsoid density of 1

            NUM_MCSIM = (int)(ELLDENS * ( 16 * m * n * pow ( csz, 2 ))
                / ( PI * corr_topo * ( LENGTHMAX + LENGTHMIN ) * ( WIDTHMAX + WIDTHMIN )) + 0.5 ); // number of ellipsoids to be tested
        }

        else NUM_MCSIM = 1; // if only 1 ellipsoid is to be tested, setting number of simulations to 1
        if ( ctrl_content == 0 ) NUM_MCSIM = 0; // performing no simulation if no valid content was detected


    // Starting ellipsoid-based slope stability calculation

        while ( ctr_mcsim < NUM_MCSIM ) { // loop over all ellipsoids to be tested:

           ctr_mcsim += 1; // updating number of ellipsoids
           for ( j=1; j <= jmax+1; j++ ) validity[j] = 1; // initializing validity of layers


        // Resetting values

            time_stepstart = clock(); // start time

            for ( j = 1; j <= jmax; j++ ) {
                dslip_test[j] = 0; // depth of slip surface

                for ( pf_ctr = 1; pf_ctr <= pf_ctr_max; pf_ctr++ ) {
                    r_cum[j][pf_ctr] = 0; // shear resistance
                    t_cum[j][pf_ctr] = 0; // shear force
                }
            }


        // Randomization of ellipsoid parameters

            if ( ZBMAX > ZBMIN ) {

                randmax = ( int ) ( 1000 * ( ZBMAX - ZBMIN ) + 0.5 );
                ell_zb_rat = ( ( ( float ) ( rand()%randmax ) ) ) / 1000 + ZBMIN;
            }
            else ell_zb_rat = ZBMAX;

            ELLAMAX = LENGTHMAX / ( 2 * pow( 1 - ell_zb_rat, 0.5 ) );
            ELLAMIN = LENGTHMIN / ( 2 * pow( 1 - ell_zb_rat, 0.5 ) );
            ELLBMAX = WIDTHMAX / ( 2 * pow( 1 - ell_zb_rat, 0.5 ) );
            ELLBMIN = WIDTHMIN / ( 2 * pow( 1 - ell_zb_rat, 0.5 ) );

            if ( ELLAMAX > ELLAMIN ) {

                randmax = ( int ) ( 1000 * ( ELLAMAX - ELLAMIN ) + 0.5 );
                ell_a = ( ( float ) ( ( rand()%randmax ) ) ) / 1000 + ELLAMIN;
            }
            else ell_a = ELLAMAX;

            if ( ELLBMAX > ell_a / RAT_GEOM_MIN ) ELLBMAX = ell_a / RAT_GEOM_MIN;
            if ( ELLBMIN < ell_a / RAT_GEOM_MAX ) ELLBMIN = ell_a / RAT_GEOM_MAX;

            if ( WIDTHMAX < 0 ) {
                ell_b = ell_a;
            }
            else {

                if ( ELLBMAX > ELLBMIN ) {

                    randmax = ( int ) ( 1000 * ( ELLBMAX - ELLBMIN ) + 0.5 );
                    ell_b = ( ( float ) ( ( rand()%randmax ) ) ) / 1000 + ELLBMIN;
                }
                else ell_b = ELLBMAX;
            }

            ell_length = 2 * ell_a * pow( 1 - ell_zb_rat, 0.5 );
            ell_width = 2 * ell_b * pow( 1 - ell_zb_rat, 0.5 );
            ell_area = ell_length * ell_width * PI / 4;

            if ( ell_a > ell_b ) ell_max = ell_a; else ell_max = ell_b;
            if ( CENTERX == UNDEFMIN ) {
                randmax = m - 2 * ( int ) ( ( ell_max / csz ) + 0.5 );
                ell_centerx = ( int ) ( ( rand()%randmax ) + 0.5 ) + ( int ) ( ( ell_max / csz ) + 0.5 );
            }
            else ell_centerx = (int) ( ( ( bd_north - CENTERY ) / csz ) + 0.5 );
            ell_centery_metric = bd_north - ( float ) ell_centerx * csz;

            if ( CENTERY == UNDEFMIN ) {
                randmax = n - 2 * ( int ) ( ( ell_max / csz ) + 0.5 );
                ell_centery = ( int ) ( ( rand()%randmax ) + 0.5 ) + ( int ) ( ( ell_max / csz ) + 0.5 );
            }
            else ell_centery = (int) ( ( ( CENTERX - bd_west ) / csz ) + 0.5 );
            ell_centerx_metric = ( float ) ell_centery * csz + bd_west;

            if ( ELL_ALPHADEF < 0 || ELL_BETADEF < 0 ) {

                randmax = 500;
                fit = (float)(( rand()%randmax ) + 500 ) / 1000;
            }
            else fit = 1;

            if ( NEWMARK == 1 ) Segment_get(&seg_distance, &distance, ell_centerx, ell_centery);
            if ( NEWMARKREF == 1 ) Segment_get(&seg_newmarkref, &newmarkref, ell_centerx, ell_centery);
            else newmarkref = UNDEFPLUS;


        // Fitting ellipsoid into slope

            xp = ell_centerx + ( int ) ( ( fit * ell_max / csz ) + 0.5 ); // maximum extent
            xm = ell_centerx - ( int ) ( ( fit * ell_max / csz ) + 0.5 );
            yp = ell_centery + ( int ) ( ( fit * ell_max / csz ) + 0.5 );
            ym = ell_centery - ( int ) ( ( fit * ell_max / csz ) + 0.5 );

            if ( xp > m - 1 || xm < 1 || yp > n - 1 || ym < 1 ) ctrl_continue = 0; else ctrl_continue = 1; // testing validity of ellipsoid

            if ( ctrl_continue == 1 ) { // if ellipsoid is still valid:

                if ( ELL_ALPHADEF < 0 || ELL_BETADEF < 0 ) { // if orientation of ellipsoid is not defined:

                    Segment_get(&seg_elev, &efnh[1], ell_centerx, ym); // aspect and slope of main inclination direction of slip surface
                    Segment_get(&seg_elev, &efnh[2], xp, ym);
                    Segment_get(&seg_elev, &efnh[3], xp, ell_centery);
                    Segment_get(&seg_elev, &efnh[4], xp, yp);
                    Segment_get(&seg_elev, &efnh[5], ell_centerx, yp);
                    Segment_get(&seg_elev, &efnh[6], xm, yp);
                    Segment_get(&seg_elev, &efnh[7], xm, ell_centery);
                    Segment_get(&seg_elev, &efnh[8], xm, ym);

                    for ( o=1; o<9; o++ ) if ( efnh[o] == UNDEFMIN ) ctrl_continue = 0; // testing validity of ellipsoid:

                    dz1 = (( efnh[4] + 2 * efnh[3] + efnh[2] ) - ( efnh[6] + 2 * efnh[7] + efnh[8] ));
                    dz2 = (( efnh[8] + 2 * efnh[1] + efnh[2] ) - ( efnh[6] + 2 * efnh[5] + efnh[4] ));

                    dz_dx = dz1 / ( 4 * csz * (xp - xm) );
                    dz_dy =  dz2 / ( 4 * csz * (yp - ym) );

                    if ( ( dz_dx == 0 && dz_dy == 0 ) || ( xp == xm || yp == ym ) ) ctrl_continue = 0;

                    if ( ctrl_continue == 1 ) { // if ellipsoid is still valid:

                        beta = atan ( pow ( pow ( dz_dx, 2 ) + pow ( dz_dy, 2 ) , 0.5 ) );
                        beta_xz = atan ( dz_dx );
                        beta_yz = atan ( dz_dy );
                        alpha = atan ( dz_dy / -dz_dx );

                        if ( beta_xz == 0 && beta_yz > 0 ) alpha = PI / 2;
                        else if ( beta_xz == 0 && beta_yz < 0 ) alpha = 3 * PI / 2;
                        else if ( beta_xz < 0 && beta_yz == 0 ) alpha = 0;
                        else if ( beta_xz > 0 && beta_yz == 0 ) alpha = PI;
                        else if ( beta_xz == 0 && beta_yz == 0 ) alpha = 0;
                        else if ( beta_xz < 0 && beta_yz > 0) alpha = -atan (tan ( beta_yz ) / tan ( beta_xz ) );
                        else if ( beta_xz > 0 && beta_yz > 0) alpha = PI - atan ( tan ( beta_yz ) / tan ( beta_xz ) );
                        else if ( beta_xz > 0 && beta_yz < 0) alpha = PI - atan ( tan ( beta_yz ) / tan ( beta_xz ) );
                        else alpha = 2 * PI - atan ( tan ( beta_yz ) / tan ( beta_xz ) );
                    }
                }
                else { // if orientation of ellipsoid is defined by the user:

                    alpha = ELL_ALPHADEF * PI / 180;
                    beta = ELL_BETADEF * PI / 180;
                    beta_xz = atan ( tan(beta) * cos (alpha));
                    beta_yz = atan ( tan(beta) * sin (alpha));
                }

                if ( ctrl_continue == 1 ) { // if ellipsoid is still valid:

                    ELLCMAX = DEPTHMAX * cos ( beta ) / ( 1 - ell_zb_rat ); // constraints for c axis of ellipsoid
                    ELLCMIN = DEPTHMIN * cos ( beta ) / ( 1 - ell_zb_rat );

                    if ( WIDTHMAX < 0 ) ell_c = ell_a;

                    else {

                        if ( ELLCMAX > ell_a ) ELLCMAX = ell_a;
                        if ( TRUNCATED == 0 ) {

                            ell_centerxd = ell_centerx + (int) ( ELLCMAX * ( 1 - ell_zb_rat ) * sin ( beta_xz ) / csz + 0.5 );
                            ell_centeryd = ell_centery - (int) ( ELLCMAX * ( 1 - ell_zb_rat ) * sin ( beta_yz ) / csz + 0.5 );

                            if ( ell_centerxd < 1 || ell_centerxd > m-1 || ell_centeryd < 1 || ell_centeryd > n-1 ) ctrl_continue = 0;

                            if ( CLASSINF == 1 ) {

                                Segment_get(&seg_depth[nly], &arr_depthr, ell_centerxd, ell_centeryd); // reading depth of layer bottom
                                if ( ELLCMAX > arr_depthr * cos ( beta ) / ( 1 - ell_zb_rat ) )
                                    ELLCMAX = arr_depthr * cos ( beta ) / ( 1 - ell_zb_rat );
                            }
                        }
                        if ( ELLCMAX > ELLCMIN  + 0.001 ) {
                            randmax = ( int ) ( 1000 * ( ELLCMAX - ELLCMIN ) + 0.5 );
                            ell_c = ( ( float ) ( rand()%randmax ) ) / 1000 + ELLCMIN;
                        }
                        else ell_c = ELLCMAX;
                    }
                    if ( ctrl_continue == 1 ) {

                        ell_zb = ell_c * ell_zb_rat;

                        ell_centerx = ell_centerx - (int) ( ell_zb * sin ( beta_xz ) / csz + 0.5 );
                        ell_centery = ell_centery + (int) ( ell_zb * sin ( beta_yz ) / csz + 0.5 ); 

                        ell_centery_metric = bd_north - ( float ) ell_centerx * csz;
                        ell_centerx_metric = ( float ) ell_centery * csz + bd_west;

                        Segment_get(&seg_elev, &param_elev, ell_centerx, ell_centery); // z coordinate of ellipsoid centre
                        zell_0 = param_elev + ell_zb / cos ( beta );


                    // Preparing ellipsoid mask

                        xp = ell_centerx + ( int ) ( ( ell_max / csz ) + 0.5 ); // area of interest
                        xm = ell_centerx - ( int ) ( ( ell_max / csz ) + 0.5 );
                        yp = ell_centery + ( int ) ( ( ell_max / csz ) + 0.5 );
                        ym = ell_centery - ( int ) ( ( ell_max / csz ) + 0.5 );

                        for ( x = xm-1; x <= xp+1; x++ ) {
	                    for ( y = ym-1; y <= yp+1; y++ ) {
                               if ( x > m - 1 || x < 1 || y > n - 1 || y < 1 ) ctrl_continue = 0;

                                if ( ctrl_continue == 1 ) {
                                    ex = x - xm + 1;
                                    ey = y - ym + 1;

                                    Segment_get(&seg_soilclass, &pscs, x, y); // reading soil class

                                    if ( USESEG == 1 ) Segment_put(&seg_esoilclass, &pscs, ex, ey); // writing soil class
                                    else eparam_soilclass[ex][ey] = pscs;

                                    Segment_get(&seg_elev, &param_elev, x, y);
                                    epes = (int)(param_elev * CHANGEDM + 0.5);
                                    ezslips = epes;
                                    Segment_get(&seg_gwdepth, &param_gwdepth, x, y);
                                    dgws = (int)(param_gwdepth * CHANGEDM + 0.5);
                                    if ( USESEG == 1 ) Segment_put(&seg_epe, &epes, ex, ey);
                                    else epe[ex][ey] = epes;
                                    if ( USESEG == 1 ) Segment_put(&seg_dgw, &dgws, ex, ey);
                                    else dgw[ex][ey] = dgws;
                                    if ( USESEG == 1 ) Segment_put(&seg_ezslip, &ezslips, ex, ey);
                                    else ezslip[ex][ey] = ezslips;

                                    if ( PSEUDOSTATIC == 1 ) {

                                        Segment_get(&seg_pga, &pga, x, y);
                                        if ( USESEG == 1 ) Segment_put(&seg_epga, &pga, ex, ey);
                                        else epga[ex][ey] = pga;
                                    }

                                    if ( CLASSINF == 1 ) { // if soil class-based mode applies:

                                        for ( j = 0; j <= nly; j++ ) {
                                            Segment_get(&seg_depth[j], &arr_depthr, x, y); // reading depth of layer bottom from array
                                            epds = (int)( arr_depthr * CHANGEDM ); // converting to integer and writing to array
                                            if ( USESEG == 1 ) Segment_put(&seg_epd[j], &epds, ex, ey);
                                            else epd[j][ex][ey] = epds;
                                        }
                                    }
                                }
                            }
                        }

                        if ( USESEG == 1 ) {

                            Segment_flush(&seg_esoilclass);
                            Segment_flush(&seg_epe);
                            Segment_flush(&seg_dgw);
                            Segment_flush(&seg_ezslip);
                            for ( j = 0; j <= nly; j++ ) Segment_flush(&seg_epd[j]);

                            if ( PSEUDOSTATIC == 1 ) Segment_flush(&seg_epga);

                            Segment_get(&seg_esoilclass, &epsc, ell_centerx-xm+1, ell_centery-ym+1);
                            Segment_get(&seg_esoilclass, &epsc_nh[1], xp-xm+1, yp-ym+1);
                            Segment_get(&seg_esoilclass, &epsc_nh[2], xp-xm+1, 1);
                            Segment_get(&seg_esoilclass, &epsc_nh[3], 1, yp-ym+1);
                            Segment_get(&seg_esoilclass, &epsc_nh[4], 1, 1);
                        }
                        else {
                            epsc = eparam_soilclass[ell_centerx-xm+1][ell_centery-ym+1];
                            epsc_nh[1] = eparam_soilclass[xp-xm+1][yp-ym+1];
                            epsc_nh[2] = eparam_soilclass[xp-xm+1][1];
                            epsc_nh[3] = eparam_soilclass[1][yp-ym+1];
                            epsc_nh[4] = eparam_soilclass[1][1];
                        }

                        if ( epsc == UNDEFMIN ) ctrl_continue = 0; // testing whether ellipsoid is within area of interest
                        else if ( epsc_nh[1] == UNDEFMIN ) ctrl_continue = 0;
                        else if ( epsc_nh[2] == UNDEFMIN ) ctrl_continue = 0;
                        else if ( epsc_nh[3] == UNDEFMIN ) ctrl_continue = 0;
                        else if ( epsc_nh[4] == UNDEFMIN ) ctrl_continue = 0;
                    }
                }
            }

            if ( ctrl_continue == 1 ) { // if ellipsoid is still valid:


            // Resetting layer variables for local reference system

                if ( LAYERS == 1 ) { // for layer mode:

                    for (j=0; j<=nly; j++) ldata_trans[j]=UNDEFMIN; // resetting array for relating local to global layer id

                    for ( ex = 0; ex < me; ex++ ) { // loop over all cells in local reference system:
                        for ( ey = 0; ey < ne; ey++ ) {
                            for ( j=0; j <= nly+1; j++ ) { // loop over all layers:

                                epds = UNDEFLYR; // resetting depth of layer in local reference system
                                if ( USESEG == 1 ) Segment_put(&seg_epd[j], &epds, ex, ey);
                                else epd[j][ex][ey] = epds;
                            }
                        }
                    }
                    if ( USESEG == 1 ) { for ( j = 0; j <= nly; j++ ) Segment_flush(&seg_epd[j]); }


                // Converting layers into local reference system

                    linit = 1; // resetting control for top layer
                    lglob = 0; // resetting global layer id
                    lloc = 0; // resetting local layer id

                    for ( k = 0; k < nlinesl; k++ ) { // loop over all lines of the layer array:

                        if ( ldata[k][0] >= xm && ldata[k][0] <= xp && ldata[k][1] >= ym && ldata[k][1] <= yp ) {
                            // if coordinates are relevant for ellipsoid:

                            ex = ldata[k][0] - xm + 1; // local coordinates
                            ey = ldata[k][1] - ym + 1;

                            if ( linit == 1 || ldata[k][2] > lglob ) { // if layer is top layer or global layer id has increased:

                                lloc += 1; // increasing local layer id

                                lglob = ldata[k][2]; // updating global layer id
                                linit = 0; // setting control for top layer to negative

                                if ( lloc < nly ) {

                                    pgamma[0][lloc] = pgamma0[0][lglob]; // storing soil parameters in local reference system
                                    pcsoil[0][lloc] = pcsoil0[0][lglob];
                                    pphi[0][lloc] = pphi0[0][lglob];
                                    ptheta[0][lloc] = ptheta0[0][lglob];

                                    if ( PF == 1 ) { // for slope failure probability:

                                        pcsoils[0][lloc] = pcsoils0[0][lglob];
                                        pcsoilmin[0][lloc] = pcsoilmin0[0][lglob];
                                        pcsoilmax[0][lloc] = pcsoilmax0[0][lglob];
                                        pphis[0][lloc] = pphis0[0][lglob];
                                        pphimin[0][lloc] = pphimin0[0][lglob];
                                        pphimax[0][lloc] = pphimax0[0][lglob];

                                        pf_pdf_peak_phi[0][lloc] = pf_pdf_peak_phi0[0][lglob];
                                        pf_pdf_peak_csoil[0][lloc] = pf_pdf_peak_csoil0[0][lglob];
                                        for ( pf_ctr = 1; pf_ctr <= pf_nsteps_phi; pf_ctr++ )
                                            pf_phi[pf_ctr][0][lloc] = pf_phi0[pf_ctr][0][lglob];
                                        for ( pf_ctr = 1; pf_ctr <= pf_nsteps_csoil; pf_ctr++ )
                                            pf_csoil[pf_ctr][0][lloc] = pf_csoil0[pf_ctr][0][lglob];
                                    }
                                }
                            }

                            if ( lloc < nly ) ldata_trans[lloc] = ldata[k][2]; // relating local to global layer id
                            if ( lloc < nly ) {
                                if ( USESEG == 1 ) Segment_put(&seg_epd[lloc], &ldata[k][3], ex, ey);
                                    // storing depth of layer in local reference system
                                else epd[lloc][ex][ey] = ldata[k][3];
                            }
                        }
                    }
                    if ( USESEG == 1 ) { for ( j = 0; j <= lloc; j++ ) Segment_flush(&seg_epd[j]); }

                    lbtm = lloc + 1; // storing bottom layer id
                }

                if ( LAYERS == 0 ) lbtm = nly; // if soil class-based mode is applied, setting layer bottom to total number of layers
                if ( PFDEPTH == 1 ) jmax = pf_nsteps_depth;
                    // for slope failure probability with depth sampling, setting number of layers to number of sampled depth values
                else jmax = lbtm;

                if ( lbtm > lbtm_max) lbtm_max = lbtm; // updating maximum number of layers in local reference system
                if ( LAYERS == 1 && nly <= lbtm_max ) ctrl_continue = 0; // updating validity of ellipsoid
                if ( LAYERS == 0 && nly < lbtm_max ) ctrl_continue = 0;

                if ( ctrl_continue == 1 ) { // if ellipsoid is still valid:

                    if ( LAYERS == 1 ) {
                        for ( j=1; j<lbtm; j++ ) validity[j] = 0; // resetting validity of layer
                    }

                    for ( x = xm; x <= xp; x++ ) { // loop over all pixels of ellipsoid mask:
                        for ( y = ym; y <= yp; y++ ) {

                            ex = x - xm + 1; // local coordinates
                            ey = y - ym + 1;

                            if ( LAYERS == 1 ) { // if layer-based mode is applied:

                                epds = UNDEFPLUS * CHANGEDM; // storing depth of artificial bottom layer in local reference system
                                if ( USESEG == 1 ) Segment_put(&seg_epd[lbtm], &epds, ex, ey);
                                else epd[lbtm][ex][ey] = epds;
                                if ( USESEG == 1 ) Segment_flush(&seg_epd[lbtm]);


                            // Building consistent layer depth dataset in local reference system

                                linit = 1; // resetting control for top layer

                                for ( j=1; j < lbtm; j++ ) { // loop over all layers (local reference system):

                                    area_eff[j] = 0; // resetting slip surface area

                                    if ( USESEG == 1 ) Segment_get(&seg_epd[j], &epds, ex, ey);
                                    else epds = epd[j][ex][ey];
                                    if ( epds != UNDEFLYR ) linit = 0; // setting control for top layer to negative if layer depth is valid
                                    else if ( epds == UNDEFLYR && linit == 1 ) epds = 0;
                                        // setting depth to zero if top layer has not been reached
                                    else epds = UNDEFPLUS * CHANGEDM;
                                        // setting depth to high undefined if layer depth is not valid, but top layer has been reached

                                    if ( j > 1 ) {
                                        if ( USESEG == 1 ) Segment_get(&seg_epd[j-1], &epd1s, ex, ey);
                                        else epd1s = epd[j-1][ex][ey];
                                        if ( epds < epd1s ) epds = epd1s; // avoiding decreasing depth with increasing layer id
                                    }
                                    if ( USESEG == 1 ) Segment_put(&seg_epd[j], &epds, ex, ey);
                                    else epd[j][ex][ey] = epds;
                                    if ( USESEG == 1 ) Segment_flush(&seg_epd[j]);
                                }
                            }


                        // Transformation of coordinate system

                            xdiff = csz * ( ell_centerx - x );
                            ydiff = csz * ( ell_centery - y );

                            xloc1_ell = -xdiff * cos ( alpha ) - ydiff * sin ( alpha );
                            yloc1_ell = -xdiff * sin ( alpha ) + ydiff * cos ( alpha );

                            xloc2_ell = xloc1_ell * ( cos (beta ) + tan ( beta ) * sin ( beta ) );
                            yloc2_ell = yloc1_ell;


                        // Parameters for quadratic equation

                            qeq_a = 1 / pow ( ell_a , 2 ) + 1 / ( pow ( ell_c , 2 ) * pow ( tan ( beta ) , 2 ) );
                            qeq_b = 2 * xloc2_ell / pow ( ell_a , 2 );
                            qeq_c = pow ( xloc2_ell , 2 ) / pow ( ell_a , 2 ) + pow ( yloc2_ell , 2 ) / pow ( ell_b , 2 ) - 1;


                        // Solution of quadratic equation and determination of slip surface elevation

                            if ( USESEG == 1 ) Segment_get(&seg_epe, &epes, ex, ey);
                            else epes = epe[ex][ey];

                            if ( pow ( qeq_b , 2 ) > 4 * qeq_a * qeq_c ) {

                                delta_xloc2 = ( - qeq_b + pow ( pow ( qeq_b , 2 ) - 4 * qeq_a * qeq_c , 0.5 ) ) / ( 2 * qeq_a );
                                zloc2_ell = -delta_xloc2 / tan ( beta );

                                ezslips = (int)( CHANGEDM * ( zell_0 + ( zloc2_ell - xloc1_ell * sin ( beta ) ) / cos ( beta ) ) + 0.5 );
                            }
                            else ezslips=epes;

                            if ( ezslips > epes ) ezslips = epes;
                            epds = epes - ezslips;
                            if ( USESEG == 1 ) Segment_put(&seg_epd[lbtm+1], &epds, ex, ey);
                            else epd[lbtm+1][ex][ey] = epds;

                            if ( TRUNCATED == 0 ) {
                                if ( USESEG == 1 ) Segment_get(&seg_epd[lbtm], &epds, ex, ey);
                                else epds = epd[lbtm][ex][ey];
                                if ( ezslips < epes - epds ) ctrl_continue = 0;
                            }
                            if ( USESEG == 1 ) Segment_put(&seg_ezslip, &ezslips, ex, ey);
                            else ezslip[ex][ey] = ezslips;

                            if ( LAYERS == 1 ) {
                                for ( j=1; j<lbtm; j++ ) { // loop over all layers (local reference system):

                                    if ( USESEG == 1 ) Segment_get(&seg_epd[j], &epds, ex, ey);
                                    else epds = epd[j][ex][ey];
                                    if ( epds < epes - ezslips && epds > 0 ) validity[j] = 1; // updating validity of layer
                                }
                            }
                        }
                    }
                    if ( USESEG == 1 ) Segment_flush(&seg_epd[lbtm+1]);
                    if ( USESEG == 1 ) Segment_flush(&seg_ezslip);
                }


            // Computing factor of safety

                if ( ctrl_continue == 1 ) { // if ellipsoid is still valid:


                // Defining sample of parameters needed for slope failure probability

                    if ( PF == 1 && PFSAMPLING <= 1 ) { // for random sampling:

                        pf_ctr0 = 0;
                        randmax = 1000000;


                    // Probability density for truncated depth

                        if ( PFDEPTH == 1 ) {

                            my_depth = log ( pow ( pf_depth_mean, 2 ) / pow ( pow ( pf_depth_mean, 2 ) + pow ( pf_depth_stdev, 2 ) , 0.5 ) );
                            sigma_depth = pow ( log( 1 + pow ( pf_depth_stdev, 2 ) / pow ( pf_depth_mean, 2 ) ) , 0.5 );

                            pf_ctr = 0;
                            while ( pf_ctr < pf_nsteps_depth ) {

                                pf_ctr0 += 1;
                                pf_depth0 = ( float )( rand()%randmax );
                                pf_depth_test = pf_depth_min + pf_depth_int * pf_depth0 / (float)randmax;

                                if ( pdftype_depth == 0 ) pf_pdf_depth = 1; // rectangular distribution
                                else if ( pdftype_depth == 1 ) pf_pdf_depth = 1 / ( pf_depth_stdev * pow ( 2 * PI, 0.5 ) )
                                   * exp ( -pow( pf_depth_test - pf_depth_mean, 2 ) / ( 2 * pow ( pf_depth_stdev, 2 )));
                                   // normal distribution
                                else if ( pdftype_depth == 2 ) // log-normal distribution
                                    pf_pdf_depth = 1 / ( pf_depth_test * sigma_depth * pow( 2 * PI, 0.5 ))
                                        * exp ( -( pow( log( pf_depth_test ) - my_depth, 2 ) / ( 2 * pow( sigma_depth , 2 ))));
                                else pf_pdf_depth = ( 1 / pf_depth_mean * exp ( -( 1 / pf_depth_mean ) * pf_depth_test ));
                                    // exponential distribution
                                pf_ctr0 += 1; // random probability density
                                pf_pdf0 = ( float )( rand()%randmax );
                                pf_pdf_test = pf_pdf0 * pf_pdf_peak_depth / (float)randmax;

                                if ( pf_pdf_test <= pf_pdf_depth ) {
                                    // comparison of random probability density and computed probability density
                                    pf_ctr += 1;
                                    pf_depth[pf_ctr] = (int)( CHANGEDM * pf_depth_test + 0.5 );
                                }
                            }
                            qsort(pf_depth, pf_nsteps_depth+1, sizeof pf_depth[0], pf_comp); // sorting sample by ascending depth
                        }

                        if ( CLASSINF == 1 ) lbtm2 = lbtm;
                        else lbtm2 = lbtm - 1;

                        for ( i=0; i<=nsoilclass; i++ ) { // loop over all soil classes and layers:
                            for ( j=0; j<=lbtm2; j++ ) {

                                if ( ( CLASSINF == 1 && i > 0 && j > 0 ) || ( LAYERS == 1 && ( i > 0 || j > 0 ) && ( i == 0 || j == 0 ) ) ) {

                                    pf_phi_mean = pphi[i][j]; // statistics of internal friction angle (phi)
                                    pf_phi_stdev = pphis[i][j];
                                    pf_phi_min = pphimin[i][j];
                                    pf_phi_max = pphimax[i][j];

                                    pf_csoil_mean = pcsoil[i][j]; // statistics of cohesion (csoil)
                                    pf_csoil_stdev = pcsoils[i][j];
                                    pf_csoil_min = pcsoilmin[i][j];
                                    pf_csoil_max = pcsoilmax[i][j];

                                    pf_phi_int = pf_phi_max - pf_phi_min;
                                    pf_csoil_int = pf_csoil_max - pf_csoil_min;

                                    my_phi = log ( pow ( pf_phi_mean, 2 ) / pow ( pow ( pf_phi_mean, 2 ) + pow ( pf_phi_stdev, 2 ) , 0.5 ) );
                                    sigma_phi = pow ( log( 1 + pow ( pf_phi_stdev, 2 ) / pow ( pf_phi_mean, 2 ) ) , 0.5 );
                                    my_csoil = log ( pow ( pf_csoil_mean, 2 )
                                       / pow ( pow ( pf_csoil_mean, 2 ) + pow ( pf_csoil_stdev, 2 ) , 0.5 ) );
                                    sigma_csoil = pow ( log( 1 + pow ( pf_csoil_stdev, 2 ) / pow ( pf_csoil_mean, 2 ) ) , 0.5 );


                                // Probability density for internal friction angle (phi)

                                    if ( pf_nsteps_phi > 1 ) {

                                        pf_ctr = 0;
                                        while ( pf_ctr < pf_nsteps_phi ) {

                                            pf_ctr0 += 1;
                                            pf_phi00 = ( float )( rand()%randmax );
                                            pf_phi_test = pf_phi_min + pf_phi_int * pf_phi00 / (float)randmax;

                                            if ( pdftype_phi == 0 ) pf_pdf_phi = 1; // rectangular distribution
                                            else if ( pdftype_phi == 1 ) pf_pdf_phi = 1 / ( pf_phi_stdev * pow ( 2 * PI, 0.5 ) )
                                                * exp ( -pow( pf_phi_test - pf_phi_mean, 2 ) / ( 2 * pow ( pf_phi_stdev, 2 )));
                                                // normal distribution
                                            else if ( pdftype_phi == 2 ) // log-normal distribution
                                                pf_pdf_phi = 1 / ( pf_phi_test * sigma_phi * pow( 2 * PI, 0.5 ))
                                                    * exp ( -( pow( log( pf_phi_test ) - my_phi, 2 ) / ( 2 * pow( sigma_phi , 2 ))));
                                            else pf_pdf_phi = ( 1 / pf_phi_mean * exp ( -( 1 / pf_phi_mean ) * pf_phi_test ));
                                                // exponential distribution
                                            pf_ctr0 += 1; // random probability density
                                            pf_pdf0 = ( float )( rand()%randmax );
                                            pf_pdf_test = pf_pdf0 * pf_pdf_peak_phi[i][j] / (float)randmax;

                                            if ( pf_pdf_test <= pf_pdf_phi ) {
                                                // comparison of random probability density and computed probability density
                                                pf_ctr += 1;
                                                pf_phi[pf_ctr][i][j] = (int)( 10000 * pf_phi_test + 0.5 );
                                            }
                                        }
                                        
                                    } else pf_phi[1][i][j] = (int)( 10000 * pf_phi_mean + 0.5 );


                                // Probability density for cohesion (csoil)

                                    if ( pf_nsteps_csoil > 1 ) {

                                        pf_ctr = 0;
                                        while ( pf_ctr < pf_nsteps_csoil ) {

                                            pf_ctr0 += 1;
                                            pf_csoil00 = ( float )( rand()%randmax );
                                            pf_csoil_test = pf_csoil_min + pf_csoil_int * pf_csoil00 / (float)randmax;

                                            if ( pdftype_csoil == 0 ) pf_pdf_csoil = 1; // rectangular distribution
                                            else if ( pdftype_csoil == 1  || pdftype_csoil == 4)
                                                pf_pdf_csoil = 1 / ( pf_csoil_stdev * pow ( 2 * PI, 0.5 ) )
                                                * exp ( -pow( pf_csoil_test - pf_csoil_mean, 2 ) / ( 2 * pow ( pf_csoil_stdev, 2 )));
                                                // normal distribution
                                            else if ( pdftype_csoil == 2 ) // log-normal distribution
                                                pf_pdf_csoil = 1 / ( pf_csoil_test * sigma_csoil * pow( 2 * PI, 0.5 ))
                                                    * exp ( -( pow( log( pf_csoil_test ) - my_csoil, 2 ) / ( 2 * pow( sigma_csoil , 2 ))));
                                            else pf_pdf_csoil = ( 1 / pf_csoil_mean * exp ( -( 1 / pf_csoil_mean ) * pf_csoil_test ));
                                                // exponential distribution
                                            pf_ctr0 += 1; // random probability density
                                            pf_pdf0 = ( float )( rand()%randmax );
                                            pf_pdf_test = pf_pdf0 * pf_pdf_peak_csoil[i][j] / (float)randmax;

                                            if ( pf_pdf_test <= pf_pdf_csoil ) {
                                                // comparison of random probability density and computed probability density
                                                pf_ctr += 1;
                                                pf_csoil[pf_ctr][i][j] = (int)( pf_csoil_test + 0.5 );
                                            }
                                        }
                                        
                                    } else pf_csoil[1][i][j] = (int)( pf_csoil_mean + 0.5 );
                                }
                            }
                        }
                    }


                // Defining profile for profile file

                    if ( ELLDENS == 0 ) {

                        s2_ctrl = 0;

                        if ( alpha == PI/2 || alpha == 3*PI/2 ) {

                            s2_x0 = ell_centerx;
                            s2_y0 = ym;
                            s2_dprex = 0;
                            s2_dprey = 1;
                            s2_d = 1;
                            s2_x[s2_ctrl] = s2_x0; s2_y[s2_ctrl] = s2_y0;
                        }
                        else if ( alpha == 0 || alpha == PI ) {

                            s2_x0 = xm;
                            s2_y0 = ell_centery;
                            s2_dprex = 1;
                            s2_dprey = 0;
                            s2_d = 1;
                            s2_x[s2_ctrl] = s2_x0; s2_y[s2_ctrl] = s2_y0;
                        }
                        else {
                            s2_dxtest = fabs( (float)( ell_centerx - xm ) / cos ( alpha ));
                            s2_dytest = fabs( (float)( ym - ell_centery ) / sin ( alpha ));

                            if ( s2_dxtest < s2_dytest ) {

                                s2_x0 = xm;

                                if (( alpha > 0 && alpha < PI/2 ) || ( alpha > PI && alpha < 3 * PI/2 ))
                                    s2_y0 = ell_centery - (int)( s2_dxtest * fabs(sin (alpha)) + 0.5 );
                                else
                                    s2_y0 = ell_centery + (int)( s2_dxtest * fabs(sin (alpha)) + 0.5 );

                                s2_dprex = 1;
                                if (( alpha > 0 && alpha < PI/2 ) || ( alpha > PI && alpha < 3 * PI/2 )) s2_dprey = 1;
                                else s2_dprey = -1;
                            }
                            else {

                                if (( alpha > 0 && alpha < PI/2 ) || ( alpha > PI && alpha < 3 * PI/2 ))
                                    s2_x0 = ell_centerx - (int)( s2_dytest * fabs(cos (alpha)) + 0.5 );
                                else
                                    s2_x0 = ell_centerx + (int)( s2_dytest * fabs(cos (alpha)) + 0.5 );
                                s2_y0 = ym;
                                if (( alpha > 0 && alpha < PI/2 ) || ( alpha > PI && alpha < 3 * PI/2 )) s2_dprex = 1;
                                else s2_dprex = -1;
                                s2_dprey = 1;
                            }
                            if ( fabs(1 / cos (alpha)) < fabs(1 / sin (alpha)) )
                                s2_d = fabs(1 / cos (alpha)); else s2_d = fabs(1 / sin (alpha));

                            s2_x[s2_ctrl] = s2_x0; s2_y[s2_ctrl] = s2_y0;
                        }
                        s2_dcum[0] = 0;

                        while ( s2_x[s2_ctrl] <= xp && s2_y[s2_ctrl] <= yp && s2_x[s2_ctrl] >= xm && s2_y[s2_ctrl] >= ym ) {

                            s2_ctrl += 1;
                            s2_x[s2_ctrl] = s2_x0 + s2_dprex * (int)( s2_ctrl * s2_d * fabs(cos (alpha)) + 0.5 );
                            s2_y[s2_ctrl] = s2_y0 + s2_dprey * (int)( s2_ctrl * s2_d * fabs(sin (alpha)) + 0.5 );
                            s2_dcum[s2_ctrl] = s2_ctrl * s2_d;
                        }
                    }


                // Pixel-specific shear resistances and shear forces for each layer and, for slope failure probability, for each sample

                    for ( ex = 1; ex <= (xp-xm+1); ex++ ) { // loop over all pixels of the ellipsoid area:
	                for ( ey = 1; ey <= (yp-ym+1); ey++ ) {

                            x = ex + xm - 1; // transforming coordinates
                            y = ey + ym - 1;

                            xnh[0]=ex; xnh[1]=ex; xnh[2]=ex+1; xnh[3]=ex+1; xnh[4]=ex+1; xnh[5]=ex; xnh[6]=ex-1; xnh[7]=ex-1; xnh[8]=ex-1;
                            ynh[0]=ey; ynh[1]=ey-1; ynh[2]=ey-1; ynh[3]=ey; ynh[4]=ey+1; ynh[5]=ey+1; ynh[6]=ey+1; ynh[7]=ey; ynh[8]=ey-1;
                                // defining neighbourhood

                            if ( USESEG == 1 ) { // defining elevation of terrain and slip surface, and soil class:
                                Segment_get(&seg_epe, &epes, ex, ey);
                                Segment_get(&seg_dgw, &dgws, ex, ey);
                                Segment_get(&seg_ezslip, &ezslips, ex, ey);
                                Segment_get(&seg_esoilclass, &epsc, ex, ey);
                            }
                            else {
                                epes = epe[ex][ey];
                                dgws = dgw[ex][ey];
                                ezslips = ezslip[ex][ey];
                                epsc = eparam_soilclass[ex][ey];
                            }

                            if ( PSEUDOSTATIC == 1 ) {

                                if ( USESEG == 1 ) Segment_get(&seg_epga, &epgas, ex, ey);
                                else epgas = epga[ex][ey];
                            }

                            if ( epsc == UNDEFMIN ) ctrl_continue = 0; // updating validity of ellipsoid
                            if ( ctrl_continue == 1 ) {
                                i = epsc;


                            // Extracting data for profile file

                                if ( ELLDENS == 0 ) { // if only one ellipsoid is tested:

                                    for ( s2_ctrl2 = 0; s2_ctrl2 <= s2_ctrl; s2_ctrl2++ ) {

                                        if ( x == s2_x[s2_ctrl2] && y == s2_y[s2_ctrl2] ) {

                                            s2_data[s2_ctrl2][0][0] = s2_dcum[s2_ctrl2] * csz;
                                            s2_data[s2_ctrl2][0][1] = ((float)epes) / CHANGEDM;
                                            s2_data[s2_ctrl2][0][2] = 0;
                                            s2_data[s2_ctrl2][0][3] = 0;
                                            s2_data[s2_ctrl2][0][4] = 0;
                                            s2_data[s2_ctrl2][0][5] = 0;
                                            s2_data[s2_ctrl2][0][6] = 0;
                                            s2_data[s2_ctrl2][0][7] = 0;

                                            for ( j=1; j <= nly; j++ ) {

                                                if ( USESEG == 1 ) Segment_get(&seg_epd[j], &epds, ex, ey);
                                                else epds = epd[j][ex][ey];
                                                if ( epds == UNDEFPLUS * CHANGEDM || epds == 0 ) s2_data[s2_ctrl2][j][1] = 0;
                                                else s2_data[s2_ctrl2][j][1] = ((float)(epes-epds)) / CHANGEDM;
                                                s2_data[s2_ctrl2][j][2] = (float)(epes) / CHANGEDM;
                                            }

                                            break;
                                        }
                                    }
                                }

                                for ( q=1; q<9; q++ ) { // elevation, depth of groundwater table, and ellipsoid bottom of neighbouring pixels:

                                    if ( USESEG == 1 ) Segment_get(&seg_epe, &epe_nhs[q], xnh[q], ynh[q]);
                                    else epe_nhs[q] = epe[xnh[q]][ynh[q]];
                                    if ( USESEG == 1 ) Segment_get(&seg_dgw, &dgw_nhs[q], xnh[q], ynh[q]);
                                    else dgw_nhs[q] = dgw[xnh[q]][ynh[q]];
                                    if ( USESEG == 1 ) Segment_get(&seg_ezslip, &ezslip_nhs[q], xnh[q], ynh[q]);
                                    else ezslip_nhs[q] = ezslip[xnh[q]][ynh[q]];
                                }

                                if ( ezslips < epes ) { // if ellipsoid bottom is below terrain surface:


                                // Aspect and slope of groundwater table

                                    if ( SEEPAGE <= 1 ) {
                                        // x and y slopes of groundwater table (assumed parallel to the terrain surface or groundwater table):

                                        for ( q=1; q<9; q++ )  efnh[q] = ((float)(epe_nhs[q]-dgw_nhs[q]))/CHANGEDM;
                                        dz_dx = (( efnh[4] + 2 * efnh[3] + efnh[2] ) - ( efnh[6] + 2 * efnh[7] + efnh[8] )) / ( 8 * csz );
                                        dz_dy = (( efnh[8] + 2 * efnh[1] + efnh[2] ) - ( efnh[6] + 2 * efnh[5] + efnh[4] )) / ( 8 * csz );
                                    }
                                    else { // x and y slopes of groundwater table (assumed parallel to average layer surface):
                                        Segment_get(&seg_dzdx, &dz_dx, x, y);
                                        Segment_get(&seg_dzdy, &dz_dy, x, y);
                                    }

                                    if ( dz_dx == 0 && dz_dy == 0 ) { beta_w = 0; alpha_w = 0; }
                                    else {
                                        beta_xz = atan ( dz_dx );
                                        beta_yz = atan ( dz_dy );
                                        beta_w = atan ( pow ( pow ( dz_dx, 2 ) + pow ( dz_dy, 2 ) , 0.5 ) );

                                        if ( beta_xz == 0 && beta_yz > 0 ) alpha_w = PI / 2;
                                        else if ( beta_xz == 0 && beta_yz < 0 ) alpha_w = 3 * PI / 2;
                                        else if ( beta_xz < 0 && beta_yz == 0 ) alpha_w = 0;
                                        else if ( beta_xz > 0 && beta_yz == 0 ) alpha_w = PI;
                                        else if ( beta_xz == 0 && beta_yz == 0 ) alpha_w = 0;
                                        else if ( beta_xz < 0 && beta_yz > 0) alpha_w = - atan (tan ( beta_yz ) / tan ( beta_xz ) );
                                        else if ( beta_xz > 0 && beta_yz > 0) alpha_w = PI - atan ( tan ( beta_yz ) / tan ( beta_xz ) );
                                        else if ( beta_xz > 0 && beta_yz < 0) alpha_w = PI - atan ( tan ( beta_yz ) / tan ( beta_xz ) );
                                        else alpha_w = 2 * PI - atan ( tan ( beta_yz ) / tan ( beta_xz ) );
                                    }

                                    weight_sub = 0; // resetting submerged weight
                                    u_appr = 0; // resetting weight of water
                                    dslip_test[0] = 0; // setting depth of slip surface for terrain surface to zero

                                    for ( j = 1; j <= jmax; j++ ) { // loop over all layers:


                                    // Slip surface parameters

                                        if ( CLASSINF == 1 ) ipar = i; // applying parameters defined for soil class, if appropriate
                                        else ipar = 0; // if layer is defined, applying parameters defined for layer

                                        if ( PFDEPTH == 1 ) epds = pf_depth[j];
                                        else if ( USESEG == 1 ) Segment_get(&seg_epd[j], &epds, ex, ey);
                                        else epds = epd[j][ex][ey];

                                        if ( epes - epds < ezslips ) { // if ellipsoid bottom is above layer bottom:

                                            if ( USESEG == 1 ) Segment_get(&seg_epd[lbtm+1], &epd1s, ex, ey);
                                            else epd1s = epd[lbtm+1][ex][ey];

                                            dslip_test[j] = epd1s; // setting slip surface depth to depth of ellipsoid bottom

                                            for ( j2 = 1; j2 <= j; j2++ ) { // loop over all layers:

                                                if ( PFDEPTH == 1 ) epd2s = pf_depth[j2];
                                                else if ( USESEG == 1 ) Segment_get(&seg_epd[j2], &epd2s, ex, ey);
                                                else epd2s = epd[j2][ex][ey];
                                                if ( epes - epd1s >= epes - epd2s ) {  // searching for slip surface layer:
                                                    j3 = j2; // layer for gamma and theta
                                                    j4 = j2; // layer for c and phi
                                                    break;
                                                }
                                            }

                                            if ( LAYERS == 1 && epds == UNDEFPLUS * CHANGEDM && j3 == lbtm ) {
                                                j3 = 0; // if layer is undefined, applying parameters for soil class
                                                j4 = 0;
                                                ipar=i;
                                            }
                                            //jctrl = 0;
                                        }
                                        else { // if layer bottom is above ellipsoid bottom:

                                            dslip_test[j] = epds; // setting slip surface depth to layer bottom
                                            j3 = j; // layer for gamma and theta
                                            j4 = j; // layer for c and phi
                                            //jctrl = 1;
                                        }
                                        if ( PFDEPTH == 1 ) {
                                            j3 = 1;
                                            j4 = 1;
                                        }


                                    // Aspect and slope of slip surface

                                        for ( q=1; q<9; q++ ) { // slip surface elevation of neighbourhood:

                                            if ( PFDEPTH == 1 ) epd_nhs[q] = pf_depth[j];
                                            else if ( USESEG == 1 ) Segment_get(&seg_epd[j], &epd_nhs[q], xnh[q], ynh[q]);
                                            else epd_nhs[q] = epd[j][xnh[q]][ynh[q]];

                                            if ( USESEG == 1 ) Segment_get(&seg_epd[lbtm+1], &epd_nhbs[q], xnh[q], ynh[q]);
                                            else epd_nhbs[q] = epd[lbtm+1][xnh[q]][ynh[q]];

                                            enh[q] = epe_nhs[q] - epd_nhs[q];
                                            if ( enh[q] < ezslip_nhs[q] ) enh[q] = epe_nhs[q] - epd_nhbs[q];
                                        }

                                        dz1 = ((float)((( enh[4] + 2 * enh[3] + enh[2] ) - ( enh[6] + 2 * enh[7] + enh[8] )))) / CHANGEDM;
                                        dz2 = ((float)((( enh[8] + 2 * enh[1] + enh[2] ) - ( enh[6] + 2 * enh[5] + enh[4] )))) / CHANGEDM;

                                        dz_dx = dz1 / ( 8 * csz);
                                        dz_dy = dz2 / ( 8 * csz);

                                        if ( dz_dx == 0 && dz_dy == 0 ) { beta_c = 0; alpha_c = 0; beta_m = 0; }
                                        else {
                                            beta_xz = atan ( dz_dx );
                                            beta_yz = atan ( dz_dy );
                                            beta_c = atan ( pow ( pow ( dz_dx, 2 ) + pow ( dz_dy, 2 ) , 0.5 ) );

                                            if ( beta_xz == 0 && beta_yz > 0 ) alpha_c = PI / 2;
                                            else if ( beta_xz == 0 && beta_yz < 0 ) alpha_c = 3 * PI / 2;
                                            else if ( beta_xz < 0 && beta_yz == 0 ) alpha_c = 0;
                                            else if ( beta_xz > 0 && beta_yz == 0 ) alpha_c = PI;
                                            else if ( beta_xz == 0 && beta_yz == 0 ) alpha_c = 0;
                                            else if ( beta_xz < 0 && beta_yz > 0) alpha_c = - atan (tan ( beta_yz ) / tan ( beta_xz ) );
                                            else if ( beta_xz > 0 && beta_yz > 0) alpha_c = PI - atan ( tan ( beta_yz ) / tan ( beta_xz ) );
                                            else if ( beta_xz > 0 && beta_yz < 0) alpha_c = PI - atan ( tan ( beta_yz ) / tan ( beta_xz ) );
                                            else alpha_c = 2 * PI - atan ( tan ( beta_yz ) / tan ( beta_xz ) );

                                            beta_m = atan ( tan ( beta_c ) * cos ( alpha_c - alpha ) );
                                        }

                                    // Area of slip surface

                                        area = pow ( csz , 2 ) * pow ( 1 - pow ( sin ( beta_xz ) , 2 ) * pow ( sin ( beta_yz ) , 2 ) , 0.5 )
                                            / ( cos ( beta_xz ) * cos ( beta_yz ) );


                                    // Submerged weight and approximation of pore water pressure related to slip surface

                                        if ( dgws >= dslip_test[j] ) depth_layer = 0; // saturated depth of layer:
                                        else if ( dgws > dslip_test[j-1] ) depth_layer = ((float)(dslip_test[j] - dgws)) / CHANGEDM;
                                        else depth_layer = ((float)(dslip_test[j] - dslip_test[j-1])) / CHANGEDM;
                                        
                                        weight_sub += pgamma[ipar][j3] * ((float)(dslip_test[j] - dslip_test[j-1])) / CHANGEDM * pow ( csz , 2 );

                                        if ( ptheta[ipar][j3] > 0 && depth_layer > 0 ) {

                                            weight_sub += GAMMA_WATER * ( ptheta[ipar][j3] - 1 ) * ( depth_layer ) * pow ( csz , 2 );
                                            u_appr += GAMMA_WATER * ( depth_layer ) * pow ( csz , 2 );
                                        }
                                        else u_appr = 0;


                                    // Seepage forces

                                        if ( ptheta[ipar][j3] > 0 && beta_w != 0 ) { // for moist soil:

                                            f_s = u_appr * sin ( beta_w );

                                            if ( f_s > 0 ) { // if seepage force is larger than zero:

                                                f_sh = f_s * cos ( beta_w );
                                                f_sv = f_s * sin ( beta_w );

                                                f_sch = f_sh * cos ( alpha_w - alpha_c );
                                                f_smh = f_sh * cos ( alpha_w - alpha );

                                                f_sc = pow ( pow ( f_sv, 2 ) + pow ( f_sch, 2 ) , 0.5 );
                                                f_sm = pow ( pow ( f_sv, 2 ) + pow ( f_smh, 2 ) , 0.5 );

                                                beta_fsc = acos ( f_sch / f_sc );
                                                beta_fsm = acos ( f_smh / f_sm ); 

                                                n_s = f_sc * sin ( beta_fsc - beta_c ); // contribution of seepage force to shear resistance
                                                t_s = f_sm * cos ( beta_fsm - beta_m ); // contribution of seepage force to shear force
                                            }
                                            else { n_s = 0; t_s = 0; } // else, setting seepage forces to zero
                                        }
                                        else { n_s = 0; t_s = 0; } // setting seepage forces to zero for dry soil


                                    // Stabilizing and destabilizing forces

                                        if ( dslip_test[j] > 0 ) { // if depth of slip surface is larger than zero:

                                            if ( PF == 0 ) { // for factor of safety:
                                                phi=pphi[ipar][j4]; // defining phi and csoil
                                                csoil=pcsoil[ipar][j4];
                                            }
                                            else { // for slope failure probability:
                                                pf_ctr_phi = 1; // initializing counters for looping over samples of phi and csoil
                                                pf_ctr_csoil = 1;
                                            }
                                            for ( pf_ctr = 1; pf_ctr <= pf_ctr_max; pf_ctr++ ) {
                                                // starting loop over samples of phi and csoil:

                                                if ( PF == 1 && PFSAMPLING >= 1 ) {
                                                    // for slope failure probability with sampling of parameters:

                                                    phi = ( (float)(pf_phi[pf_ctr_phi][ipar][j4]) ) / 10000; // defining phi

                                                    if ( pdftype_csoil == 4 ) { // defining csoil:
                                                        csoil = regpar1 * phi * 180 / PI + regpar2 + (float)pf_csoil[pf_ctr_csoil][ipar][j4];
                                                        if ( csoil < 0 ) csoil = 0;
                                                    }
                                                    else csoil = (float)pf_csoil[pf_ctr_csoil][ipar][j4];

                                                    //if ( jctrl == 1 ) { csoil = 0; phi -= 10; }

                                                    pf_ctr_phi += 1; // updating counters of phi and csoil
                                                    if ( pf_ctr_phi > pf_nsteps_phi ) {
                                                        pf_ctr_phi = 1;
                                                        pf_ctr_csoil += 1;
                                                    }
                                                }
                                                else if ( PF == 1 ) {
                                                    // for slope failure probability with sampling of parameter combinations:

                                                    phi = ( (float)(pf_phi[j][ipar][j4]) ) / 10000; // defining phi

                                                    if ( pdftype_csoil == 4 ) { // defining csoil:
                                                        csoil = regpar1 * phi * 180 / PI + regpar2 + (float)pf_csoil[j][ipar][j4];
                                                        if ( csoil < 0 ) csoil = 0;
                                                    }
                                                    else csoil = (float)pf_csoil[j][ipar][j4];
                                                }
                                                
                                                if ( REVISED == 1 ) fact_revised = cos ( beta_m );
                                                else fact_revised = 1;

                                                if ( PSEUDOSTATIC == 1 ) { // for pseudostatic seismic slope stability model:

                                                    seisforce = weight_sub / 9.81 * seiscoef * epgas; //!!!CHECK dry or saturated density

                                                    r_cumcell = ( csoil * area + ( weight_sub * cos ( beta_c ) - seisforce * sin( beta_c ) + n_s ) * tan ( phi )) * fact_revised;
                                                        // shear resistance

                                                    t_cumcell = ( weight_sub * sin ( beta_m ) + seisforce * cos( beta_m ) + t_s ) * fact_revised; // shear force

                                                } else {

                                                    r_cumcell = ( csoil * area + ( weight_sub * cos ( beta_c ) + n_s ) * tan ( phi )) * fact_revised;
                                                        // shear resistance

                                                    t_cumcell = ( weight_sub * sin ( beta_m ) + t_s ) * fact_revised; // shear force
                                                }

                                                r_cum[j][pf_ctr] += r_cumcell;
                                                    // updating shear resistance for layer and, for pf, sampled values of csoil and phi
                                                t_cum[j][pf_ctr] += t_cumcell;
                                                    // updating shear force for layer and, for pf, sampled values of csoil and phi

                                            }
                                            if ( CLASSINF == 1 ) area_eff[j] = 9999; // updating control for minimum slip surface area
                                            else area_eff[j] += 1;
                                        }


                                    // Extracting data for profile file

                                        if ( ELLDENS == 0 && validity[j] == 1 ) { // if only one ellipsoid is tested and layer is valid:

                                            for ( s2_ctrl2 = 0; s2_ctrl2 <=s2_ctrl; s2_ctrl2++ ) { // loop over all profile coordinates:

                                                x = ex + xm - 1; // transforming pixel coordinates
                                                y = ey + ym - 1;

                                                if ( x == s2_x[s2_ctrl2] && y == s2_y[s2_ctrl2] ) { // if coordinates match:

                                                    if ( epds == UNDEFPLUS * CHANGEDM || epds == 0 ) s2_data[s2_ctrl2][j][1] = 0; 
                                                    else s2_data[s2_ctrl2][j][1] = ((float)(epes-epds)) / CHANGEDM;
                                                    s2_data[s2_ctrl2][j][2] = ((float)(epes-dslip_test[j])) / CHANGEDM;
                                                    if ( s2_data[s2_ctrl2][j][2] < depthslipmin ) depthslipmin = s2_data[s2_ctrl2][j][2];
                                                    if ( dslip_test[j] > 0 ) s2_data[s2_ctrl2][j][3] = n_s / 1000;
                                                    if ( dslip_test[j] > 0 ) s2_data[s2_ctrl2][j][4] = -t_s / 1000;
                                                    if ( dslip_test[j] > 0 ) s2_data[s2_ctrl2][j][5] = r_cumcell / 1000;
                                                    if ( dslip_test[j] > 0 ) s2_data[s2_ctrl2][j][6] = -t_cumcell / 1000;
                                                    s2_data[s2_ctrl2][j][7] = ((float)(epes-dgws)) / CHANGEDM;
                                                        // writing relevant parameters to profile array

                                                    break;
                                                }
                                            }
                                        }
                                    }
		                }
                            }
	                }
                    }


                // Stability of slip surface

                    fosmin = UNDEFPLUS; // initializing variables: minimum factor of safety
                    fosdeep = UNDEFPLUS; // deepest factor of safety lower than 1

                    if ( NEWMARK == 1 ) dnewmark = UNDEFMIN;
                    if ( PFDEPTH == 1 || ( PF == 1 && TRUNCATED == 0 )) {
                        pfail = 0; // slope failure probability
                        fos_mean = 0; // mean of fos

                        if ( NEWMARK == 1 ) {
                            dnewmark_mean = 0;
                            pnewmark = 0;
                        }
                    }

                    jcrit2 = UNDEFMIN; // layer with overall minimum factor of safety
                    k = 0; // layer with lowest factor of safety
                    l = 0; // deepest layer with factor of safety lower than 1

                    if ( TRUNCATED == 1 ) { // for truncated ellipsoids:

                        for ( j = 1; j <= jmax; j++ ) { // loop over all layers:

                            if ( validity[j] == 1 ) { // if layer is valid:

                                if ( PF == 1 && PFDEPTH == 0 ) {
                                    pfaill[j] = 0;
                                    fos_meanl[j] = 0; 
                                    // initializing variables for layer-specific slope failure probability and sum of fos

                                    if ( NEWMARK == 1 ) {
                                        dnewmark_meanl[j] = 0;
                                        pnewmarkl[j] = 0;
                                    }
                                }

                                for ( pf_ctr = 1; pf_ctr <= pf_ctr_max; pf_ctr++ ) { // loop over all samples (for slope failure probability)

                                    if ( r_cum[j][pf_ctr] > 0 && t_cum[j][pf_ctr] > 0 && area_eff[j] >= 10 ) { // if parameters are valid:

                                        fos_surf[j][pf_ctr] = r_cum[j][pf_ctr] / t_cum[j][pf_ctr]; // layer-specific factor of safety

                                        if ( PFDEPTH == 1 && PFSAMPLING <= 1 ) fos_mean += fos_surf[j][pf_ctr] / (float)pf_nsteps;
                                            // updating mean fos
                                        else if ( PFDEPTH == 1 && PFSAMPLING == 2 )
                                            fos_mean += fos_surf[j][pf_ctr] * pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                        else if ( PF == 1 && PFSAMPLING <= 1 ) fos_meanl[j] += fos_surf[j][pf_ctr] / (float)pf_nsteps;
                                        else if ( PF == 1 && PFSAMPLING == 2 ) fos_meanl[j] += fos_surf[j][pf_ctr] * pf_weight_geo[pf_ctr];

                                        if ( NEWMARK == 1 ) {

                                            gnewmark = fnewmark( fos_surf[j][pf_ctr], beta, magnitude, distance, eqdepth, pars1, pars2, newmarkcoef, newmarkref, iariasmin, iariasmax, 
                                                pparias, narias );
                                            if ( PF == 0 && PFDEPTH == 0 ) dnewmark_surf[j][pf_ctr] = gnewmark[0];
                                            else {
                                                dnewmark_surf[j][pf_ctr] = gnewmark[2];
                                                pnewmark_surf[j][pf_ctr] = gnewmark[1];
                                            }
                                            free( gnewmark );

                                            if ( PFDEPTH == 1 && PFSAMPLING <= 1 ) dnewmark_mean += dnewmark_surf[j][pf_ctr] / (float)pf_nsteps;
                                            else if ( PFDEPTH == 1 && PFSAMPLING == 2 )
                                                dnewmark_mean += dnewmark_surf[j][pf_ctr] * pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                            else if ( PF == 1 && PFSAMPLING <= 1 ) dnewmark_meanl[j] += dnewmark_surf[j][pf_ctr] / (float)pf_nsteps;
                                            else if ( PF == 1 && PFSAMPLING == 2 ) dnewmark_meanl[j] += dnewmark_surf[j][pf_ctr] * pf_weight_geo[pf_ctr];
                                        }
                                    }
                                    else {

                                        fos_surf[j][pf_ctr] = UNDEFPLUS;
                                        if ( NEWMARK == 1 ) dnewmark_surf[j][pf_ctr] = UNDEFMIN;
                                        validity[j] = 0;
                                    }

                                    if ( fos_surf[j][pf_ctr] < fosmin ) { // updating lowest factor of safety and associated layer
                                        fosmin = fos_surf[j][pf_ctr];
                                        if ( NEWMARK == 1 ) dnewmark = dnewmark_surf[j][pf_ctr];
                                        k = j;
                                    }

                                    if ( fos_surf[j][pf_ctr] < 1 ) { // updating deepest factor of safety lower than 1 and associated layer
                                        fosdeep = fos_surf[j][pf_ctr];
                                        l = j;

                                        if ( PFDEPTH == 1 && PFSAMPLING <= 1 ) pfail += 1 / (float)pf_nsteps;
                                            // updating slope failure probability
                                        else if ( PFDEPTH == 1 && PFSAMPLING == 2 ) pfail += pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                        else if ( PF == 1 && PFSAMPLING <= 1 ) pfaill[j] += 1 / (float)pf_nsteps;
                                        else if ( PF == 1 && PFSAMPLING == 2 ) pfaill[j] += pf_weight_geo[pf_ctr];
                                    }

                                    if ( NEWMARK == 1 ) {

                                        if ( PFDEPTH == 1 && PFSAMPLING <= 1 ) pnewmark += pnewmark_surf[j][pf_ctr] / (float)pf_nsteps;
                                        else if ( PFDEPTH == 1 && PFSAMPLING == 2 ) pnewmark += pnewmark_surf[j][pf_ctr] * pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                        else if ( PF == 1 && PFSAMPLING <= 1 ) pnewmarkl[j] += pnewmark_surf[j][pf_ctr] / (float)pf_nsteps;
                                        else if ( PF == 1 && PFSAMPLING == 2 ) pnewmarkl[j] += pnewmark_surf[j][pf_ctr] * pf_weight_geo[pf_ctr];
                                    }

                                    if ( PF == 0 && fos_surf[j][pf_ctr] < fos_critabs ) { // updating complementary parameters

                                        ell_centerx_crit = ell_centerx_metric;
                                        yell0_crit = ell_centery_metric;
                                        ella_crit = ell_a;
                                        ellb_crit = ell_b;

                                        ellc_crit = ell_c;
                                        zb_crit = ell_zb;
                                        fit_crit = fit;
                                        alpha_crit = alpha * 180 / PI;
                                        beta_crit = beta * 180 / PI;
                                        ell_width_crit = ell_width;
                                        ell_length_crit = ell_length;
                                        ell_area_crit = ell_area;
                                        fos_critabs = fos_surf[j][pf_ctr];
                                        jcrit = j;
                                        jcrit2 = j;
                                    }
                                }
                            }
                        }
                    }

                    else { // for non-truncated ellipsoids:

                        for ( pf_ctr = 1; pf_ctr <= pf_ctr_max; pf_ctr++ ) { // loop over all samples (for slope failure probability)

                            if ( r_cum[lbtm][pf_ctr] > 0 && t_cum[lbtm][pf_ctr] > 0 ) { // if parameters are valid:
                                fos_surf[lbtm][pf_ctr] = r_cum[lbtm][pf_ctr] / t_cum[lbtm][pf_ctr]; // factor of safety
                                fosmin = fos_surf[lbtm][pf_ctr]; // updating lowest factor of safety
                                fosdeep = fos_surf[lbtm][pf_ctr]; // updating deepest factor of safety lower than 1
                                k = lbtm; // setting associated layers to bottom layer (for integrity of code)
                                l = lbtm;

                                if ( NEWMARK == 1 ) {

                                    gnewmark = fnewmark( fos_surf[lbtm][pf_ctr], beta, magnitude, distance, eqdepth, pars1, pars2, newmarkcoef, newmarkref, iariasmin, iariasmax, 
                                        pparias, narias );
                                    if ( PF == 0 && PFDEPTH == 0 ) dnewmark_surf[lbtm][pf_ctr] = gnewmark[0];
                                    else {
                                        dnewmark_surf[lbtm][pf_ctr] = gnewmark[2];
                                        pnewmark_surf[lbtm][pf_ctr] = gnewmark[1];
                                    }
                                    free( gnewmark );
                                }
                            }

                            if ( PF == 1 && PFSAMPLING <= 1 )
                                fos_mean += fos_surf[lbtm][pf_ctr] / (float)pf_nsteps; // updating mean fos
                            else if ( PF == 1 && PFSAMPLING == 2 ) fos_mean += fos_surf[lbtm][pf_ctr] * pf_weight_geo[pf_ctr];


                            if ( NEWMARK == 1 ) {

                                if ( PF == 1 && PFSAMPLING <= 1 ) dnewmark_mean += dnewmark_surf[lbtm][pf_ctr] / (float)pf_nsteps;
                                else if ( PF == 1 && PFSAMPLING == 2 )
                                    dnewmark_mean += dnewmark_surf[lbtm][pf_ctr] * pf_weight_geo[pf_ctr];
                            }

                            if ( fos_surf[lbtm][pf_ctr] < 1 && PF == 1 && PFSAMPLING <= 1 )
                                pfail += 1 / (float)pf_nsteps; // updating slope failure probability
                            else if ( fos_surf[lbtm][pf_ctr] < 1 && PF == 1 && PFSAMPLING == 2 ) pfail += pf_weight_geo[pf_ctr];

                            if ( NEWMARK == 1 ) {

                                if ( PF == 1 && PFSAMPLING <= 1 )
                                    pnewmark += pnewmark_surf[lbtm][pf_ctr] / (float)pf_nsteps;
                                else if ( PF == 1 && PFSAMPLING == 2 ) pnewmark += pnewmark_surf[lbtm][pf_ctr] * pf_weight_geo[pf_ctr];

                            } 

                            if ( PF == 0 && fos_surf[lbtm][pf_ctr] < fos_critabs ) { // updating complementary parameters
                                ell_centerx_crit = ell_centerx_metric;
                                yell0_crit = ell_centery_metric;
                                ella_crit = ell_a;
                                ellb_crit = ell_b;
                                ellc_crit = ell_c;
                                zb_crit = ell_zb;
                                fit_crit = fit;
                                alpha_crit = alpha * 180 / PI;
                                beta_crit = beta * 180 / PI;
                                ell_width_crit = ell_width;
                                ell_length_crit = ell_length;
                                ell_area_crit = ell_area;
                                fos_critabs = fos_surf[lbtm][pf_ctr];
                                jcrit = lbtm;
                                jcrit2 = lbtm;
                            }
                        }
                    }


                // Computing standard deviation of fos

                    if ( PF == 1 ) { // for slope failure probability mode:

                       fos_stdev = 0;
                       if ( NEWMARK == 1 ) dnewmark_stdev = 0;
                       if ( TRUNCATED == 1 ) { // for truncated ellipsoids:

                           for ( j = 1; j <= jmax; j++ ) { // loop over all layers:

                                if ( validity[j] == 1 ) { // if layer is valid:

                                    if ( PFDEPTH == 0 ) fos_stdevl[j] = 0; // initializing standard deviation of fos
                                    if ( PFDEPTH == 0 && NEWMARK == 1 ) dnewmark_stdevl[j] = 0;

                                    for ( pf_ctr = 1; pf_ctr <= pf_ctr_max; pf_ctr++ ) { // loop over all samples

                                        if ( PFDEPTH == 1 && PFSAMPLING <= 1 )
                                            fos_stdev += pow( fos_surf[j][pf_ctr] - fos_mean, 2 ) / (float)pf_nsteps;
                                            // updating standard deviation of fos
                                        else if ( PFDEPTH == 1 && PFSAMPLING == 2 )
                                            fos_stdev += pow( fos_surf[j][pf_ctr] - fos_mean, 2 )
                                                * pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                        else if ( PFSAMPLING <= 1 )
                                            fos_stdevl[j] += pow( fos_surf[j][pf_ctr] - fos_meanl[j], 2 ) / (float)pf_nsteps;
                                        else if ( PFSAMPLING == 2 )
                                            fos_stdevl[j] += pow( fos_surf[j][pf_ctr] - fos_meanl[j], 2 ) * pf_weight_geo[pf_ctr];

                                        if ( NEWMARK == 1 ) {

                                            if ( PFDEPTH == 1 && PFSAMPLING <= 1 )
                                                dnewmark_stdev += pow( dnewmark_surf[j][pf_ctr] - dnewmark_mean, 2 ) / (float)pf_nsteps;
                                            else if ( PFDEPTH == 1 && PFSAMPLING == 2 )
                                                dnewmark_stdev += pow( dnewmark_surf[j][pf_ctr] - dnewmark_mean, 2 )
                                                    * pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                            else if ( PFSAMPLING <= 1 )
                                                dnewmark_stdevl[j] += pow( dnewmark_surf[j][pf_ctr] - dnewmark_meanl[j], 2 ) / (float)pf_nsteps;
                                            else if ( PFSAMPLING == 2 )
                                                dnewmark_stdevl[j] += pow( dnewmark_surf[j][pf_ctr] - dnewmark_meanl[j], 2 ) * pf_weight_geo[pf_ctr];
                                        }
                                    }
                                }
                            }
                        }

                        else { // for non-truncated ellipsoids:

                            for ( pf_ctr = 1; pf_ctr <= pf_ctr_max; pf_ctr++ ) { // loop over all samples (for slope failure probability)

                                if ( PF == 1 && PFSAMPLING <= 1 )
                                    fos_stdev += pow( fos_surf[lbtm][pf_ctr] - fos_mean, 2 ) / (float)pf_nsteps;
                                        // updating standard deviation of fos
                                else if ( PF == 1 && PFSAMPLING == 2 )
                                    fos_stdev += pow( fos_surf[lbtm][pf_ctr] - fos_mean, 2 ) * pf_weight_geo[pf_ctr];

                                if ( NEWMARK == 1 ) {

                                    if ( PF == 1 && PFSAMPLING <= 1 )
                                        dnewmark_stdev += pow( dnewmark_surf[lbtm][pf_ctr] - dnewmark_mean, 2 ) / (float)pf_nsteps;
                                            // updating standard deviation of fos
                                    else if ( PF == 1 && PFSAMPLING == 2 )
                                        dnewmark_stdev += pow( dnewmark_surf[lbtm][pf_ctr] - dnewmark_mean, 2 ) * pf_weight_geo[pf_ctr];
                                }
                            }
                        }
                    }
                }


            // Storing results

                if ( fosmin < 0.001 || fosmin == UNDEFPLUS ) ctrl_continue = 0; // updating validity of ellipsoid

                if ( ctrl_continue == 1 ) { // if ellipsoid is still valid:

                    ell_depth = 0; // initializing control for maximum depth of ellipsoid
                    for ( x = xm; x < xp; x++ ) { // loop over all pixels of ellipsoid mask:
	                for ( y = ym; y < yp; y++ ) {

                            ex = x - xm + 1; // transforming coordinates
                            ey = y - ym + 1;

                            Segment_get(&seg_depthfail_deep, &depthfail_deep, x, y); // reading depths
                            if ( USESEG == 1 ) Segment_get(&seg_epd[lbtm+1], &epd1s, ex, ey);
                            else epd1s = epd[lbtm+1][ex][ey];

                            if ( epd1s > ell_depth ) ell_depth = epd1s; // updating maximum depth of ellipsoid

                            if ( epd1s > 0 ) { // if depth of ellipsoid bottom is larger than zero:

                                Segment_get(&seg_fos_min, &fos_min, x, y); // reading old minimum factor of safety
                                if ( PF == 1 ) Segment_get(&seg_pfail, &pfail0, x, y); // reading old slope failure probability

                                if ( NEWMARK == 1 ) {

                                    Segment_get(&seg_d_newmark, &d_newmark, x, y);
                                    if ( PF == 1 ) Segment_get(&seg_pnewmark, &pnewmark0, x, y);

                                }

                                if ( ADDITIONAL == 1 ) { // for additional output, updating number of tested ellipsoids for pixel 
                                    Segment_get(&seg_cumsurf, &cumsurf, x, y);
                                    cumsurf += 1;
                                    Segment_put(&seg_cumsurf, &cumsurf, x, y);
                                }

                                if ( TRUNCATED == 1 && LAYERS == 1 ) { // for layer mode with truncated ellipsoids:

                                    if ( PF == 0 ) cumfos_critmax = 0; // initializing control for factor of safety lower than 1

                                    for ( j = 1; j <= jmax; j++ ) { // loop over all layers:

                                        output_add = 0; // initializing control for additional output

                                        if ( validity[j] == 1 ) { // if layer is valid:

                                            if ( PFDEPTH == 1 ) epd2s = pf_depth[j]; // reading depth of layer
                                            else if ( USESEG == 1 ) Segment_get(&seg_epd[j], &epd2s, ex, ey);
                                            else epd2s = epd[j][ex][ey];

                                            for ( pf_ctr = 1; pf_ctr <= pf_ctr_max; pf_ctr++ ) {
                                                // loop over all samples (relevant for slope failure probability):

                                                if ( ADDITIONAL == 1 && epd2s > 0 ) {
                                                    // for additional output if depth of layer is larger than zero:

                                                    if ( fos_surf[j][pf_ctr] < 1.0 ) cumfos_critmax = 1;
                                                        // updating control for factor of safety lower than 1
                                                }

                                                if ( PF == 0 && fos_surf[j][pf_ctr] < fos_min && epd2s > 0 ) {
                                                    // if factor of safety is lower than old value:
                                                    fos_min = fos_surf[j][pf_ctr]; // updating and storing lowest factor of safety for pixel
                                                    Segment_put(&seg_fos_min, &fos_min, x, y);
                                                    output_add = 1; // setting control for additional output to positive
                                                }

                                                if ( NEWMARK == 1 ) {

                                                    if ( PF == 0 && dnewmark_surf[j][pf_ctr] > d_newmark && epd2s > 0 ) {

                                                        d_newmark = dnewmark_surf[j][pf_ctr];
                                                        Segment_put(&seg_d_newmark, &d_newmark, x, y);
                                                        Segment_flush(&seg_d_newmark);
                                                    }
                                                }

                                                if ( fos_surf[j][pf_ctr] < 1 && epd2s > 0 ) {
                                                    // if factor of safety is lower than 1, additional output:
                                                    if ( epd2s >= epd1s ) depth_control = ((float)(epd1s))/CHANGEDM;
                                                    else depth_control = ((float)(epd2s))/CHANGEDM;

		                                    if ( depth_control > depthfail_deep ) {
                                                        fos_deep = fos_surf[j][pf_ctr];
                                                        depthfail_deep = depth_control;
                                                        Segment_put(&seg_fos_deep, &fos_deep, x, y);
                                                        Segment_put(&seg_depthfail_deep, &depthfail_deep, x, y);
                                                    }
                                                }
                                            }
                                            if ( PF == 1 && pfaill[j] > pfail0 && epd2s > 0 ) {
                                                // storing layer-specific slope failure probability:
                                                Segment_put(&seg_pfail, &pfaill[j], x, y);
                                                pfail0 = pfaill[j]; // updating maximum slope failure probability for pixel
                                            }

                                            if ( NEWMARK == 1 ) {

                                                if ( PF == 1 && pnewmarkl[j] > pnewmark0 && epd2s > 0 ) {

                                                    Segment_put(&seg_pnewmark, &pnewmarkl[j], x, y);
                                                    Segment_flush(&seg_pnewmark);
                                                    pnewmark0 = pnewmarkl[j];
                                                }
                                            }

                                            if ( PF == 1 && fos_meanl[j] < fos_min && epd2s > 0 ) {
                                                // if mean factor of safety is lower than old value:
                                                fos_min = fos_meanl[j]; // updating and storing lowest factor of safety for pixel
                                                fos_stdev0 = pow ( fos_stdevl[j] * pf_ctr_max / ( pf_ctr_max - 1 ), 0.5 );
                                                    // updating standard deviation of fos
                                                Segment_put(&seg_fos_min, &fos_min, x, y);
                                                Segment_put(&seg_fos_stdev, &fos_stdev0, x, y);
                                                output_add = 1; // setting control for additional output to positive

                                            }

                                            if ( NEWMARK == 1 ) {

                                                if ( PF == 1 && dnewmark_meanl[j] > d_newmark && epd2s > 0 ) {

                                                    d_newmark = dnewmark_meanl[j]; // updating and storing lowest factor of safety for pixel
                                                    dnewmark_stdev0 = pow ( dnewmark_stdevl[j]  * pf_ctr_max / ( pf_ctr_max - 1 ), 0.5 );
                                                    Segment_put(&seg_d_newmark, &d_newmark, x, y);
                                                    Segment_put(&seg_d_newmark_stdev, &dnewmark_stdev0, x, y);
                                                    Segment_flush(&seg_d_newmark);
                                                    Segment_flush(&seg_d_newmark_stdev);
                                                }
                                            }

                                            if ( output_add == 1 && ADDITIONAL == 1 ) {
                                                // for additional output, updating and storing results

                                                if ( ELLDENS == 0 ) {
                                                    // for additional output or only one ellipsoid tested, updating and storing results
                                                    if ( epd2s >= epd1s ) depthfail_min = ((float)(epd1s))/CHANGEDM;
                                                    else depthfail_min = ((float)(epd2s))/CHANGEDM;
                                                    Segment_put(&seg_depthfail_min, &depthfail_min, x, y);
                                                }
                                                lyr_crit = ldata_trans[j];
                                                amin = ell_length;
                                                bmin = ell_width;
                                                zbmin = ell_zb / ell_c;
                                                Segment_put(&seg_lyr_crit, &lyr_crit, x, y);
                                                Segment_put(&seg_amin, &amin, x, y);
                                                Segment_put(&seg_bmin, &bmin, x, y);
                                                Segment_put(&seg_zbmin, &zbmin, x, y);
                                            }

                                            if ( ADDITIONAL == 1 ) {
                                                // for additional output, updating and storing number of tested ellipsoids
                                                // with factor of safety lower than one for pixel:
                                                Segment_get(&seg_cumfos_crit, &cumfos_crit, x, y);
                                                cumfos_crit = cumfos_crit + cumfos_critmax;
                                                Segment_put(&seg_cumfos_crit, &cumfos_crit, x, y);
                                            }
                                        }
                                    }
                                }
                                else { // for soil class mode or mode with non-truncated ellipsoids:

                                    output_add = 0; // initializing control for additional output

                                    if ( PFDEPTH == 1 ) epd2s = pf_depth[k]; // reading depth of most unstable layer
                                    else if ( USESEG == 1 ) Segment_get(&seg_epd[k], &epd2s, ex, ey);
                                    else epd2s = epd[k][ex][ey];

                                    if ( ADDITIONAL == 1 && epd2s > 0 && fosmin < 1.0 ) {
                                        // for additional output if depth of layer is larger than zero and factor of safety is lower than 1:
                                        Segment_get(&seg_cumfos_crit, &cumfos_crit, x, y);
                                        cumfos_crit += 1;
                                        Segment_put(&seg_cumfos_crit, &cumfos_crit, x, y);
                                            // updating and storing control for factor of safety lower than 1
                                    }

                                    if ( PF == 0 && fosmin < fos_min && epd2s > 0 ) { // if factor of safety is lower than old value:
                                        fos_min = fosmin; // updating and storing lowest factor of safety for pixel
                                        Segment_put(&seg_fos_min, &fos_min, x, y);
                                        output_add = 1; // setting control for additional output to positive
                                    }

                                    if ( NEWMARK == 1 ) {

                                        if ( PF == 0 && dnewmark > d_newmark && epd2s > 0 ) {

                                            d_newmark = dnewmark;
                                            Segment_put(&seg_d_newmark, &d_newmark, x, y);
                                            Segment_flush(&seg_d_newmark);
                                        }
                                    }

                                    if ( fosdeep < 1 ) { // if additional output:

                                        if ( USESEG == 1 ) Segment_get(&seg_epd[l], &epd2s, ex, ey);
                                            // reading depth of deepest layer with safety factor lower than 1
                                        else epd2s = epd[l][ex][ey];
                                        if ( epd2s >= epd1s ) depth_control = ((float)(epd1s))/CHANGEDM; // defining depth of slip surface
                                        else depth_control = ((float)(epd2s))/CHANGEDM;

		                        if ( epd2s > 0 && depth_control > depthfail_deep ) {
                                            fos_deep = fosdeep;
                                            depthfail_deep = depth_control;
                                            Segment_put(&seg_fos_deep, &fos_deep, x, y);
                                            Segment_put(&seg_depthfail_deep, &depthfail_deep, x, y);
                                        }
                                    }

                                    if ( PF == 1 && ( PFDEPTH == 1 || TRUNCATED == 0 ) ) {
                                        // updating and storing slope failure probability:
                                        if ( pfail > pfail0 && epd2s > 0 ) {

                                            Segment_put(&seg_pfail, &pfail, x, y);

                                        }

                                        if ( NEWMARK == 1 ) {

                                            if ( pnewmark > pnewmark0 && epd2s > 0 ) {

                                                Segment_put(&seg_pnewmark, &pnewmark, x, y);
                                                Segment_flush(&seg_pnewmark);
                                            }
                                        }

                                        if ( fos_mean < fos_min && epd2s > 0 ) {
                                            // if mean factor of safety is lower than old value:
                                            fos_min = fos_mean; // updating and storing lowest factor of safety for pixel
                                            fos_stdev0 = pow ( fos_stdev  * pf_ctr_max / ( pf_ctr_max - 1 ), 0.5 );
                                                // updating standard deviation of fos
                                            Segment_put(&seg_fos_min, &fos_min, x, y);
                                            Segment_put(&seg_fos_stdev, &fos_stdev0, x, y);
                                            output_add = 1; // setting control for additional output to positive
                                        }

                                        if ( NEWMARK == 1 ) {

                                            if ( dnewmark_mean > d_newmark && epd2s > 0 ) {

                                                d_newmark = dnewmark_mean; // updating and storing lowest factor of safety for pixel
                                                dnewmark_stdev0 = pow ( dnewmark_stdev  * pf_ctr_max / ( pf_ctr_max - 1 ), 0.5 );
                                                Segment_put(&seg_d_newmark, &d_newmark, x, y);
                                                Segment_put(&seg_d_newmark_stdev, &dnewmark_stdev0, x, y);
                                                Segment_flush(&seg_d_newmark);
                                                Segment_flush(&seg_d_newmark_stdev);
                                            }
                                        }
                                    }
                                    else if ( PF == 1 ) { // updating and storing layer-specific slope failure probability:

                                        for ( j = 1; j <= jmax; j++ ) {

                                            if ( pfaill[j] > pfail0 && epd2s > 0 && validity[j] == 1 ) {
                                                Segment_put(&seg_pfail, &pfaill[j], x, y);
                                                pfail0 = pfaill[j];
                                            }

                                            if ( NEWMARK == 1 ) {

                                                if ( pnewmarkl[j] > pnewmark0 && epd2s > 0 && validity[j] == 1 ) {

                                                    Segment_put(&seg_pnewmark, &pnewmarkl[j], x, y);
                                                    Segment_flush(&seg_pnewmark);
                                                    pnewmark0 = pnewmarkl[j];
                                                }
                                            }

                                            if ( fos_meanl[j] < fos_min && epd2s > 0 && validity[j] == 1 ) {
                                                // if mean factor of safety is lower than old value:
                                                fos_min = fos_meanl[j]; // updating and storing lowest factor of safety for pixel
                                                fos_stdev0 = pow ( fos_stdevl[j]  * pf_ctr_max / ( pf_ctr_max - 1 ), 0.5 );
                                                    // updating standard deviation of fos
                                                Segment_put(&seg_fos_min, &fos_min, x, y);
                                                Segment_put(&seg_fos_stdev, &fos_stdev0, x, y);
                                                output_add = 1; // setting control for additional output to positive
                                            }

                                            if ( NEWMARK == 1 ) {

                                                if ( dnewmark_meanl[j] > d_newmark && epd2s > 0 && validity[j] == 1 ) {

                                                     d_newmark = dnewmark_meanl[j]; // updating and storing lowest factor of safety for pixel
                                                     dnewmark_stdev0 = pow ( dnewmark_stdevl[j]  * pf_ctr_max / ( pf_ctr_max - 1 ), 0.5 );
                                                     Segment_put(&seg_d_newmark, &d_newmark, x, y);
                                                     Segment_put(&seg_d_newmark_stdev, &dnewmark_stdev0, x, y);
                                                     Segment_flush(&seg_d_newmark);
                                                     Segment_flush(&seg_d_newmark_stdev);
                                                }
                                            }
                                        }
                                    }

                                    if ( output_add == 1 && ADDITIONAL == 1 ) { // for additional output, updating and storing results

                                        if ( ELLDENS == 0 ) {
                                            // for additional output or only one ellipsoid tested, updating and storing results
                                            if ( epd2s >= epd1s ) depthfail_min = ((float)(epd1s))/CHANGEDM;
                                            else depthfail_min = ((float)(epd2s))/CHANGEDM;
                                            Segment_put(&seg_depthfail_min, &depthfail_min, x, y);
                                        }
                                        lyr_crit = ldata_trans[k];
                                        amin = ell_length;
                                        bmin = ell_width;
                                        zbmin = ell_zb / ell_c;
                                        Segment_put(&seg_lyr_crit, &lyr_crit, x, y);
                                        Segment_put(&seg_amin, &amin, x, y);
                                        Segment_put(&seg_bmin, &bmin, x, y);
                                        Segment_put(&seg_zbmin, &zbmin, x, y);
		                    }
                                }
		            }
                        }
                    }
                    ell_depth = ell_depth / CHANGEDM; // maximum depth of ellipsoid


                // Storing slip surface with absolutely lowest factor of safety

                    if ( jcrit2 > UNDEFMIN && ADDITIONAL == 1 ) { // if ellipsoid shows the absolutely lowest factor of safety:
                        for ( x = 0; x < m; x++ ) { // loop over ellipsoid mask:
	                    for ( y = 0; y < n; y++ ) {
                                depthfail_crit = UNDEFMIN; // resetting depth of slip surface with absolutely lowest factor of safety
                                Segment_put(&seg_depthfail_crit, &depthfail_crit, x, y);
                            }
	                }
                        Segment_flush(&seg_depthfail_crit);

                        for ( x = xm; x < xp; x++ ) { // another loop over the ellipsoid mask:
	                    for ( y = ym; y < yp; y++ ) {

                                ex = x - xm + 1; // transforming coordinates
                                ey = y - ym + 1;

                                ell_depth_crit = ell_depth; // updating maximum depth of ellipsoid with lowest factor of safety
                                if ( USESEG == 1 ) Segment_get(&seg_epd[lbtm+1], &epd1s, ex, ey); // reading depth of ellipsoid bottom
                                else epd1s = epd[lbtm+1][ex][ey];

                                if ( USESEG == 1 ) Segment_get(&seg_epd[jcrit2], &epd2s, ex, ey); // reading depth of layer
                                else epd2s = epd[jcrit2][ex][ey];

                                if ( epd2s >= epd1s ) depthfail_crit = ((float)(epd1s))/CHANGEDM; // defining depth of slip surface
                                else depthfail_crit = ((float)(epd2s))/CHANGEDM;
                                Segment_put(&seg_depthfail_crit, &depthfail_crit, x, y); // storing result
                            }
	                }
                    }


                // Forcing to write results to segment files

                    Segment_flush(&seg_fos_min);
                    if ( PF == 1 ) {
                        Segment_flush(&seg_pfail);
                        Segment_flush(&seg_fos_stdev);
                    }
                    Segment_flush(&seg_depthfail_deep);
                    Segment_flush(&seg_fos_deep);
                    if ( ADDITIONAL == 1 ) {
                        Segment_flush(&seg_cumfos_crit);
                        Segment_flush(&seg_cumsurf);
                        Segment_flush(&seg_depthfail_min);
                        Segment_flush(&seg_lyr_crit);
                        Segment_flush(&seg_depthfail_crit);
                        Segment_flush(&seg_amin);
                        Segment_flush(&seg_bmin);
                        Segment_flush(&seg_zbmin);
                    }
                    if ( ELLDENS == 0 ) Segment_flush(&seg_depthfail_min);

                    if ( ctrl_continue == 1 ) ctr_valid += 1; // updating number of valid ellipsoids

                    time_stepend = clock(); // end time of processing ellipsoid
                    time_elapsed = ( ( float ) ( time_stepend - time_stepstart ) ) / CLOCKS_PER_SEC; // time needed for processing ellipsoid


                // Updating summary file

                    if ( SUMMARY == 1 ) { // if mode includes summary, writing summary of ellipsoid to file:

                        fprintf ( fsummary1, "%i\t%i\t%.0f\t%.0f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f\t%.0f\t%.2f\t%.3f",
                            ctr_valid, ctr_mcsim, ell_centerx_metric, ell_centery_metric, ell_a, ell_b, ell_c, ell_zb,
                            fit, alpha * 180 / PI, beta * 180 / PI, ell_length, ell_width, ell_area, ell_depth, fosmin );

                        if ( TRUNCATED == 1 ) for ( j = 1; j <= nly; j++ ) fprintf ( fsummary1, "\t%.3f", fos_surf[j][1] );
                        fprintf ( fsummary1, "\t%.2f\n", time_elapsed );
                    }


                // Updating profile file

                    if ( ELLDENS == 0 ) { // if only one ellipsoid is tested, writing data to profile file:

                        for ( k=1; k<=7; k++ ) {
                            for ( j=1; j <= lbtm; j++ ) {
                                if ( k == 1 ) {
                                    for ( s2_ctrl2 = 0; s2_ctrl2 <=s2_ctrl; s2_ctrl2++ ) {
                                        fprintf ( fprofile, "%.2f\t", s2_data[s2_ctrl2][0][0]);
                                    }
                                    fprintf ( fprofile, "\n");
                                }
                            }
                            for ( j=1; j <= lbtm; j++ ) {
                                if ( k == 1 ) {
                                    for ( s2_ctrl2 = 0; s2_ctrl2 <=s2_ctrl; s2_ctrl2++ ) {
                                        if ( s2_data[s2_ctrl2][0][1] == 0 ) fprintf ( fprofile, "NA\t");
                                        else fprintf ( fprofile, "%.2f\t", s2_data[s2_ctrl2][0][1]);
                                    }
                                    fprintf ( fprofile, "\n");
                                }
                            }
                            for ( j=1; j <= lbtm; j++ ) {
                                for ( s2_ctrl2 = 0; s2_ctrl2 <=s2_ctrl; s2_ctrl2++ ) {

                                    if ( k <= 2 && ( s2_data[s2_ctrl2][j][k] == 0 || s2_data[s2_ctrl2][j][k] < depthslipmin - 5 ))
                                        fprintf ( fprofile, "NA\t");
                                    else if ( k > 2 && s2_data[s2_ctrl2][j][k] == 0 && j==jcrit) fprintf ( fprofile, "NA\t");
                                    else if ( k <= 2 || j==jcrit ) fprintf ( fprofile, "%.2f\t", s2_data[s2_ctrl2][j][k]);
                                }
                                if ( k <= 2 || j==jcrit ) fprintf ( fprofile, "\n");
                            }
                        }
                    }
                }
            }
            interval_ctr += 1; // updating counter for number of tested ellipsoids since start or latest output


        // Updating evolution file, counters and terminal output at an interval of ellipsoid density of 1

            if ( ELLDENS > 0 && interval_ctr == elldens_interval ) { // writing evolution of results to file:
                cell_change = 0;
                for ( x = 1; x < m-1; x++ ) { // loop over all pixels:
	            for ( y = 1; y < n-1; y++ ) {

                        if ( PF == 0 ) { // for factor of safety:
                            Segment_get(&seg_fos_min, &fos_min, x, y);
                            if ( fos_min < 1 ) cell_change += 1;
                        }
                        else { // for slope failure probability:

                            if ( NEWMARK == 1 ) {

                                Segment_get(&seg_pnewmark, &pnewmark, x, y);
                                if ( pnewmark > 0 ) cell_change += pnewmark;

                            }

                            else {

                                Segment_get(&seg_pfail, &pfail, x, y);
                                if ( pfail > 0 ) cell_change += pfail;

                            }
                        }
                    }
                }
                num_elldens_interval += 1; // updating tested ellipsoid density
                fprintf (fchange, "%i\t%i\t%.3f\n", num_elldens_interval, ctr_mcsim, cell_change); // writing data to evolution file
                interval_ctr = 0; // resetting counter for number of tested ellipsoids since start or latest output

                pccomp = 100 * ((float)ctr_mcsim/(float)NUM_MCSIM); // percentage of already tested ellipsoids out of all ellipsoids
                if ( PF == 0 ) printf("Computing factor of safety: %.0f%% of %i ellipsoids completed ...\r", pccomp, NUM_MCSIM);
                else printf("Computing slope failure probability: %.0f%% of %i ellipsoids completed ...\r", pccomp, NUM_MCSIM);
                fflush(stdout);
            }
        }


    // Completing summary files and updating terminal output

        if ( SUMMARY == 1 && ctrl_content == 1 ) { // if mode includes summary:

            fprintf ( fsummary1, "\nCRIT\t\t%.0f\t%.0f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f\t%.0f\t%.2f\t%.3f",
                ell_centerx_crit, yell0_crit, ella_crit, ellb_crit, ellc_crit, zb_crit,
                fit_crit, alpha_crit, beta_crit, ell_length_crit, ell_width_crit, ell_area_crit, ell_depth_crit, fos_critabs );

            if ( TRUNCATED == 1 ) { //!!!CHECK brackets
            
                for ( j = 1; j <= nly; j++ ) {
                
                    if ( j == jcrit ) fprintf ( fsummary1, "\tCRIT" ); 
                    else fprintf ( fsummary1, "\t----" );
                    
                }
            }
            fprintf ( fsummary1, "\n" ); // parameters of critical slip surface
        }

        fprintf ( fsummary2, "%i\n", lbtm_max ); // writing maximum number of layers in local reference system to second summary file
        fprintf ( fsummary2, "%i\n", nly ); // writing maximum allowed number of layers in local reference system to second summary file

        if ( ELLDENS == 0 ) { // if profile plot shall be drawn, write parameters necessary for graphic design to second summary file:

            fprintf ( fsummary2, "%i\n", lbtm );
            fprintf ( fsummary2, "%i\n", jcrit );
            fprintf ( fsummary2, "%.3f\n", fosmin );
        }

        time_end = clock(); // time at end of calculation
        time_elapsed = ( ( float ) ( time_end - time_start ) ) / CLOCKS_PER_SEC; // time needed for calculation

        if ( PF == 0 ) printf("Computing factor of safety ... completed in %.2f seconds.                                    \n",
            time_elapsed);
        else printf("Computing slope failure probability ... completed in %.2f seconds.                                    \n",
            time_elapsed);
        printf("Number of valid ellipsoids: %i out of %i.\n", ctr_valid, ctr_mcsim );
    }


// Computing factor of safety with infinite slope stability model

    if ( INFMOD == 1 ) { // for infinite slope stability mode:

        time_start = clock(); // start time

        if ( PF == 0 ) printf("Computing factor of safety for infinite slope ...\r" );
        else printf("Computing slope failure probability for infinite slope ...\r" );

        for ( x = 1; x < m-1; x++ ) { // loop over all pixels:
	    for ( y = 1; y < n-1; y++ ) {

                ctrl_continue = 1; // initializing control for validity of pixel
                xnh[0]=x; xnh[1]=x; xnh[2]=x+1; xnh[3]=x+1; xnh[4]=x+1; xnh[5]=x; xnh[6]=x-1; xnh[7]=x-1; xnh[8]=x-1;
                ynh[0]=y; ynh[1]=y-1; ynh[2]=y-1; ynh[3]=y; ynh[4]=y+1; ynh[5]=y+1; ynh[6]=y+1; ynh[7]=y; ynh[8]=y-1;
                    // defining neighbourhood of pixel

                for ( q=0; q<9; q++ ) { // testing whether neighbourhood of cell is valid:
                    Segment_get(&seg_soilclass, &pscs, xnh[q], ynh[q]);

                    if ( pscs == UNDEFMIN ) { // testing validity of pixel
                        ctrl_continue = 0;
                        break;
                    }
                }

                depthfail_deep=-9999;
                Segment_put(&seg_depthfail_deep, &depthfail_deep, x, y);

                if ( ctrl_continue == 1 ) { // if pixel is still valid, initializing variables

                    Segment_get(&seg_soilclass, &i, x, y);
                    j=1;
                    weight_sub=0;
                    u_appr=0;
                    fos_min=9998;
                    fos_deep=9998;
                    depthfail_min=-9999;
                    if ( NEWMARK == 1 ) d_newmark = 0;

                    if ( PF == 1 ) { 

                        pfail = 0; // initializing slope failure probability
                        if ( PFDEPTH == 0 ) fos_mean = 9998; else fos_mean = 0; 
                        if ( NEWMARK == 1 ) { pnewmark = 0; dnewmark_mean = 0; }

                    }

                    for ( q=0; q<9; q++ ) {

                        Segment_get(&seg_elev, &ise[q], xnh[q], ynh[q]); // reading elevation of neighbourhood
                        Segment_get(&seg_gwdepth, &isg[q], xnh[q], ynh[q]); // reading depth of groundwater table
                    }

                    while ( j <= jmax ) { // loop over all layers:

                        if ( PFDEPTH == 1 ) j2 = 1;
                        else j2 = j;


                // Slope of the layer bottom

                        for ( q=0; q<9; q++ ) {

                            if ( PFDEPTH == 1 ) isdn[q] = ( (float)pf_depth[j] ) / CHANGEDM;
                            else Segment_get(&seg_depth[j], &isdn[q], xnh[q], ynh[q]);
                            ienh[q] = ise[q];// - isdn[q]; !!!Surface slope is applied
                        }

                        dz_dx = (( ienh[4] + 2 * ienh[3] + ienh[2] ) - ( ienh[6] + 2 * ienh[7] + ienh[8] )) / ( 8 * csz );
                        dz_dy = (( ienh[8] + 2 * ienh[1] + ienh[2] ) - ( ienh[6] + 2 * ienh[5] + ienh[4] )) / ( 8 * csz );
                        beta_xz = atan ( dz_dx );
                        beta_yz = atan ( dz_dy );
                        slope = atan ( pow ( pow ( dz_dx, 2 ) + pow ( dz_dy, 2 ) , 0.5 ) );
                        if ( dz_dx == 0 && dz_dy == 0 ) slope = 0;


                    // Slope of the groundwater table

                        for ( q=0; q<9; q++ ) ignh[q] = ise[q] - isg[q];

                        dz_dx = (( ignh[4] + 2 * ignh[3] + ignh[2] ) - ( ignh[6] + 2 * ignh[7] + ignh[8] )) / ( 8 * csz );
                        dz_dy = (( ignh[8] + 2 * ignh[1] + ignh[2] ) - ( ignh[6] + 2 * ignh[5] + ignh[4] )) / ( 8 * csz );
                        gslope = atan ( pow ( pow ( dz_dx, 2 ) + pow ( dz_dy, 2 ) , 0.5 ) );
                        if ( dz_dx == 0 && dz_dy == 0 ) gslope = 0;
                        gslope = slope; //!!!Surface terrain slope is applied

                        if ( PFDEPTH == 1 ) {
                            isd[j] = ( (float)pf_depth[j] ) / CHANGEDM;
                            isd11 = ( (float)pf_depth[j-1] ) / CHANGEDM;
                        }
                        else {
                            Segment_get(&seg_depth[j], &isd[j], x, y);
                            Segment_get(&seg_depth[j-1], &isd11, x, y);
                        }

                        if ( isg[0] >= isd[j] ) depth_layer = 0; // depth of saturated layer
                        else if ( isg[0] > isd11 ) depth_layer = isd[j] - isg[0];
                        else depth_layer = isd[j] - isd11;


                    // Soil parameters and pore water pressure related to slip surface

                        weight_sub += pgamma[i][j2] * ( isd[j] - isd11 ) * pow ( csz , 2 );
                        if ( ptheta[i][j2] > 0 && depth_layer > 0 ) {
                            weight_sub += GAMMA_WATER * ( ptheta[i][j2] - 1 ) * ( depth_layer ) * pow ( csz , 2 );
                            u_appr += GAMMA_WATER * ( depth_layer ) * pow ( csz , 2 ); 
                        }
                        else u_appr = 0;


                    // Stabilizing and destabilizing forces

                        if ( slope > 0 && isd[j] > 0 ) {

                            if ( PF == 0 ) { // for factor of safety:
                                phi=pphi[i][j2]; // defining phi and csoil
                                csoil=pcsoil[i][j2];
                            }
                            else { // for slope failure probability:

                                pf_ctr_phi = 1; // initializing counters for looping over samples of phi and csoil
                                pf_ctr_csoil = 1;
                                if ( PFDEPTH == 0 ) { 

                                    pfaill[j] = 0; // initializing layer-specific slope failure probability
                                    fos_meanl[j] = 0; 
                                    if ( NEWMARK == 1 ) { pnewmarkl[j] = 0; dnewmark_meanl[j] = 0; }
                                }
                            }
                            for ( pf_ctr = 1; pf_ctr <= pf_ctr_max; pf_ctr++ ) {
                                // starting loop over samples of phi and csoil:

                                if ( PF == 1 ) { // for slope failure probability:

                                    phi = ( (float)(pf_phi[pf_ctr_phi][i][j2]) ) / 10000; // defining phi

                                    if ( pdftype_csoil == 4 ) { // defining csoil:
                                        csoil = regpar1 * phi * 180 / PI + regpar2 + (float)pf_csoil[pf_ctr_csoil][i][j2];
                                        if ( csoil < 0 ) csoil = 0;
                                    }
                                    else csoil = (float)pf_csoil[pf_ctr_csoil][i][j2];

                                    pf_ctr_phi += 1; // updating counters of phi and csoil
                                    if ( pf_ctr_phi > pf_nsteps_phi ) {
                                        pf_ctr_phi = 1;
                                        pf_ctr_csoil += 1;
                                    }
                                }

                                area = pow ( csz , 2 ) * pow ( 1 - pow ( sin ( beta_xz ) , 2 ) * pow ( sin ( beta_yz ) , 2 ) , 0.5 )
                                    / ( cos ( beta_xz ) * cos ( beta_yz ) );

                                if ( PSEUDOSTATIC == 1 ) { // for pseudostatic seismic slope stability model:

                                    Segment_get(&seg_pga, &pga, x, y);
                                    seisforce = weight_sub / 9.81 * seiscoef * pga; //!!!CHECK dry or saturated density

                                    r_cum[j][pf_ctr] = csoil * area + ( weight_sub * cos ( slope ) - seisforce * sin( slope )) * tan ( phi ); // shear resistance
                                    t_cum[j][pf_ctr] = weight_sub * sin ( slope ) + seisforce * cos( beta_m ) + u_appr * sin ( gslope ); // shear force

                                } else {

                                    r_cum[j][pf_ctr] = csoil * area + weight_sub  * cos ( slope ) * tan ( phi );
                                    t_cum[j][pf_ctr] = weight_sub * sin ( slope ) + u_appr * sin ( gslope );
                                }

                                if ( r_cum[j][pf_ctr] > 0 && t_cum[j][pf_ctr] > 0 ) {

                                    fos_surf[j][pf_ctr] = r_cum[j][pf_ctr] / t_cum[j][pf_ctr];

                                    if ( NEWMARK == 1 ) {

                                        Segment_get(&seg_distance, &distance, x, y);
                                        if ( NEWMARKREF == 1 ) Segment_get(&seg_newmarkref, &newmarkref, x, y);
                                        else newmarkref = UNDEFPLUS;

                                        gnewmark = fnewmark( fos_surf[j][pf_ctr], slope, magnitude, distance, eqdepth, pars1, pars2, newmarkcoef, newmarkref, iariasmin, iariasmax, 
                                            pparias, narias );
                                        if ( PF == 0 && PFDEPTH == 0 ) dnewmark_surf[j][pf_ctr] = gnewmark[0];
                                        else {
                                            dnewmark_surf[j][pf_ctr] = gnewmark[2];
                                            pnewmark_surf[j][pf_ctr] = gnewmark[1];
                                        }
                                        free( gnewmark );
                                    }

                                    if ( PFDEPTH == 1 && PFSAMPLING <= 1 ) fos_mean += fos_surf[j][pf_ctr] / (float)pf_nsteps;
                                        // updating mean fos
                                    else if ( PFDEPTH == 1 && PFSAMPLING == 2 )
                                        fos_mean += fos_surf[j][pf_ctr] * pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                    else if ( PF == 1 && PFSAMPLING <= 1 ) fos_meanl[j] += fos_surf[j][pf_ctr] / (float)pf_nsteps;
                                    else if ( PF == 1 && PFSAMPLING == 2 ) fos_meanl[j] += fos_surf[j][pf_ctr] * pf_weight_geo[pf_ctr];

                                    if ( NEWMARK == 1 ) {
                                        if ( PFDEPTH == 1 && PFSAMPLING <= 1 ) dnewmark_mean += dnewmark_surf[j][pf_ctr] / (float)pf_nsteps;
                                            // updating mean fos
                                        else if ( PFDEPTH == 1 && PFSAMPLING == 2 )
                                            dnewmark_mean += dnewmark_surf[j][pf_ctr] * pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                        else if ( PF == 1 && PFSAMPLING <= 1 ) dnewmark_meanl[j] += dnewmark_surf[j][pf_ctr] / (float)pf_nsteps;
                                        else if ( PF == 1 && PFSAMPLING == 2 ) dnewmark_meanl[j] += dnewmark_surf[j][pf_ctr] * pf_weight_geo[pf_ctr];
                                    }
                                }

                                else { 
                                    fos_surf[j][pf_ctr] = UNDEFPLUS;
                                    if ( NEWMARK == 1 ) dnewmark_surf[j][pf_ctr] = 0;
                                }

                                if ( fos_surf[j][pf_ctr] < fos_min ) {
                                    fos_min = fos_surf[j][pf_ctr];
                                    depthfail_min = isd[j];

                                }
                                if ( NEWMARK == 1 ) {

                                    if ( dnewmark_surf[j][pf_ctr] > d_newmark ) {

                                        d_newmark = dnewmark_surf[j][pf_ctr];
                                        Segment_put(&seg_d_newmark, &d_newmark, x, y);
                                        Segment_flush(&seg_d_newmark);
                                    }
                                }

                                if ( fos_surf[j][pf_ctr] < 1 ) {
                                    fos_deep = fos_surf[j][pf_ctr];
                                    depthfail_deep = isd[j];
                                    Segment_put(&seg_depthfail_deep, &depthfail_deep, x, y);

                                    if ( PFDEPTH == 1 ) pfail += pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                    else if ( PF == 1 ) pfaill[j] += pf_weight_geo[pf_ctr]; // updating slope failure probability
                                }

                                if ( NEWMARK == 1 ) {

                                    if ( PFDEPTH == 1 ) pnewmark += pnewmark_surf[j][pf_ctr] * pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                    else if ( PF == 1 ) pnewmarkl[j] += pnewmark_surf[j][pf_ctr] * pf_weight_geo[pf_ctr];
                                }
                            }

                            if ( PF == 1 && PFDEPTH == 0 ) { // updating slope failure probability:

                                if ( pfaill[j] > pfail ) pfail = pfaill[j];
                                if ( fos_meanl[j] < fos_mean ) fos_mean = fos_meanl[j];
                                if ( NEWMARK == 1 ) {
                                    if ( pnewmarkl[j] > pnewmark ) pnewmark = pnewmarkl[j];
                                    if ( dnewmark_meanl[j] > dnewmark_mean ) dnewmark_mean = dnewmark_meanl[j];
                                }
                            }
                        }
                        j += 1;
		    }


               // Computing standard deviation of fos

                    if ( PF == 1 ) { // for slope failure probability mode:

                       fos_stdev = 0;
                       if ( NEWMARK == 1 ) dnewmark_stdev = 0;

                       for ( j = 1; j <= jmax; j++ ) { // loop over all layers:

                            if ( PFDEPTH == 0 ) fos_stdevl[j] = 0; // initializing standard deviation of fos
                            if ( PFDEPTH == 0 && NEWMARK == 1 ) dnewmark_stdevl[j] = 0;

                            for ( pf_ctr = 1; pf_ctr <= pf_ctr_max; pf_ctr++ ) { // loop over all samples

                                if ( PFDEPTH == 1 && PFSAMPLING <= 1 )
                                    fos_stdev += pow( fos_surf[j][pf_ctr] - fos_mean, 2 ) / (float)pf_nsteps;
                                    // updating standard deviation of fos
                                else if ( PFDEPTH == 1 && PFSAMPLING == 2 )
                                    fos_stdev += pow( fos_surf[j][pf_ctr] - fos_mean, 2 )
                                        * pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                else if ( PFSAMPLING <= 1 )
                                    fos_stdevl[j] += pow( fos_surf[j][pf_ctr] - fos_meanl[j], 2 ) / (float)pf_nsteps;
                                else if ( PFSAMPLING == 2 )
                                    fos_stdevl[j] += pow( fos_surf[j][pf_ctr] - fos_meanl[j], 2 ) * pf_weight_geo[pf_ctr];

                                if ( NEWMARK == 1 ) {

                                    if ( PFDEPTH == 1 && PFSAMPLING <= 1 )
                                        dnewmark_stdev += pow( dnewmark_surf[j][pf_ctr] - dnewmark_mean, 2 ) / (float)pf_nsteps;
                                    else if ( PFDEPTH == 1 && PFSAMPLING == 2 )
                                        dnewmark_stdev += pow( dnewmark_surf[j][pf_ctr] - dnewmark_mean, 2 )
                                            * pf_weight_geo[pf_ctr] * pf_weight_depth[j];
                                    else if ( PFSAMPLING <= 1 )
                                        dnewmark_stdevl[j] += pow( dnewmark_surf[j][pf_ctr] - dnewmark_meanl[j], 2 ) / (float)pf_nsteps;
                                    else if ( PFSAMPLING == 2 )
                                        dnewmark_stdevl[j] += pow( dnewmark_surf[j][pf_ctr] - dnewmark_meanl[j], 2 ) * pf_weight_geo[pf_ctr];
                                }
                            }

                            if ( fos_meanl[j] < fos_mean ) // if mean factor of safety is lower than old value:
                                fos_stdev = pow ( fos_stdevl[j]  * pf_ctr_max / ( pf_ctr_max - 1 ), 0.5 ); // updating standard deviation of fos

                            if ( NEWMARK == 1 ) {
                                if ( dnewmark_meanl[j] == dnewmark_mean ) dnewmark_stdev = pow ( dnewmark_stdevl[j]  * pf_ctr_max / ( pf_ctr_max - 1 ), 0.5 );
                            }
                        }

                        if ( NEWMARK == 1 ) {
                            Segment_put(&seg_pnewmark, &pnewmark, x, y);
                            Segment_put(&seg_d_newmark, &dnewmark_mean, x, y);
                            Segment_put(&seg_d_newmark_stdev, &dnewmark_stdev, x, y);
                        }

                        Segment_put(&seg_pfail, &pfail, x, y);
                        Segment_put(&seg_fos_min, &fos_mean, x, y);
                        Segment_put(&seg_fos_stdev, &fos_stdev, x, y);
                    }
                    else Segment_put(&seg_fos_min, &fos_min, x, y);

	        }
                if ( PF == 1) {
                    Segment_flush(&seg_pfail);
                    Segment_flush(&seg_fos_min);
                    Segment_flush(&seg_fos_stdev);

                    if ( NEWMARK == 1 ) {

                        Segment_flush(&seg_pnewmark);
                        Segment_flush(&seg_d_newmark);
                        Segment_flush(&seg_d_newmark_stdev);
                    }
                }
                else Segment_flush(&seg_fos_min);
                Segment_flush(&seg_depthfail_deep);
            }
        }


    // Printing results

        time_end = clock(); // end time of calculation
        time_elapsed = ( ( float ) ( time_end - time_start ) ) / CLOCKS_PER_SEC; // time needed for calculation

        if ( PF == 0 ) printf("Computing factor of safety for infinite slope ... completed in %.2f seconds. \n", time_elapsed);
        else printf("Computing slope failure probability for infinite slope ... completed in %.2f seconds. \n", time_elapsed);

    }


// Writing output raster maps

    foutrastd( mfos_min, seg_fos_min, 3 );
    foutrastd( mdepthfail_deep, seg_depthfail_deep, 3 );

    if ( PSEUDOSTATIC == 1 ) foutrastd( mfos_pseudo, seg_fos_pseudo, 3 );
    if ( NEWMARK == 1 ) foutrastd( mdnewmark, seg_d_newmark, 3 );

    if ( PF == 1 ) {

        if ( NEWMARK == 1 ) foutrastd( mpfail, seg_pnewmark, 3 );
        else foutrastd( mpfail, seg_pfail, 3 );

        if ( NEWMARK == 1 ) foutrastd( mfos_stdev, seg_d_newmark_stdev, 3 );
        else foutrastd( mfos_stdev, seg_fos_stdev, 3 );
    }
    if ( ADDITIONAL == 1 ) {

        foutrasti( mcumfos_crit, seg_cumfos_crit );
        foutrasti( mcumsurf, seg_cumsurf );
        foutrastd( mfos_deep, seg_fos_deep, 3 );
        foutrastd( mdepthfail_min, seg_depthfail_min, 3 );
        foutrastd( mdepthfail_crit, seg_depthfail_crit, 3 );
        foutrasti( mlyr_crit, seg_lyr_crit );
        foutrastd( mamin, seg_amin, 2 );
        foutrastd( mbmin, seg_bmin, 2 );
        foutrastd( mzbmin, seg_zbmin, 2 );
    }
    if ( ELLDENS == 0 ) foutrastd( mdepthfail_min, seg_depthfail_min, 3 );


// Closing output files

    fclose ( fsummary2 );
    if ( SUMMARY == 1 ) fclose ( fsummary1 );
    if ( ELLDENS == 0 ) fclose ( fprofile );
    else fclose ( fchange );


// Releasing segment data

    Segment_release(&seg_soilclass);
    Segment_release(&seg_elev);
    Segment_release(&seg_gwdepth);

    if ( PSEUDOSTATIC == 1 ) { Segment_release(&seg_pga); Segment_release(&seg_fos_pseudo); }
    if ( NEWMARK == 1 ) { Segment_release(&seg_distance); Segment_release(&seg_d_newmark); Segment_release(&seg_pnewmark); Segment_release(&seg_d_newmark_stdev); }
    if ( NEWMARKREF == 1 ) Segment_release(&seg_newmarkref);

    if ( CLASSINF == 1 ) for ( j = 0; j <= nly; j++ ) Segment_release(&seg_depth[j]);

    if ( LAYERS == 1 && SEEPAGE == 2 ) {
        Segment_release(&seg_dzdx);
        Segment_release(&seg_dzdy);
    }

    if ( ELLIPS == 1 && USESEG == 1 ) { 

        for ( j = 0; j <= nly+1; j++ ) Segment_release(&seg_epd[j]);
        Segment_release(&seg_epe);
        if ( PSEUDOSTATIC == 1 ) Segment_release(&seg_epga);
        Segment_release(&seg_ezslip);
        Segment_release(&seg_esoilclass);
    }

    Segment_release(&seg_fos_min);
    Segment_release(&seg_pfail);
    Segment_release(&seg_fos_stdev);
    Segment_release(&seg_depthfail_deep);
    Segment_release(&seg_fos_deep);

    if ( ADDITIONAL == 1 ) {
        Segment_release(&seg_cumfos_crit);
        Segment_release(&seg_cumsurf);
        Segment_release(&seg_depthfail_min);
        Segment_release(&seg_depthfail_crit);
        Segment_release(&seg_lyr_crit);
        Segment_release(&seg_amin);
        Segment_release(&seg_bmin);
        Segment_release(&seg_zbmin);
    }
    if ( ELLDENS == 0 ) Segment_release(&seg_depthfail_min);


// Freeing memory

    free(TEMPDIR); free(MAPSET); free(ldata[0]); free(ldata); free(r_cum[0]); free(t_cum[0]); free(fos_surf[0]); free(r_cum); free(t_cum); free(fos_surf);

    free(pgamma0[0]); free(pcsoil0[0]); free(pphi0[0]); free(ptheta0[0]); free(pgamma[0]); free(pcsoil[0]); free(pphi[0]); free(ptheta[0]); free(param_depth0[0]);  
    free(pgamma0); free(pcsoil0); free(pphi0); free(ptheta0); free(pgamma); free(pcsoil); free(pphi); free(ptheta); free(param_depth0);
        
    if ( ELLIPS == 1 && USESEG == 0 ) {

        free(epe[0]), free(epe); if ( PSEUDOSTATIC == 1 ) { free(epga[0]), free(epga); }
        free(ezslip[0]); free(ezslip); free(eparam_soilclass[0]); free(eparam_soilclass);
        for ( j=0; j <= jmax+1; j++) { free(epd[j][0]); free(epd[j]); } free(efnh); free(epsc_nh);
    }
    
    if ( PF == 1) {

        free(pcsoils0[0]); free(pcsoilmin0[0]); free(pcsoilmax0[0]); free(pphis0[0]); free(pphimin0[0]); free(pphimax0[0]); free(pcsoils[0]); free(pcsoilmin[0]); free(pcsoilmax[0]); 
        free(pphis[0]); free(pphimin[0]); free(pphimax[0]);
        free(pcsoils0); free(pcsoilmin0); free(pcsoilmax0); free(pphis0); free(pphimin0); free(pphimax0); free(pcsoils); free(pcsoilmin); free(pcsoilmax); free(pphis); free(pphimin); free(pphimax);

        free(pf_pdf_peak_phi0[0]); free(pf_pdf_peak_phi0); free(pf_pdf_peak_phi[0]); free(pf_pdf_peak_phi);
        free(pf_pdf_peak_csoil0[0]); free(pf_pdf_peak_csoil0); free(pf_pdf_peak_csoil[0]); free(pf_pdf_peak_csoil);
        for ( pf_ctr = 0; pf_ctr <= pf_nsteps_phi; pf_ctr++ ) { free(pf_phi0[pf_ctr][0]); free(pf_phi0[pf_ctr]); free(pf_phi[pf_ctr][0]); free(pf_phi[pf_ctr]); }
        for ( pf_ctr = 0; pf_ctr <= pf_nsteps_csoil; pf_ctr++ ) { free(pf_csoil0[pf_ctr][0]); free(pf_csoil0[pf_ctr]); free(pf_csoil[pf_ctr][0]); free(pf_csoil[pf_ctr]); }
    }

    if ( NEWMARK == 1 ) { free( dnewmark_surf[0] ); free( dnewmark_surf ); }}

    return 0;
}
