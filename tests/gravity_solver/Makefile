SRC = ../../source
APPNAME = app

OPTS = sdims_make=3
OPTS += nx1_make=128
OPTS += nx2_make=128
OPTS += nx3_make=128
OPTS += ddx1_make=1 
OPTS += ddx2_make=1 
OPTS += ddx3_make=1 
OPTS += ngc_make=2
OPTS += nas_make=0

OPTS += GEOMETRY_CARTESIAN_UNIFORM

OPTS += X1U_OUTFLOW
OPTS += X2L_OUTFLOW
OPTS += X2U_OUTFLOW
OPTS += X3L_OUTFLOW
OPTS += X3U_OUTFLOW

OPTS += SKIP_HYDRO

OPTS += cfl_make=0.4_rp
OPTS += tmax_make=3.6e4_rp
OPTS += stepmax_make=1
OPTS += wctmax_make=85000.0_rp

OPTS += LIM2ND_REC
OPTS += HLLC_FLUX
OPTS += Mach_cutoff_make=1.0e-4_rp
OPTS += RK3_STEPPER

OPTS += USE_GRAVITY
OPTS += USE_GRAVITY_SOLVER
OPTS += gs_tol_make=1.000e-10_rp 
OPTS += CONST_GRAV_UNITY

OPTS += GMG_PRECONDITIONER
OPTS += gmg_niter_coarse_make=5
OPTS += gmg_niter_presmooth_make=3
OPTS += gmg_niter_postsmooth_make=5
OPTS += gmg_max_level_make=7
OPTS += gmg_omega_make=0.85_rp

OPTS += OUTPUT_DT
OPTS += dt_dump_make=100.0_rp
OPTS += nsteps_dump_make=200
OPTS += info_terminal_rate_make=100
OPTS += dt_restart_make=1.0e4_rp

OPTS += USE_DOUBLE_PRECISION
OPTS += LITTLE_ENDIAN

#--------------------------------#

OPTS +=nvars_make=sdims_make+2+nas_make

# Include machine-specific configuration
include ../../Makefile.config



