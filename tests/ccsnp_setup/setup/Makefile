SRC = ../../../source
DATA = ../../../data
APPNAME = app

#=================================#
#GRID OPTIONS
#=================================#

OPTS = sdims_make=3
OPTS += nx1_make=128
OPTS += nx2_make=128
OPTS += nx3_make=128
OPTS += ddx1_make=2
OPTS += ddx2_make=2
OPTS += ddx3_make=1
OPTS += ngc_make=2
OPTS += nas_make=12
OPTS += nspecies_make=12
OPTS += nreacs_make=38

#=================================#
#GRID GEOMETRY
#=================================#

OPTS += GEOMETRY_CARTESIAN_UNIFORM

#=================================#
#BOUNDARY CONDITIONS
#=================================#

OPTS += X1_PERIODIC #overwritten by USE_INTERNAL_BOUNDARIES
OPTS += X2_PERIODIC #overwritten by USE_INTERNAL_BOUNDARIES
OPTS += X3_PERIODIC #overwritten by USE_INTERNAL_BOUNDARIES
OPTS += USE_INTERNAL_BOUNDARIES

#=================================#
#TIME STEPPING
#=================================#

OPTS += cfl_make=0.8_rp
OPTS += RK3_STEPPER

#=================================#
#GODUNOV
#=================================#

OPTS += PPH_REC
OPTS += LHLLD_FLUX
OPTS += Mach_cutoff_make=1.0e-4_rp

#=================================#
#MHD
#=================================#

OPTS += USE_MHD
OPTS += Mach_ct_make=1.0e-4_rp

#=================================#
#GRAVITY
#=================================#

OPTS += USE_GRAVITY
OPTS += EVOLVE_ETOT
OPTS += USE_WB

#=================================#
#EOS
#=================================#

OPTS += HELMHOLTZ_EOS
OPTS += USE_COULOMB_CORRECTIONS
OPTS += path_to_helm_table=\"$(DATA)/helm_table_timmes_x1.dat\"
OPTS += helm_nT_make=101
OPTS += helm_nrho_make=271
OPTS += helm_ltlo_make=3.0_rp
OPTS += helm_lthi_make=13.0_rp
OPTS += helm_ldlo_make=-12.0_rp
OPTS += helm_ldhi_make=15.0_rp
OPTS += USE_FASTEOS

#=================================#
#NUCLEAR NETWORK
#=================================#

OPTS += ADVECT_SPECIES
OPTS += USE_NUCLEAR_NETWORK
OPTS += USE_ELECTRON_SCREENING
OPTS += NUCLEAR_NETWORK_GODUNOV
OPTS += NUCLEAR_NETWORK_BE
OPTS += nn_Tmin_make=1e6_rp
OPTS += nn_nltol_make=1.0e-12_rp
OPTS += nn_ltol_make=1.0e-12_rp

#=================================#
#THERMAL NEUTRINOS
#=================================#

OPTS += USE_NEULOSS

#=================================#
#STOPPING CRITERIA
#=================================#

OPTS += tmax_make=4000.0_rp
OPTS += wctmax_make=85000.0_rp
OPTS += stepmax_make=100000000

#=================================#
#OUTPUT OPTIONS
#=================================#

OPTS += OUTPUT_DT
OPTS += dt_dump_make=25.0_rp
OPTS += info_terminal_rate_make=1
OPTS += dt_restart_make=40000.0_rp

#planes
OPTS += SAVE_PLANES
OPTS += nplanes_x1_make=0
OPTS += nplanes_x2_make=0
OPTS += nplanes_x3_make=1
OPTS += planes_dt_dump_make=0.2_rp

#compression
OPTS += compression_factor_make=2

#ra-iles
OPTS += SAVE_RPROFS
OPTS += rprofs_dt_dump_make=0.1_rp

#=================================#
#PRECISION AND ENDIANITY
#=================================#

OPTS += USE_DOUBLE_PRECISION
OPTS += LITTLE_ENDIAN

#=================================#
#RESTART
#=================================#

OPTS += RESTART_LAST

#--------------------------------#

OPTS +=nvars_make=sdims_make+2+nas_make

# Include machine-specific configuration
include ${SRC}/../Makefile.config 
