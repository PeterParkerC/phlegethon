import matplotlib.pyplot as plt
import numpy as np
import math
import os
import time
from scipy.integrate import cumulative_trapezoid as cumtrapz
import subprocess as sp
import sys
import fileinput
import scipy.interpolate as interp1d
from scipy.interpolate import CubicSpline as interp1d
from scipy.optimize import curve_fit

part = 'genoa.p'
X_nint= 5
NT = 200
Nrho = 200
ddx1 = 5
ddx2 = 5
log10_T_min= 1.0
log10_T_max= 8.0
log10_rho_min= -20.0
log10_rho_max= 0.0

#-------------------------------------------------------------------------------#

run = True
compile_code = True

for idx in range(X_nint+1):

        path = './id_%d'%(idx)

        os.makedirs(path,exist_ok=True)

        cmd = 'cp -r setup/* %s/'%(path)
        process = sp.Popen(cmd,shell=True)
        process.wait()
       
        phl_conf = fileinput.input('%s/Makefile'%(path),inplace=1)
        for line in phl_conf:

           if("OPTS = Nrho" in line):
             line = "OPTS = Nrho=%d \n"%(Nrho)

           if("OPTS += NT" in line):
             line = "OPTS += NT=%d \n"%(NT)
 
           if("OPTS += log10_T_min" in line):
             line = "OPTS += log10_T_min=%.3f_rp \n"%(log10_T_min)
  
           if("OPTS += log10_T_max" in line):
             line = "OPTS += log10_T_max=%.3f_rp \n"%(log10_T_max)
  
           if("OPTS += log10_rho_min" in line):
             line = "OPTS += log10_rho_min=%.3f_rp \n"%(log10_rho_min)
 
           if("OPTS += log10_rho_max" in line):
             line = "OPTS += log10_rho_max=%.3f_rp \n"%(log10_rho_max)

           if("OPTS += ddx1" in line):
             line = "OPTS += ddx1=%d \n"%(ddx1)

           if("OPTS += ddx2" in line):
             line = "OPTS += ddx2=%d \n"%(ddx2)

           if("OPTS += X_nint" in line):
             line = "OPTS += X_nint=%d \n"%(X_nint)

           if("OPTS += sim_index" in line):
             line = "OPTS += sim_index=%d \n"%(idx)

           sys.stdout.write(line)

        phl_conf.close()

        ncores = ddx1*ddx2

        job_spt = fileinput.input('%s/pig_table.job'%(path),inplace=1)
       
        for line in job_spt:
         if("#SBATCH -p" in line):
          line = "#SBATCH -p %s \n"%(part)
         if("#SBATCH -N 1" in line):
          line = "#SBATCH -N %d \n"%(1)
         if("#SBATCH -n 1" in line):
          line = "#SBATCH -n %d \n"%(ncores)
         if("mpirun -n 1 ./run.generate_pig_table" in line):
          line = "mpirun -n %d ./run.generate_pig_table \n"%(ncores)
         sys.stdout.write(line)
        job_spt.close()

        if(compile_code):
            os.chdir(path)
            cmd = 'make clean; make'
            process = sp.Popen(cmd,shell=True)
            process.wait()
            os.chdir('../')

        if(run):
            os.chdir(path)
            cmd = 'sbatch pig_table.job'
            process = sp.Popen(cmd,shell=True)
            process.wait()
            os.chdir('../')

