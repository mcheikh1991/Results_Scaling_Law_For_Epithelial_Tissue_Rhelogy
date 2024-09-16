import numpy as np
import matplotlib as mpl
import os

Location = os.getenv('HOME')+ "/Dropbox/python/Tools/"
exec(compile(open(Location+"Functions.py", "rb").read(), Location+"Functions.py", 'exec'))


def writeSbatchInfor(f, N, n):
  f.write('#!/bin/bash\n')
  f.write('#SBATCH -J LE%d\n' % N)
  f.write('#SBATCH --partition=32GB\n')
  f.write('#SBATCH --output=log%d.out\n' %N)
  f.write('#SBATCH --error=log%d.out\n' %N)
  f.write('#SBATCH --nodes=1          # number of nodes requested by user\n')
  f.write('#SBATCH --ntasks=%d        # number of total tasks\n' % n)
  f.write('#SBATCH --time=1-23:00:00   # run time, format: D-H:M:S (max wallclock time)\n')
  f.write('#SBATCH -A Dobrovinski\n')
  f.write('##SBATCH --mail-user=mohamadibrahim.cheikh@utsouthwestern.edu\n')
  f.write('##SBATCH --mail-type=all\n')
  f.write('\n')
  f.write('module purge\n')
  f.write('module add shared slurm gcc/8.3.0 \n')


# system properties:
#====================================================
platform = 'linux-biohpc' #['windows','linux','linux-biohpc']
omp_threads = 10
dim         = 3
csv_output_file   = "output_0_Loading_Tri/track.csv"
vtk_output_folder = "output_0_Loading_Tri"
mesh_input_folder = "solid/tri_mesh"

# time properties:
#====================================================
time_dt               = 0.025  # [ms]
time_end              = 60000. # [ms]
use_real_time_to_save = 'true'

# if use_real_time_to_save == true
time_dtime_save_vtk   =   1000. # [ms]
time_dtime_save_csv   =     20. # [ms]
time_dtime_print      =    100. # [ms]

# if use_real_time_to_save == false
time_ntimesteps_save_vtk = 100
time_ntimesteps_save_csv = 100
time_ntimesteps_print =    100

# boundary properties:
#====================================================
boundary_n       = 1 # number of boundaries

# Boundary Type;
# no_solid_bc=0, force_add=1, force_insert=2, dist_fixed=3, dist_speed=4, force_speed=5, 
# symmetry=6 , force_speed_ellipse=7, dist_speed_ellipse=8, dist_speed_x_linear_y=9, 
# axis_fixed=10, force_add_ellipse=11

boundary_name, boundary_type, = [], []
boundary_dir, boundary_fm = [], []
boundary_a, boundary_b = [], []
boundary_dtheta0, boundary_dthetadt = [], []
boundary_k = []

boundary_name.append("forced_point")
boundary_type.append(7)
boundary_a.append(251.5)
boundary_b.append( 88.0)
boundary_dir.append(1)
boundary_dtheta0.append(0.0)
boundary_dthetadt.append(1.84E-06)
boundary_k.append(1.0)

# box properties:
#====================================================
boxes_max_points      = 100
boxes_max_lines       = 100
boxes_max_triangles   = 100
boxes_max_squares     = 100
boxes_m               = 1.0

# le force properties:
#====================================================
# [force_le_null=0, force_le_fixed_disp=1, force_le_initial_disp=2, force_le_initial_disp_multi=3, force_le_initial_disp_cond=4]
le_type = 2 
le_k    = 1.0

# barrier force properties:
#====================================================
bf_type                = 1     # [ force_bf_null=0, force_bf_all_points=1, force_bf_boundary_points=2]

bf_bool_force_line     = 'true'
bf_type_force_line     = 1     # 0=const, 1=linear, 2=exp
bf_d_min_line          = 0.75  # [um]
bf_c1_line             = 0.0
bf_c2_line             = 5.0

bf_bool_force_point    = 'true'
bf_type_force_point    = 1     # 0=const, 1=linear, 2=exp
bf_d_min_point         = 0.5  # [um]
bf_c1_point            = 0.0
bf_c2_point            = 5.0

bf_bool_force_triangle = 'false'
bf_type_force_triangle = 1     # 0=const, 1=linear, 2=exp
bf_d_min_triangle      = 0.75  # [um]
bf_c1_triangle         = 0.0
bf_c2_triangle         = 5.0

bf_bool_force_square   = 'false'
bf_type_force_square   = 1     # 0=const, 1=linear, 2=exp
bf_d_min_square        = 0.75  # [um]
bf_c1_square           = 0.0
bf_c2_square           = 5.0

# l0-dynamicsproperties:
#====================================================
ld_type               = 0      # 0=null, 1=alaways, 2=conditional
ld_lambda             = 0.0
ld_tau                = 2000
 
# soft force properties:
#====================================================
sf_type               = 1      # 0=null, 1=ellipsoid, 2=twoellipsoid
sf_a                  = 251.5
sf_b                  =  88.0
sf_c                  =  88.0
sf_bd                 =  0.0
sf_bv                 =  0.0
sf_ds                 =  1.0
sf_sf                 = 100.0
print_sf_points       = 'false'

# surface tension properties:
#====================================================
st_type               = 0      # 0=null, 1=ellipsoid, 2=twoellipsoid
st_gamma              = 0.0

#creating the input file
if platform == 'windows':
  filename = 'windows_batch.bat'
elif platform == 'linux':
  filename = 'linux_batch.sh'
elif platform == 'linux-biohpc':
  filename = 'sbatch_petsc.sh'

f = open(filename, 'w')
if platform == 'windows':
  f.write('@echo off\n')
  f.write('Solid_Solver.exe ')
elif platform == 'linux':
  f.write('#!/bin/bash\n')
  f.write('./Solid_Solver.o ')
elif platform == 'linux-biohpc':
  writeSbatchInfor(f, 0, omp_threads)
  f.write('./Solid_Solver.o ')

f.write("-print_sf_points %s " % print_sf_points)

f.write('-omp %d ' % omp_threads)
f.write('-dim %d ' % dim)
f.write('-csv_output_file %s ' % csv_output_file)
f.write('-vtk_output_folder %s ' % vtk_output_folder)
f.write('-mesh_input_folder %s ' % mesh_input_folder)

f.write("-time_dt %lf " % time_dt)
f.write("-time_end %lf " % time_end)   
f.write("-use_real_time_to_save %s " % use_real_time_to_save)
f.write("-time_dtime_save_vtk %lf " % time_dtime_save_vtk)
f.write("-time_dtime_save_csv %lf " % time_dtime_save_csv)
f.write("-time_dtime_print %lf " % time_dtime_print)

f.write('-n_bound %d ' % boundary_n)

for i in range(boundary_n):
  if (boundary_type[i] == 11):
    f.write('-bound_%d_name %s ' % (i, boundary_name[i]))
    f.write('-bound_%d_type %d ' % (i, boundary_type[i]))
    f.write('-bound_%d_a %lf ' % (i, boundary_a[i]))
    f.write('-bound_%d_b %lf ' % (i, boundary_b[i]))
    f.write('-bound_%d_dir %d ' % (i, boundary_dir[i]))
    f.write('-bound_%d_fm %lf ' % (i, boundary_fm[i]))
  elif (boundary_type[i] == 7):
    f.write('-bound_%d_name %s ' % (i, boundary_name[i]))
    f.write('-bound_%d_type %d ' % (i, boundary_type[i]))
    f.write('-bound_%d_a %lf ' % (i, boundary_a[i]))
    f.write('-bound_%d_b %lf ' % (i, boundary_b[i]))
    f.write('-bound_%d_dir %d ' % (i, boundary_dir[i]))
    f.write('-bound_%d_dtheta0 %lf ' % (i, boundary_dtheta0[i]))
    f.write('-bound_%d_dthetadt %lf ' % (i, boundary_dthetadt[i]))
    f.write('-bound_%d_k %lf ' % (i, boundary_k[i]))
  else:
    raise NameError("Boundary not defined yet")


f.write("-boxes_max_points %d " % boxes_max_points)
f.write("-boxes_max_lines %d " % boxes_max_lines)
f.write("-boxes_max_triangles %d " % boxes_max_triangles)
f.write("-boxes_max_squares %d " % boxes_max_squares)
f.write("-boxes_m %lf " % boxes_m)

f.write("-le_type %d " % le_type)
if le_type == 1 or le_type == 2:
  f.write("-le_k %lf " % le_k)

f.write("-bf_type %d " % bf_type)

f.write('-bf_bool_force_line %s ' % bf_bool_force_line)
if bf_bool_force_line == 'true':
  f.write("-bf_type_force_line %d " % bf_type_force_line) 
  f.write("-bf_d_min_line %lf " % bf_d_min_line)
  f.write("-bf_c1_line %lf " % bf_c1_line)
  f.write("-bf_c2_line %lf " % bf_c2_line)

f.write('-bf_bool_force_point %s ' % bf_bool_force_point)
if bf_bool_force_point == 'true':
  f.write("-bf_type_force_point %lf " % bf_type_force_point)
  f.write("-bf_d_min_point %lf " % bf_d_min_point)
  f.write("-bf_c1_point %lf " % bf_c1_point)
  f.write("-bf_c2_point %lf " % bf_c2_point)

f.write('-bf_bool_force_triangle %s ' % bf_bool_force_triangle)
if bf_bool_force_triangle == 'true':
  f.write("-bf_type_force_triangle %d " % bf_type_force_triangle) 
  f.write("-bf_d_min_triangle %lf " % bf_d_min_triangle)
  f.write("-bf_c1_triangle %lf " % bf_c1_triangle)
  f.write("-bf_c2_triangle %lf " % bf_c2_triangle)

f.write('-bf_bool_force_square %s ' % bf_bool_force_square)
if bf_bool_force_square == 'true':
  f.write("-bf_type_force_square %lf " % bf_type_force_square)
  f.write("-bf_d_min_square %lf " % bf_d_min_square)
  f.write("-bf_c1_square %lf " % bf_c1_square)
  f.write("-bf_c2_square %lf " % bf_c2_square)

f.write("-ld_type %d " % ld_type)
if ld_type != 0.0:
  f.write("-ld_lambda %lf " % ld_lambda)
  f.write("-ld_tau %lf " % ld_tau)

f.write("-sf_type %d " % sf_type)
if sf_type != 0.0:
  f.write("-sf_a %lf " % sf_a)
  f.write("-sf_b %lf " % sf_b)
  f.write("-sf_c %lf " % sf_c)
  f.write("-sf_bd %lf " % sf_bd)
  f.write("-sf_bv %lf " % sf_bv)
  f.write("-sf_ds %lf " % sf_ds)
  f.write("-sf_sf %lf " % sf_sf)

f.write("-st_type %d " % st_type)
if st_type != 0.0:
  f.write("-st_gamma %lf " % st_gamma)

f.write("\n")
f.close()
