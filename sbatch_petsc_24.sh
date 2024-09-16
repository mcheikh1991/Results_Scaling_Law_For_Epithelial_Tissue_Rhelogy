#!/bin/bash
#SBATCH -J LE_LD24
#SBATCH --partition=super
#SBATCH --output=log24.out
#SBATCH --error=log24.out
#SBATCH --nodes=1          # number of nodes requested by user
#SBATCH --ntasks=5         # number of total tasks
#SBATCH --time=7-23:00:00  # run time, format: D-H:M:S (max wallclock time)
#SBATCH -A Dobrovinski
##SBATCH --mail-user=mohamadibrahim.cheikh@utsouthwestern.edu
##SBATCH --mail-type=all

module purge
module add shared slurm intel/2021.3.0
./Solid_Solver.i -print_sf_points false -omp 5 -dim 3 -csv_output_file output_24/track.csv -vtk_output_folder output_24 -mesh_input_folder solid/hex_mesh_a125.75-b44-c44 -time_dt 0.0001 -time_end 100001.0 -use_real_time_to_save true -time_dtime_save_vtk 1000.0 -time_dtime_save_csv 1.000000 -time_dtime_print 100.000000 -n_bound 1 -bound_0_name moving_point -bound_0_type 11 -bound_0_a 125.75 -bound_0_b 44.0 -bound_0_dir -1 -bound_0_fm 0.01 -le_type 2 -le_k 0.1 -bf_type 2 -bf_bool_force_line true -bf_type_force_line 1 -bf_d_min_line 0.75 -bf_c1_line 0.0 -bf_c2_line 5.0 -bf_bool_force_point true -bf_type_force_point 1 -bf_d_min_point 0.75 -bf_c1_point 0.0 -bf_c2_point 5.0 -bf_bool_force_triangle false -bf_bool_force_square false -bf_file_set_points solid/hex_mesh_a125.75-b44-c44/bound_moving_point.dat -ld_type 1 -ld_lambda 0.0 -ld_tau 2000 -sf_type 1 -sf_a 125.75 -sf_b 44.0 -sf_c 44.0 -sf_ds 1.0 -sf_sf 5500.0 -st_type 0 -reload false -reload_time 54000 -reload_vtk_file output_24/S_54.vtk -drag_c 0.001

