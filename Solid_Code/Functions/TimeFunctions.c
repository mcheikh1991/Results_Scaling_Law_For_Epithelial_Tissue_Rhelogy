double myclock() 
{
  static time_t t_start = 0;  // Save and subtract off each time

  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  if( t_start == 0 ) t_start = ts.tv_sec;

  return (double) (ts.tv_sec - t_start) + ts.tv_nsec * 1.0e-9;
}

void DefaultTimeParameters()
{  
  // Default Time
  times.real_time =          0.0; // [ms]
  times.dt =                 0.1; // [ms]
  times.end_time =      100000.0; // [ms]

  times.use_real_time_to_save = true;
  times.dtime_save_vtk =  1000.0; // [ms]
  times.dtime_save_csv =   100.0; // [ms]
  times.dtime_print =      100.0; // [ms]

  times.ntimeSteps_save_vtk = 1000;
  times.ntimeSteps_save_csv =   10;
  times.ntimeSteps_print =      10;

  times.reload      =      false;
  times.reload_time =        0.0; // [ms]
  times.reload_step =          0; // [ms]
}

bool ReadTimeParameters(const char* str_argv1, const char* str_argv2)
{
  if (strcmp(str_argv1,"-time_dt") == 0){     
    times.dt = atof(str_argv2); 
    return true;}

  if ((strcmp(str_argv1,"-time_end") == 0)|| (strcmp(str_argv1,"-time_end_time") == 0)){    
    times.end_time = atof(str_argv2);
    return true;}

  if (strcmp(str_argv1,"-use_real_time_to_save") == 0){ 
    if (strcmp(str_argv2,"false") == 0) {
      times.use_real_time_to_save = false; }
    else if (strcmp(str_argv2,"true") == 0) {
      times.use_real_time_to_save = true; }
    else{
      printf("    ERROR!! Unknown bool command '%s'; should be either 'false' or 'true'\n", str_argv1);
      exit(EXIT_FAILURE);}
    return true;}

  if (strcmp(str_argv1,"-time_dtime_save_vtk") == 0){    
    times.dtime_save_vtk = atof(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-time_dtime_save_csv") == 0){    
    times.dtime_save_csv = atof(str_argv2);
    return true;}

  if (strcmp(str_argv1,"-time_dtime_print") == 0){    
    times.dtime_print = atof(str_argv2);
    return true;}

  if (strcmp(str_argv1,"-time_ntimesteps_save_vtk") == 0){    // Boundary or moving point index
    times.ntimeSteps_save_vtk = atol(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-time_ntimesteps_save_csv") == 0){    // Boundary or moving point index
    times.ntimeSteps_save_csv = atol(str_argv2);
    return true;}

  if (strcmp(str_argv1,"-time_ntimesteps_print") == 0){       // Boundary or moving point index
    times.ntimeSteps_print = atol(str_argv2); 
    return true;}
  return false;
}

bool ReadReloadParameters(const char* str_argv1, const char* str_argv2)
{
  if (strcmp(str_argv1,"-reload") == 0){ 
    if (strcmp(str_argv2,"false") == 0) {       times.reload = false; }
    else if (strcmp(str_argv2,"true") == 0) {   times.reload = true; }
    else{
      printf("    ERROR!! Unknown bool command '%s'; should be either 'false' or 'true'\n", str_argv1);
      exit(EXIT_FAILURE);}
    return true;}

  if (strcmp(str_argv1,"-reload_time") == 0){    
    times.reload_time = atof(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-reload_step") == 0){     
    times.reload_step = atoi(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-reload_folder") == 0){     
    strcpy(reloadFolder, str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-reload_vtk_file")==0){      // vtk_reload_file
    strcpy(reloadVTKFile, str_argv2); 
    return true;}

  return false;
}

void ReportStatus()
{
  times.ttotal = myclock() - times.tstart;
  if (VC.Type == FORCE_VC_CONST_CYLIN || VC.Type == FORCE_VC_CONST_CELLS)
    printf("Iteration: %d - Time: %lf ms - Computational time: %lf sec - Volume: %lf\n", times.time_step, times.real_time, times.ttotal, VC.volume);
  else
    printf("Iteration: %d - Time: %lf ms - Computational time: %lf sec\n", times.time_step, times.real_time, times.ttotal);
  times.tstart = myclock(); 
  times.time_print += times.dtime_print;
  fflush(stdout);
}
