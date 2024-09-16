/*
// Everything to lower case
for (int i=1; i<argc; i++)
  for(int j = 0; argv[i][j]; j++)
    argv[i][j] = tolower(argv[i][j]);
*/
for (int i=1; i<argc; i++)
{
  //printf("Reading command: %s\n", argv[i]);
  
  // Optional commands:
  //====================================================
  if (strcmp(argv[i],"-print_sf_points") == 0){ 
    if (strcmp(argv[i+1],"false") == 0) {
      print_sf_points = false; i++;}
    else if (strcmp(argv[i+1],"true") == 0) {
      print_sf_points = true; i++;}
    else{
      printf("    ERROR!! Unknown bool command '%s'; should be either 'false' or 'true'\n", argv[i]);
      exit(EXIT_FAILURE);}
    continue;} 

  if (strcmp(argv[i],"-print_mesh_points") == 0){ 
    if (strcmp(argv[i+1],"false") == 0) {
      print_mesh_points = false; i++;}
    else if (strcmp(argv[i+1],"true") == 0) {
      print_mesh_points = true; i++;}
    else{
      printf("    ERROR!! Unknown bool command '%s'; should be either 'false' or 'true'\n", argv[i]);
      exit(EXIT_FAILURE);}
    continue;} 

  if (strcmp(argv[i],"-print_mesh_lines") == 0){ 
    if (strcmp(argv[i+1],"false") == 0) {
      print_mesh_lines = false; i++;}
    else if (strcmp(argv[i+1],"true") == 0) {
      print_mesh_lines = true; i++;}
    else{
      printf("    ERROR!! Unknown bool command '%s'; should be either 'false' or 'true'\n", argv[i]);
      exit(EXIT_FAILURE);}
    continue;} 

  // System Properties:
  //====================================================
  if (strcmp(argv[i],"-omp") == 0) {                // OMP Threads
    omp_nthreads = atol(argv[i+1]); i++;     continue;}

  if (strcmp(argv[i],"-dim") == 0) {                // DIM
    DIM = atol(argv[i+1]); i++;              continue;}

  if (strcmp(argv[i],"-csv_output_file") == 0) {    // csv_output_file
    strcpy(csvOutputFile, argv[i+1]); i++;   continue;}

  if (strcmp(argv[i],"-vtk_output_folder") == 0) {  // vtk_output_folder
    strcpy(vtkOutputFolder, argv[i+1]); i++; continue;}

  if (strcmp(argv[i],"-mesh_input_folder") == 0) {  // mesh_input_folder
    strcpy(meshInputFolder, argv[i+1]); i++; continue;}

  if (strcmp(argv[i],"-vtk_input_file") == 0) {     // vtk_input_file
    strcpy(vtkInputFile, argv[i+1]); i++; continue;}
    
  // Time & Reload Properties:
  //====================================================
  if (ReadTimeParameters(argv[i],argv[i+1])) {i++; continue;}
  if (ReadReloadParameters(argv[i],argv[i+1])) {i++; continue;}

  // Boundary Properties:
  //====================================================
  if (strcmp(argv[i],"-n_bound") == 0)    // Number of boundaries
  {
    boundaries.N = atol(argv[i+1]); i++; 
    boundaries.groups = (Boundary_Points_Struct *)malloc(boundaries.N*sizeof(Boundary_Points_Struct));
    
    if (boundaries.N>0) i++;
    else if (boundaries.N<0){ printf("    ERROR!! Number of boundaries can not be %d\n", boundaries.N);   exit(EXIT_FAILURE);}

    char boundaryName[256];
    // Now read type of each boundary
    for (int j=0; j<boundaries.N; j++)
    {  
      // Boundary Name
      sprintf(boundaryName, "-bound_%d_name", j);
      if (strcmp(argv[i],boundaryName) == 0){  
        boundaries.groups[j].Name = malloc(256);
        sprintf(boundaries.groups[j].Name, "%s", argv[i+1]);  i+=2;}
      else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
      memset(boundaryName,0,256);

      // Boundary points index
      sprintf(boundaryName, "./%s/bound_%s.dat", meshInputFolder, boundaries.groups[j].Name);
      boundaries.groups[j].N = GetFileSize(boundaryName);
      boundaries.groups[j].I = (int *)malloc(boundaries.groups[j].N*sizeof(int));
      ReadIntFile(boundaryName, boundaries.groups[j].N, boundaries.groups[j].I);
      memset(boundaryName,0,256);
      if (boundaries.groups[j].N < 1){  
        printf("    ERROR!! Boundary '%s' has '%d' points. It should be 1 or more\n", boundaries.groups[j].Name, boundaries.groups[j].N);   exit(EXIT_FAILURE);}

      // Boundary Type
      sprintf(boundaryName, "-bound_%d_type", j);
      if (strcmp(argv[i],boundaryName) == 0){          
        boundaries.groups[j].Type = (SOLID_BC_TYPE) atol(argv[i+1]); i+=2;}
      else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
      memset(boundaryName,0,256);

      if (boundaries.groups[j].Type == FORCE_ADD)               // Add forces linearly
      {
        sprintf(boundaryName, "-bound_%d_fx", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].Fx = atof(argv[i+1]); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);

        sprintf(boundaryName, "-bound_%d_fy", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].Fy = atof(argv[i+1]); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);

        if (DIM==3)
        {
          sprintf(boundaryName, "-bound_%d_fz", j);
          if (strcmp(argv[i],boundaryName) == 0){          
            boundaries.groups[j].Fz = atof(argv[i+1]); i+=2; }
          else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
          memset(boundaryName,0,256);
        }
        else
          boundaries.groups[j].Fz = 0;

      }
      else if (boundaries.groups[j].Type == FORCE_ADD_ELLIPSE)  // Add forces on an ellipsoid
      {
        sprintf(boundaryName, "-bound_%d_a", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].A = atof(argv[i+1]); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);

        sprintf(boundaryName, "-bound_%d_b", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].B = atof(argv[i+1]); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);

        sprintf(boundaryName, "-bound_%d_dir", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].Dir = atol(argv[i+1]); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);

        sprintf(boundaryName, "-bound_%d_fm", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].Fm = fabs(atof(argv[i+1])); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        if (boundaries.groups[j].Fm<= 0.0) { printf("    ERROR in '%s'!! The force magnitude = '%lf' should not be equal or below 0.0\n", boundaryName, boundaries.groups[j].Fm);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);
      }
      else if (boundaries.groups[j].Type == FORCE_SPEED_ELLIPSE)  // Add forces on an ellipsoid
      {
        sprintf(boundaryName, "-bound_%d_a", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].A = atof(argv[i+1]); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);

        sprintf(boundaryName, "-bound_%d_b", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].B = atof(argv[i+1]); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);

        sprintf(boundaryName, "-bound_%d_dir", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].Dir = atol(argv[i+1]); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);

        sprintf(boundaryName, "-bound_%d_dtheta0", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].Dtheta0 = atof(argv[i+1]); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);

        sprintf(boundaryName, "-bound_%d_dthetadt", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].DthetaDt = atof(argv[i+1]); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);

        sprintf(boundaryName, "-bound_%d_k", j);
        if (strcmp(argv[i],boundaryName) == 0){          
          boundaries.groups[j].K = atof(argv[i+1]); i+=2; }
        else{ printf("    ERROR!! '%s' command not found or is not in right order, However command '%s' is found\n", boundaryName, argv[i]);   exit(EXIT_FAILURE);}
        memset(boundaryName,0,256);
      }
      else if (boundaries.groups[j].Type == DIST_FIXED || boundaries.groups[j].Type == NO_SOLID_BC)         // Fix position
      {

      }
      else{
        printf("    ERROR!! Boundary type %s for boundary %d is undefined for now\n", SOLID_BC_TYPE_CHAR[boundaries.groups[j].Type], j);
        exit(EXIT_FAILURE);}
    }
    i-=1;
    continue;
  }

  if (strcmp(argv[i],"-drag_c") == 0) {                // drag coefficent
    drag.c = atof(argv[i+1]); i++;  continue;}

  // Box Properties:
  //====================================================
  if (ReadBoxParameters(argv[i], argv[i+1])) { i++; continue;}

  // LE Force Properties:
  //====================================================
  if (ReadLinearElasticParameters(argv[i], argv[i+1])) {i++; continue;}

  // Barrier Force Properties:
  //====================================================
  if (ReadBarrierForceParameters(argv[i], argv[i+1])) { i++; continue;}

  // L-Node Dynamics Parameters:
  //====================================================
  if (ReadL0DynamicsParameters(argv[i],argv[i+1])) {i++; continue;}

  // Soft Force Parameters:
  //====================================================
  if (ReadSoftForceParameters(argv[i],argv[i+1])) {i++; continue;}

  // Surface Tension Parameters:
  //====================================================
  if (ReadSurfaceTensionParameters(argv[i],argv[i+1])) {i++; continue;}

  // Volume Conservation Parameters:
  //====================================================
  if (ReadVolumeConservationParameters(argv[i],argv[i+1])) {i++; continue;}

  // Constrain Parameters:
  //====================================================
  if (ReadConstrainParameters(argv[i],argv[i+1])) {i++; continue;}

  // Dynamic Mesh Parameters:
  //====================================================
  if (ReadDynamicMeshParameters(argv[i], argv[i+1])) {i++; continue;}

  // Error
  //====================================================
  printf("    ERROR!! Unknown command '%s'\n", argv[i]);
  exit(EXIT_FAILURE);
}

// Create outputfolder
{
  struct stat st = {0};
  if (stat(vtkOutputFolder, &st) == -1) {
    int check;
    #if __linux == 1
      check = mkdir(vtkOutputFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    #else
      check = mkdir(vtkOutputFolder);  
    #endif 
    if (!check)  printf("Directory '%s' created\n", vtkOutputFolder);
    else {
      printf("Unable to create directory '%s'\n", vtkOutputFolder);
      exit(1);
    }
  }
}

omp_set_num_threads(omp_nthreads);

printf("Open-MP Threads: %d\n", omp_nthreads);
printf("DIM: %d\n", DIM);
printf("VTK Input File: %s\n", vtkInputFile);
printf("Mesh Input Folder: %s\n", meshInputFolder);
printf("CSV Output File: %s\n", csvOutputFile);
printf("VTK Output Folder: %s\n", vtkOutputFolder);
printf("\n");

printf("Time Parameters:\n");
printf("-----------------------------------\n");
printf("  dt:\t\t\t%lf\n", times.dt);
printf("  start_time:\t\t%lf\n", times.real_time);
printf("  end_time:\t\t%lf\n", times.end_time);
printf("  dtime_save_vtk:\t%lf\n", times.dtime_save_vtk);
printf("  dtime_save_csv:\t%lf\n", times.dtime_save_csv);
printf("  dtime_print:\t\t%lf\n", times.dtime_print);
printf("  reload:\t\t%s\n", times.reload ? "true" : "false");
printf("  reload_time:\t\t%lf\n", times.reload_time);
printf("  reload_step:\t\t%d\n", times.reload_step);
printf("\n");

printf("Boundary Parameters:\n");
printf("-----------------------------------\n");
printf("  Number of boundaries:\t\t%d\n", boundaries.N);
for (int j =0; j<boundaries.N; j++)
{
  printf("    Boundary Group:\t\t%d\n", j);
  printf("      Number of Points:\t\t%d\n", boundaries.groups[j].N);
  printf("      Type:\t\t\t%s\n", SOLID_BC_TYPE_CHAR[boundaries.groups[j].Type]);
  if (boundaries.groups[j].Type==FORCE_ADD)
    printf("      Added Force:  (%lf, %lf, %lf) \n", boundaries.groups[j].Fx, boundaries.groups[j].Fy, boundaries.groups[j].Fz);
  else if (boundaries.groups[j].Type==FORCE_ADD_ELLIPSE)
    printf("      Added Force on Ellipse:  (A=%lf, B=%lf, Dir=%d, and Mag=%lf) \n", 
      boundaries.groups[j].A, boundaries.groups[j].B, boundaries.groups[j].Dir, boundaries.groups[j].Fm);
  else if (boundaries.groups[j].Type==FORCE_SPEED_ELLIPSE)
    printf("      Added Speed Force on Ellipse:  (A=%lf, B=%lf, Dir=%d, Dtheta0=%lf, DthetaDt=%lf, and K=%lf) \n", 
      boundaries.groups[j].A, boundaries.groups[j].B, boundaries.groups[j].Dir, 
      boundaries.groups[j].Dtheta0, boundaries.groups[j].DthetaDt, boundaries.groups[j].K);
}
printf("\n");


printf("Drag Parameters:\n");
printf("-----------------------------------\n");
printf("  Drag coefficent:\t\t%lf\n", drag.c);
printf("\n");

PrintBoxProperties();
PrintLinearElasticForceProperties();
PrintBarrierForceProperties();
PrintSoftForceProperties();
PrintSurfaceTensionProperties();
PrintVolumeConservationProperties();
PrintL0DynamicsProperties();
PrintConstrainProperties();
PrintDynamicMeshProperties();

printf("Finished reading input parameters\n");
printf("-----------------------------------\n");

printf("\n");
printf("\n");