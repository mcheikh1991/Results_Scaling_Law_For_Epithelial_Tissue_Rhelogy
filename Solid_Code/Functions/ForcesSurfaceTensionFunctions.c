void DefaultSurfaceTensionParameters()
{
  ST.Type = FORCE_ST_NULL;
  ST.gamma = 0.0;
}

bool ReadSurfaceTensionParameters(const char* str_argv1, const char* str_argv2)
{
  if (strcmp(str_argv1,"-st_type") == 0){    
    ST.Type = (FORCE_ST_TYPE) atol(str_argv2);
    return true;}

  if (strcmp(str_argv1,"-st_gamma") == 0){    
    ST.gamma = atof(str_argv2); 
    return true;}

  return false;
}


void PrintSurfaceTensionProperties()
{
  printf("Surface Tension Properties:\n");
  printf("-----------------------------------\n");
  printf("  ST.Type:\t\t%s\n", FORCE_ST_CHAR[ST.Type]);
  if (ST.Type != FORCE_ST_NULL)
    printf("  ST.gamma:\t\t\t%lf\n", ST.gamma);
  if (ST.Type == FORCE_ST_MULTI)
    printf("  ST.M:  From ST-M.dat file\n");
  printf("\n");

  printf("Volume Force Properties:\n");
  printf("-----------------------------------\n");
  printf("  VC.Type:\t\t%s\n", FORCE_VC_CHAR[VC.Type]);
  if (VC.Type != FORCE_VC_NULL)
  {
    if (VC.Type == FORCE_VC_CONST_CYLIN || VC.Type == FORCE_VC_CONST_CYLIN_CELLS)
      printf("  VC.alpha_cylin:\t\t\t%lf\n", VC.alpha_cylin);
    if (VC.Type == FORCE_VC_CONST_CELLS || VC.Type == FORCE_VC_CONST_CYLIN_CELLS)
      printf("  VC.alpha_cells:\t\t\t%lf\n", VC.alpha_cells);
    printf("  VC.center:\t\t\t(%lf, %lf, %lf)\n", VC.center[0], VC.center[1], VC.center[2]);
  }
  printf("\n");
}

void ReadSurfaceTensionMeshData()
{
  char fileToRead[240];
  if (ST.Type == FORCE_ST_MULTI)
  {
    printf("  Reading the ST-M for each line\n");
    ST.M_line = (double *)malloc(lines.N*sizeof(double));

    sprintf(fileToRead, "./%s/st-m.dat", meshInputFolder);
    ReadDoubleFile(fileToRead, lines.N, ST.M_line);
  }
  else if (ST.Type == FORCE_ST_CELLS)
  {
    ST.Cell_Index  = (int *)malloc(triangles.N*sizeof(int));
    sprintf(fileToRead, "./%s/triangle_cell_index.dat", meshInputFolder);
    ReadIntFile(fileToRead, triangles.N, ST.Cell_Index);
  }
}

//=======================================================================
void SurfaceTensionConst()
{ 
  #pragma omp parallel
  {
    const double gamma = ST.gamma; // Surface Tension const
    double L, nX, nY, nZ;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();

    double *fx_loc = (double *)calloc(forces.N,sizeof(double));
    double *fy_loc = (double *)calloc(forces.N,sizeof(double));
    double *fz_loc = (double *)calloc(forces.N,sizeof(double));
    
    GetLimitsOMPLoop(lines.N, omp_nthreads, omp_rank, &j_start, &j_end);
    for(j=j_start; j<j_end; j++)
    {
      const int p0 = lines.I0[j];
      const int p1 = lines.I1[j];
      CalcNormalVector(p0, p1, &nX, &nY, &nZ, &L);

      fx_loc[p0] +=  gamma*nX ;
      fy_loc[p0] +=  gamma*nY ;
      fz_loc[p0] +=  gamma*nZ*(DIM-2) ;

      fx_loc[p1] += -gamma*nX ;
      fy_loc[p1] += -gamma*nY ;
      fz_loc[p1] += -gamma*nZ*(DIM-2) ;
    }

    #pragma omp critical
    {
      for(int p=0; p<points.N; p++)
        {
          forces.x[p] += fx_loc[p];
          forces.y[p] += fy_loc[p]; 
          forces.z[p] += fz_loc[p]; 
        }
    }

    free(fx_loc);
    free(fy_loc);
    free(fz_loc);
  }
}


void SurfaceTensionMulti()
{ 
  #pragma omp parallel
  {
    const double gamma = ST.gamma; // Surface Tension const
    double L, nX, nY, nZ;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();

    double *fx_loc = (double *)calloc(forces.N,sizeof(double));
    double *fy_loc = (double *)calloc(forces.N,sizeof(double));
    double *fz_loc = (double *)calloc(forces.N,sizeof(double));
    
    GetLimitsOMPLoop(lines.N, omp_nthreads, omp_rank, &j_start, &j_end);
    for(j=j_start; j<j_end; j++)
    {
      const int p0 = lines.I0[j];
      const int p1 = lines.I1[j];
      const double M = ST.M_line[j];
      CalcNormalVector(p0, p1, &nX, &nY, &nZ, &L);

      fx_loc[p0] +=  gamma*M*nX ;
      fy_loc[p0] +=  gamma*M*nY ;
      fz_loc[p0] +=  gamma*M*nZ*(DIM-2) ;

      fx_loc[p1] += -gamma*M*nX ;
      fy_loc[p1] += -gamma*M*nY ;
      fz_loc[p1] += -gamma*M*nZ*(DIM-2) ;
    }

    #pragma omp critical
    {
      for(int p=0; p<points.N; p++)
        {
          forces.x[p] += fx_loc[p];
          forces.y[p] += fy_loc[p]; 
          forces.z[p] += fz_loc[p]; 
        }
    }

    free(fx_loc);
    free(fy_loc);
    free(fz_loc);
  }
}



void SurfaceTensionInitializeCells()
{
  int I_max = ST.Cell_Index[0], I_min = ST.Cell_Index[0];
  for(int t=0; t<triangles.N; t++)
  {
    if (ST.Cell_Index[t] > I_max) I_max = ST.Cell_Index[t];
    if (ST.Cell_Index[t] < I_min) I_min = ST.Cell_Index[t];
  }

  // Number of Cells
  ST.N_Cells = I_max - I_min + 1;

  ST.A  = (double *)calloc(ST.N_Cells, sizeof(double));

  // Calculating Initial Area A0 
  ST.A0 = (double *)calloc(ST.N_Cells, sizeof(double));


  #pragma omp parallel
  {
    int t, t_start, t_end, omp_rank = (int) omp_get_thread_num();

    double *A0_loc = (double *)calloc(ST.N_Cells, sizeof(double));

    GetLimitsOMPLoop(triangles.N, omp_nthreads, omp_rank, &t_start, &t_end);
    for(t=t_start; t<t_end; t++)
    {
      const int c = ST.Cell_Index[t]; // The cell index of the triangle 't'
      A0_loc[c] += CalcTriangleArea2(t);
    }

    #pragma omp critical
    {
      for(int c=0; c<ST.N_Cells; c++)
        ST.A0[c] += A0_loc[c];
    }

    free(A0_loc);
  }
}


void SurfaceTensionAreaCells()
{
  //Reset the values in A to '0'
  memset(ST.A, 0, ST.N_Cells * sizeof(double));   

  // Calculate new Area of cells
  #pragma omp parallel
  {
    int t, t_start, t_end, omp_rank = (int) omp_get_thread_num();

    double *A_loc = (double *)calloc(ST.N_Cells, sizeof(double));

    GetLimitsOMPLoop(triangles.N, omp_nthreads, omp_rank, &t_start, &t_end);
    for(t=t_start; t<t_end; t++)
    {
      const int c = ST.Cell_Index[t]; // The cell index of the triangle 't'
      A_loc[c] += CalcTriangleArea2(t);
    }

    #pragma omp critical
    {
      for(int c=0; c<ST.N_Cells; c++)
        ST.A[c] += A_loc[c];
    }

    free(A_loc);
  }


  // Calculate the new surface tension force
  #pragma omp parallel
  {
    int t, t_start, t_end, omp_rank = (int) omp_get_thread_num();
    double Mx, My, Mz; // Triangle Center
    double nVx, nVy, nVz, V_mag; // Vector

    double *fx_loc = (double *)calloc(forces.N,sizeof(double));
    double *fy_loc = (double *)calloc(forces.N,sizeof(double));
    double *fz_loc = (double *)calloc(forces.N,sizeof(double));

    GetLimitsOMPLoop(triangles.N, omp_nthreads, omp_rank, &t_start, &t_end);
    for(t=t_start; t<t_end; t++)
    {
      const int p0 = triangles.I0[t];
      const int p1 = triangles.I1[t];
      const int p2 = triangles.I2[t];
      const int c = ST.Cell_Index[t]; // The cell index of the triangle 'i'

      const double ST_Mag = ST.gamma*(ST.A[c]-ST.A0[c])/ST.A[c];
      CalcTriangleCenter(t, &Mx, &My, &Mz);

      CalcNormalVector2( Mx, My, Mz, points.x[p0], points.y[p0], points.z[p0], &nVx, &nVy, &nVz, &V_mag);
      fx_loc[p0] += ST_Mag*nVx;
      fy_loc[p0] += ST_Mag*nVy;
      fz_loc[p0] += ST_Mag*nVz;

      CalcNormalVector2( Mx, My, Mz, points.x[p1], points.y[p1], points.z[p1], &nVx, &nVy, &nVz, &V_mag);
      fx_loc[p1] += ST_Mag*nVx;
      fy_loc[p1] += ST_Mag*nVy;
      fz_loc[p1] += ST_Mag*nVz;

      CalcNormalVector2( Mx, My, Mz, points.x[p2], points.y[p2], points.z[p2], &nVx, &nVy, &nVz, &V_mag);
      fx_loc[p2] += ST_Mag*nVx;
      fy_loc[p2] += ST_Mag*nVy;
      fz_loc[p2] += ST_Mag*nVz;
    }

    #pragma omp critical
    {
      for(int p=0; p<points.N; p++)
        {
          forces.x[p] += fx_loc[p];
          forces.y[p] += fy_loc[p]; 
          forces.z[p] += fz_loc[p]; 
        }
    }

    free(fx_loc);
    free(fy_loc);
    free(fz_loc);
  }
}

void ChooseSurfaceTensionFunc()
{
  printf("Choosing the surface tension force function\n"); 

  if (ST.Type == FORCE_ST_NULL){
    printf("  No surface tension force function is choosen\n"); 
    ST.Func = CalculateNullForce;}
  else if (ST.Type == FORCE_ST_CONST){
    printf("  Chosen surface tension force function is for edges with constant gammas\n"); 
    ST.Func = SurfaceTensionConst;}
  else if (ST.Type == FORCE_ST_MULTI){
    printf("  Chosen surface tension force function is for edges with multiple gammas\n"); 
    ST.Func = SurfaceTensionMulti;}
  else if (ST.Type == FORCE_ST_CELLS){
    printf("  Chosen surface tension force function is for cells made of triangules\n"); 
    SurfaceTensionInitializeCells();
    ST.Func = SurfaceTensionAreaCells;}
  else{
    printf("  ERROR!! Surface tension function is undefined\n");
    exit(EXIT_FAILURE);}
}


void ForcesAddSurfaceTensionForces(){ ST.Func();}
