void DefaultLinearElasticParameters()
{
  // LE Force
  LE.Type = FORCE_LE_INITIAL_DISP;
  LE.K = 0.0;
}

void ReadLinearElasticMeshData()
{
  char fileToRead[240];
  if (LE.Type == FORCE_LE_INITIAL_DISP_MULTI)
  {
    printf("  Reading the LE-M for each line\n");
    LE.M_line = (double *)malloc(lines.N*sizeof(double));

    sprintf(fileToRead, "./%s/le-m.dat", LE.M_file_name);
    ReadDoubleFile(fileToRead, lines.N, LE.M_line);
  }
}

void PrintLinearElasticForceProperties()
{
  printf("Linear Elastic Force Properties:\n");
  printf("-----------------------------------\n");
  printf("  LE.Type:\t\t%s\n", FORCE_LE_CHAR[LE.Type]);
  if (LE.Type == FORCE_LE_INITIAL_DISP || LE.Type == FORCE_LE_INITIAL_DISP_MULTI)
    printf("  LE.K:\t\t\t%lf\n", LE.K);
  if (LE.Type == FORCE_LE_INITIAL_DISP_MULTI)
    printf("  LE.M:\t\t\tFrom LE-M.dat file\n");
  printf("\n");
}

bool ReadLinearElasticParameters(const char* str_argv1, const char* str_argv2)
{
  if (strcmp(str_argv1,"-le_type") == 0){    
    LE.Type = (FORCE_LE_TYPE) atol(str_argv2); 
    return true;}
    
  if (LE.Type == FORCE_LE_FIXED_DISP || LE.Type == FORCE_LE_INITIAL_DISP || LE.Type == FORCE_LE_INITIAL_DISP_MULTI)
    if (strcmp(str_argv1,"-le_k") == 0){    
      LE.K = atof(str_argv2); 
      return true;}

  if (LE.Type == FORCE_LE_INITIAL_DISP_MULTI)
    if (strcmp(str_argv1,"-le_multiplier_file_loc") == 0) { 
      strcpy(LE.M_file_name, str_argv2); 
      return true;}

  return false;
}

//==============================================================================
void CalculateLinearElastic()
{
  #pragma omp parallel
  {
    const double K = LE.K;
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
      const double L0 = lines.len[j];
      CalcNormalVector(p0, p1, &nX, &nY, &nZ, &L);

      fx_loc[p0] +=  K*(L-L0)*nX ;
      fy_loc[p0] +=  K*(L-L0)*nY ;
      fz_loc[p0] +=  K*(L-L0)*nZ*(DIM-2) ;

      fx_loc[p1] += -K*(L-L0)*nX ;
      fy_loc[p1] += -K*(L-L0)*nY ;
      fz_loc[p1] += -K*(L-L0)*nZ*(DIM-2) ;
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


void CalculateLinearElasticMultiple()
{
  #pragma omp parallel
  {
    double L, nX, nY, nZ;
    const double K = LE.K;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();

    double *fx_loc = (double *)calloc(forces.N,sizeof(double));
    double *fy_loc = (double *)calloc(forces.N,sizeof(double));
    double *fz_loc = (double *)calloc(forces.N,sizeof(double));

    GetLimitsOMPLoop(lines.N, omp_nthreads, omp_rank, &j_start, &j_end);
    for(j=j_start; j<j_end; j++)
    {
      const int p0 = lines.I0[j];
      const int p1 = lines.I1[j];
      const double L0 = lines.len[j];
      const double M = LE.M_line[j];
      CalcNormalVector(p0, p1, &nX, &nY, &nZ, &L);

      fx_loc[p0] +=  K*M*(L-L0)*nX ;
      fy_loc[p0] +=  K*M*(L-L0)*nY ;
      fz_loc[p0] +=  K*M*(L-L0)*nZ*(DIM-2) ;

      fx_loc[p1] += -K*M*(L-L0)*nX ;
      fy_loc[p1] += -K*M*(L-L0)*nY ;
      fz_loc[p1] += -K*M*(L-L0)*nZ*(DIM-2) ;
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

void ChooseLinearElasticFunc()
{
  printf("Choosing the linear elastic force function\n"); 
  if (LE.Type == FORCE_LE_NULL){
    printf("  No linear elastic force function is choosen\n"); 
    LE.Func = CalculateNullForce;}
  else if (LE.Type == FORCE_LE_INITIAL_DISP) {
    printf("  Linear elastic force function is chosen with K=%lf, and L0 from initial length\n", LE.K); 
    LE.Func = CalculateLinearElastic;}
  else if (LE.Type == FORCE_LE_INITIAL_DISP_MULTI) {
    printf("  Linear elastic force function is chosen with multiple K from le-m.dat file, and L0 from initial length\n"); 
    LE.Func = CalculateLinearElasticMultiple;}
  else{
    printf("  ERROR!! Linear elastic force function '%s' is undefined\n", FORCE_LE_CHAR[LE.Type]);
    exit(EXIT_FAILURE);}
}

void ForcesAddLinearElasticForce(){ LE.Func();}
