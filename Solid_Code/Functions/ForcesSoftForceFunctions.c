void  DefaultSoftForceParameters()
{
  // Soft Force
  SF.Type= FORCE_SF_NULL;  
  SF.A = 250;
  SF.B = 77;
  SF.C = 77;
  SF.ds = 2;
  SF.sf = 0.11*50000;
}

void PrintSoftForceProperties()
{
  printf("Soft Force Properties:\n");
  printf("-----------------------------------\n");
  printf("  SF.Type:\t\t%s\n", FORCE_SF_CHAR[SF.Type]);
  if (SF.Type != FORCE_SF_NULL)
  {
    printf("  SF.A:\t\t\t%lf\n", SF.A);
    printf("  SF.B:\t\t\t%lf\n", SF.B);
    printf("  SF.C:\t\t\t%lf\n", SF.C);
    printf("  SF.Bd:\t\t%lf\n", SF.Bd);
    printf("  SF.Bv:\t\t%lf\n", SF.Bv);
    printf("  SF.ds:\t\t%lf\n", SF.ds);
    printf("  SF.sf:\t\t%lf\n", SF.sf);
  }
  printf("\n");
}
bool ReadSoftForceParameters(const char* str_argv1, const char* str_argv2)
{
  if (strcmp(str_argv1,"-sf_type") == 0){    
    SF.Type = (FORCE_SF_TYPE) atol(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-sf_a") == 0){    
    SF.A = atof(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-sf_b") == 0){    
    SF.B = atof(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-sf_c") == 0){    
    SF.C = atof(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-sf_bd") == 0){    
    SF.Bd = atof(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-sf_bv") == 0){    
    SF.Bv = atof(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-sf_ds") == 0){    
    SF.ds = atof(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-sf_sf") == 0){
    SF.sf = atof(str_argv2); 
    return true;}

  return false;
}



//==============================================================================

void SoftForceCreateOneEllipsoid(bool constrain)
{
  /* A functiona the calculates the points on an ellispoid */

  // We remove ds from the axis as a tolerance
  const double ds = SF.ds;
  const double As = SF.A - ds;
  const double Al = SF.A + ds;
  const double Bs = SF.B - ds;
  const double Bl = SF.B + ds;
  const double Cs = SF.C - ds;
  const double Cl = SF.C + ds;

  double Es, El = 0; // Small and Large Ellipsoids

  printf("  Looking for points to be add soft force on ellipsoid A:%g, B:%g, C:%g with tol:%g\n", As+ds, Bs+ds, Cs+ds, ds); 

  // Counting the number of points
  SF.N = 0;
  for (int p=0; p<points.N; p++){
    const double x = points.x[p];
    const double y = points.y[p];
    const double z = points.z[p];

    if (DIM == 3) 
      Es = x*x/(As*As) + y*y/(Bs*Bs) + z*z/(Cs*Cs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl) + z*z/(Cl*Cl);
    else          
      Es = x*x/(As*As) + y*y/(Bs*Bs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl);

    if (Es>=1.0 && El<=1.0) SF.N++;
  }

  SF.p = (int *)calloc(SF.N,sizeof(int));

  // Saving the points to receive softforce
  int j = 0;
  for (int p=0; p<points.N; p++){
    const double x = points.x[p];
    const double y = points.y[p];
    const double z = points.z[p];

    if (DIM == 3) 
      Es = x*x/(As*As) + y*y/(Bs*Bs) + z*z/(Cs*Cs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl) + z*z/(Cl*Cl);
    else          
      Es = x*x/(As*As) + y*y/(Bs*Bs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl);

    if (Es>=1.0 && El<=1.0) {
      SF.p[j] = p;
      j++;}
  }

  printf("  The number of points that will recieve soft force is is %d\n",  SF.N); 

  /* A method to constrain points using vector from origin*/
  if (constrain)
  {
    double Rho;

    const double A  = SF.A;
    const double B  = SF.B;
    const double C  = SF.C;

    for (int j=0; j<SF.N; j++)
    {
      const int p = SF.p[j];

      const double x = points.x[p];
      const double y = points.y[p];
      const double z = points.z[p];

      Rho = sqrt( (x*x)/(A*A) + (y*y)/(B*B) + (z*z)/(C*C) );

      points.x[p] = x/Rho;
      points.y[p] = y/Rho;
      points.z[p] = z/Rho*(DIM-2);
    }
    printf("  Finished constraining %d points on the ellipsoid\n",  SF.N); 
  }
}


void SoftForceCreateTwoEllipsoids(bool constrain)
{
  /* A functiona the calculates the points on an ellispoid */

  // We remove ds from the axis as a tolerance
  const double ds  = SF.ds;
  const double As  = SF.A  - ds;
  const double Al  = SF.A  + ds;
  const double Bds = SF.Bd - ds;
  const double Bdl = SF.Bd + ds;
  const double Bvs = SF.Bv - ds;
  const double Bvl = SF.Bv + ds;
  const double Cs  = SF.C  - ds;
  const double Cl  = SF.C  + ds;
  double Bs, Bl;
  double Es, El = 0; // Small and Large Ellipsoids

  printf("  Looking for points to be add soft force on ellipsoid A:%g, Bd:%g, C:%g (y>0) & ellipsoid A:%g, Bv:%g, C:%g (y<0) with tol:%g\n", 
      (double) As+ds, (double) Bds+ds, (double) Cs+ds, (double) As+ds, (double) Bvs+ds, (double) Cs+ds, (double) ds); 

  // Counting the number of points
  SF.N = 0;
  for (int p=0; p<points.N; p++){
    const double x = points.x[p];
    const double y = points.y[p];
    const double z = points.z[p];

    if (y>0) Bs = Bds, Bl = Bdl;
    else     Bs = Bvs, Bl = Bvl;

    if (DIM == 3) 
      Es = x*x/(As*As) + y*y/(Bs*Bs) + z*z/(Cs*Cs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl) + z*z/(Cl*Cl);
    else          
      Es = x*x/(As*As) + y*y/(Bs*Bs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl);

    if (Es>=1.0 && El<=1.0) SF.N++;
  }

  SF.p = (int *)calloc(SF.N,sizeof(int));

  // Saving the points on the constrain
  int j = 0;
  for (int p=0; p<points.N; p++){
    const double x = points.x[p];
    const double y = points.y[p];
    const double z = points.z[p];

    if (y>0) Bs = Bds, Bl = Bdl;
    else     Bs = Bvs, Bl = Bvl;

    if (DIM == 3) 
      Es = x*x/(As*As) + y*y/(Bs*Bs) + z*z/(Cs*Cs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl) + z*z/(Cl*Cl);
    else          
      Es = x*x/(As*As) + y*y/(Bs*Bs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl);

    if (Es>=1.0 && El<=1.0) {
      SF.p[j] = p;
      j++;}
  }

  printf("  The number of points that will recive soft force is %d\n",  SF.N); 

  if (constrain)
  {
    double Rho;

    const double A  = SF.A;
    const double Bd = SF.Bd;
    const double Bv = SF.Bv;
    const double C  = SF.C;
    double B;

    for (int j=0; j<SF.N; j++)
    {
      const int p = SF.p[j];

      const double x = points.x[p];
      const double y = points.y[p];
      const double z = points.z[p];

      if (y>0) B = Bd;
      else     B = Bv;

      Rho = sqrt( (x*x)/(A*A) + (y*y)/(B*B) + (z*z)/(C*C) );

      points.x[p] = x/Rho;
      points.y[p] = y/Rho;
      points.z[p] = z/Rho*(DIM-2);

    }
    printf("  Finished constraining %d points on two ellipsoids\n",  SF.N); 
  }
}

void SoftForceOnOneEllipsoid()
{
  #pragma omp parallel
  {    
    const double A = SF.A, B = SF.B, C = SF.C, sf = SF.sf;
    double E, dEdx, dEdy, dEdz;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();

    double *fx_loc = (double *)calloc(forces.N,sizeof(double));
    double *fy_loc = (double *)calloc(forces.N,sizeof(double));
    double *fz_loc = (double *)calloc(forces.N,sizeof(double));

    GetLimitsOMPLoop(SF.N, omp_nthreads, omp_rank, &j_start, &j_end);
    for (j=j_start; j<j_end; j++)
    {
      const int p = SF.p[j];

      const double x = points.x[p];
      const double y = points.y[p];
      const double z = points.z[p];

      E = -1.0 + x*x/(A*A) + y*y/(B*B) + (DIM-2)*z*z/(C*C); // Ellipsoid Equation 

      //Gradient of Ellipsoid Equation
      dEdx = 2.0*x/(A*A);
      dEdy = 2.0*y/(B*B);
      dEdz = 2.0*z/(C*C);

      fx_loc[p] += -sf*E*dEdx;
      fy_loc[p] += -sf*E*dEdy;
      fz_loc[p] += -sf*E*dEdz;
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

void SoftForceOnOneEllipsoidOutside()
{
  /* A method to add soft force*/
  #pragma omp parallel
  { 
    const double A = SF.A, B = SF.B, C = SF.C, sf = SF.sf;
    double E, dEdx, dEdy, dEdz;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();

    double *fx_loc = (double *)calloc(forces.N,sizeof(double));
    double *fy_loc = (double *)calloc(forces.N,sizeof(double));
    double *fz_loc = (double *)calloc(forces.N,sizeof(double));

    GetLimitsOMPLoop(SF.N, omp_nthreads, omp_rank, &j_start, &j_end);
    for (j=j_start; j<j_end; j++)
    {
      const int p = SF.p[j];

      const double x = points.x[p];
      const double y = points.y[p];
      const double z = points.z[p];

      E = -1.0 + x*x/(A*A) + y*y/(B*B) + (DIM-2)*z*z/(C*C); // Ellipsoid Equation 

      if (E>=0) // Point outside
      {
        //Gradient of Ellipsoid Equation
        dEdx = 2.0*x/(A*A);
        dEdy = 2.0*y/(B*B);
        dEdz = 2.0*z/(C*C);

        fx_loc[p] += -sf*E*dEdx;
        fy_loc[p] += -sf*E*dEdy;
        fz_loc[p] += -sf*E*dEdz;
      }
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

void SoftForceOnTwoEllipsoids()
{
  #pragma omp parallel
  {   
    /* A method to add soft force*/
    const double A = SF.A, Bd = SF.Bd, Bv = SF.Bv, C = SF.C, sf = SF.sf;
    double E, dEdx, dEdy, dEdz, B;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();

    double *fx_loc = (double *)calloc(forces.N,sizeof(double));
    double *fy_loc = (double *)calloc(forces.N,sizeof(double));
    double *fz_loc = (double *)calloc(forces.N,sizeof(double));

    GetLimitsOMPLoop(SF.N, omp_nthreads, omp_rank, &j_start, &j_end);
    for (j=j_start; j<j_end; j++)
    {
      const int p = SF.p[j];

      const double x = points.x[p];
      const double y = points.y[p];
      const double z = points.z[p];

      if (y>0) B = Bd;
      else     B = Bv;
      
      E = -1.0 + x*x/(A*A) + y*y/(B*B) + (DIM-2)*z*z/(C*C); // Ellipsoid Equation 

      //Gradient of Ellipsoid Equation
      dEdx = 2.0*x/(A*A);
      dEdy = 2.0*y/(B*B);
      dEdz = 2.0*z/(C*C);

      fx_loc[p] += -sf*E*dEdx;
      fy_loc[p] += -sf*E*dEdy;
      fz_loc[p] += -sf*E*dEdz;
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

void ChooseSoftForcesFunc(bool print_sf_points)
{
  printf("Choosing the soft force force function\n"); 

  if (SF.Type == FORCE_SF_NULL){
    printf("  No soft force function is choosen\n"); 
    SF.Func = CalculateNullForce;}
  else if (SF.Type == FORCE_SF_ELLIPSOID){
    printf("  Chosen soft force function is for an ellipsoid\n"); 
    SoftForceCreateOneEllipsoid(true);
    SF.Func = SoftForceOnOneEllipsoid;}
  else if (SF.Type == FORCE_SF_TWO_ELLIPSOIDS){
    printf("  Chosen soft force function is for two ellipsoids\n"); 
    SoftForceCreateTwoEllipsoids(true);
    SF.Func = SoftForceOnTwoEllipsoids;}
  else if (SF.Type == FORCE_SF_ELLIPSOID_OUTSIDE){
    printf("  Chosen soft force function is for an ellipsoid but only if points are out\n"); 
    SoftForceCreateOneEllipsoid(false);
    SF.Func = SoftForceOnOneEllipsoidOutside;}
  else{
    printf("  ERROR!! soft force function is undefined\n");
    exit(EXIT_FAILURE);}

  if((SF.Type != FORCE_SF_NULL) && (print_sf_points))
    WriteVTKPointFile("SF_Points", SF.N, SF.p);
}


void ForcesAddSoftForces(){ SF.Func();}

