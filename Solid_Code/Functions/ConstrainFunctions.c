void DefaultConstrainParameters()
{
  CN.Type = CONSTRAIN_NULL;
}

void PrintConstrainProperties()
{
  printf("Constrain Properties:\n");
  printf("-----------------------------------\n");
  printf("  CN.Type:\t\t%s\n", CONSTRAIN_CHAR[CN.Type]);
  if (CN.Type != CONSTRAIN_NULL)
  {
    if (CN.Type == CONSTRAIN_AXIS || CN.Type == CONSTRAIN_TWO_AXES)
    {
      printf("  CN.axis_i:\t\t%d\n", CN.axis_i);
      printf("  CN.axis_v:\t\t%lf\n", CN.axis_v);
      if (CN.Type == CONSTRAIN_TWO_AXES)
      {
        printf("  CN.axis_i2:\t\t%d\n", CN.axis_i2);
        printf("  CN.axis_v2:\t\t%lf\n", CN.axis_v2);
      }
    }

  }
  printf("\n");
}

bool ReadConstrainParameters(const char* str_argv1, const char* str_argv2)
{
  if (strcmp(str_argv1,"-cn_type") == 0){    
    CN.Type = (CONSTRAIN_TYPE) atol(str_argv2); 
    return true;}
    
  if (CN.Type == CONSTRAIN_AXIS || CN.Type == CONSTRAIN_TWO_AXES)
  {
    if (strcmp(str_argv1,"-cn_axis_i") == 0){    // index of axis (0=x, 1=y, 2=z)
      CN.axis_i = atol(str_argv2); 
      return true;}

    if (strcmp(str_argv1,"-cn_axis_v") == 0){    // value of axis
      CN.axis_v = atof(str_argv2); 
      return true;}

    if (CN.Type == CONSTRAIN_TWO_AXES)
    {
      if (strcmp(str_argv1,"-cn_axis_i2") == 0){    // index of axis (0=x, 1=y, 2=z)
        CN.axis_i2 = atol(str_argv2); 
        return true;}

      if (strcmp(str_argv1,"-cn_axis_v2") == 0){    // value of axis
        CN.axis_v2 = atof(str_argv2); 
        return true;}
    }

    if (strcmp(str_argv1,"-cn_ds") == 0){    // value of axis
      CN.ds = atof(str_argv2); 
      return true;}
  }
  else if (CN.Type!= CONSTRAIN_NULL){
    printf("    ERROR in ReadConstrainParameters!! Only constrain type '%s' and '%s' are defined for now not '%s'\n", CONSTRAIN_CHAR[3], CONSTRAIN_CHAR[4], CONSTRAIN_CHAR[CN.Type]);
    exit(EXIT_FAILURE);}

  return false;
}



// Single Axis Constrain
//------------------------------------------------------------------------------------

void ConstrainCorrectPointsOnAxis()
{
  /* A method to constrain points on chosen axis*/

  #pragma omp parallel
  {
    const double V = CN.axis_v;
    const int axis = CN.axis_i;
    const int N    = CN.N;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();

    GetLimitsOMPLoop(N, omp_nthreads, omp_rank, &j_start, &j_end);
    for(j=j_start; j<j_end; j++)
    {
      const int p = CN.p[j]; 

      if (axis==0)       points.x[p] = V;
      else if  (axis==1) points.y[p] = V;
      else if  (axis==2) points.z[p] = V;
    }
  }
}

void ConstrainCorrectForcesOnAxis()
{
  /* A method to constrain the forces on constrained points*/

  #pragma omp parallel
  {
    const int axis = CN.axis_i;
    const int N    = CN.N;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();

    GetLimitsOMPLoop(N, omp_nthreads, omp_rank, &j_start, &j_end);
    for(j=j_start; j<j_end; j++)
    {
      const int p = CN.p[j];

      if (axis==0)       forces.x[p] = 0.0;
      else if  (axis==1) forces.y[p] = 0.0;
      else if  (axis==2) forces.z[p] = 0.0;
    }
  }
}



void ConstrainCreateAxis()
{
  /* A function that finds the point to constrain with respect to the axis */
  // We remove ds from the axis as a tolerance
  const double ds    = CN.ds;
  const double Vmin  = CN.axis_v - ds;
  const double Vmax  = CN.axis_v + ds;
  const int    axis  = CN.axis_i;

  if (axis==0)      printf("  Looking for points to be constrainted on x:%g with tol:%g\n", (double) CN.axis_v, (double) ds);
  else if (axis==1) printf("  Looking for points to be constrainted on y:%g with tol:%g\n", (double) CN.axis_v, (double) ds);
  else if (axis==2) printf("  Looking for points to be constrainted on z:%g with tol:%g\n", (double) CN.axis_v, (double) ds);
  else{ 
    printf("    ERROR in ConstrainCreateAxis!! Incorrect index %d. It should either be 0 for x, 1 for y, and 2 for z\n", axis);
    exit(EXIT_FAILURE);}

  // Counting the number of points
  CN.N = 0;
  for (int p=0; p<points.N; p++)
  {
    const double x = points.x[p];
    const double y = points.y[p];
    const double z = points.z[p];
    if (axis==0 && x>=Vmin && x<=Vmax)       CN.N++;
    else if (axis==1 && y>=Vmin && y<=Vmax)  CN.N++;
    else if (axis==2 && z>=Vmin && z<=Vmax)  CN.N++;
  }
  CN.p = (int *)malloc(CN.N*sizeof(int));

  // Saving the points on the constrain
  int j = 0;
  for (int p=0; p<points.N; p++)
  {
    const double x = points.x[p];
    const double y = points.y[p];
    const double z = points.z[p];
    if (axis==0 && x>=Vmin && x<=Vmax){       CN.p[j] = p;    j++; }
    else if (axis==1 && y>=Vmin && y<=Vmax){  CN.p[j] = p;    j++; }
    else if (axis==2 && z>=Vmin && z<=Vmax){  CN.p[j] = p;    j++; }
  }
  printf("    The number of constrained points is %d\n",  CN.N); 
  ConstrainCorrectPointsOnAxis();
}

void ConstrainFreeAxis() {  free(CN.p); }


// Two Axis Constrain
//------------------------------------------------------------------------------------

void ConstrainCorrectPointsOnTwoAxes()
{
  /* A method to constrain points on chosen axis*/

  #pragma omp parallel
  {
    const double V = CN.axis_v;
    const int axis = CN.axis_i;
    const int N    = CN.N;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();

    GetLimitsOMPLoop(N, omp_nthreads, omp_rank, &j_start, &j_end);
    for (j=j_start; j<j_end; j++)
    {
      const int p = CN.p[j];

      if (axis==0)       points.x[p] = V;
      else if  (axis==1) points.y[p] = V;
      else if  (axis==2) points.z[p] = V;
    }
  }
  
  #pragma omp parallel
  {
    const double V2 = CN.axis_v2;
    const int axis2 = CN.axis_i2;
    const int N2    = CN.N2;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();

    GetLimitsOMPLoop(N2, omp_nthreads, omp_rank, &j_start, &j_end);
    for (j=j_start; j<j_end; j++)
    {
      const int p2 = CN.p2[j];

      if (axis2==0)       points.x[p2] = V2;
      else if  (axis2==1) points.y[p2] = V2;
      else if  (axis2==2) points.z[p2] = V2;
    }
  }
}


void ConstrainCorrectForcesOnTwoAxes()
{
  /* A method to constrain the forces on constrained points*/
  #pragma omp parallel
  {
    const int axis = CN.axis_i;
    const int N    = CN.N;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(N, omp_nthreads, omp_rank, &j_start, &j_end);
    for (j=j_start; j<j_end; j++)
    {
      const int p = CN.p[j];

      if (axis==0)       forces.x[p] = 0.0;
      else if  (axis==1) forces.y[p] = 0.0;
      else if  (axis==2) forces.z[p] = 0.0;
    }
  }

  #pragma omp parallel
  {  
    const int axis2 = CN.axis_i2;
    const int N2    = CN.N2;
    int j, j_start, j_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(N2, omp_nthreads, omp_rank, &j_start, &j_end);
    for (j=j_start; j<j_end; j++)
    {
      const int p2 = CN.p2[j];

      if (axis2==0)       forces.x[p2] = 0.0;
      else if  (axis2==1) forces.y[p2] = 0.0;
      else if  (axis2==2) forces.z[p2] = 0.0;
    }
  }
}

void ConstrainCreateTwoAxes()
{
  /* A function that finds the point to constrain with respect to the axis */
  // We remove ds from the axis as a tolerance
  const double ds    = CN.ds;

  const double Vmin  = CN.axis_v - ds;
  const double Vmax  = CN.axis_v + ds;
  const int    axis  = CN.axis_i;

  const double Vmin2  = CN.axis_v2 - ds;
  const double Vmax2  = CN.axis_v2 + ds;
  const int    axis2  = CN.axis_i2;

  // Counting the number of points
  CN.N = 0, CN.N2 = 0;
  for (int p=0; p<points.N; p++)
  {
    const double x = points.x[p];
    const double y = points.y[p];
    const double z = points.z[p];
    if ((axis==0 && x>=Vmin && x<=Vmax) || (axis==1 && y>=Vmin && y<=Vmax) || (axis==2 && z>=Vmin && z<=Vmax))
        CN.N++;
    else if ((axis2==0 && x>=Vmin2 && x<=Vmax2) || (axis2==1 && y>=Vmin2 && y<=Vmax2)|| (axis2==2 && z>=Vmin2 && z<=Vmax2))
        CN.N2++;
  }
  CN.p = (int *)malloc(CN.N*sizeof(int));
  CN.p2 = (int *)malloc(CN.N2*sizeof(int));

  // Saving the points on the constrain
  int j1 = 0, j2 = 0;
  for (int p=0; p<points.N; p++)
  {
    const double x = points.x[p];
    const double y = points.y[p];
    const double z = points.z[p];
    if     ((axis==0 && x>=Vmin && x<=Vmax) || (axis==1 && y>=Vmin && y<=Vmax) || (axis==2 && z>=Vmin && z<=Vmax)){
      CN.p[j1] = p;  j1++; }
    else if ((axis2==0 && x>=Vmin2 && x<=Vmax2) || (axis2==1 && y>=Vmin2 && y<=Vmax2)|| (axis2==2 && z>=Vmin2 && z<=Vmax2)){
      CN.p2[j2] = p; j2++; }
  }

  //printf("    The number of constrained points is %d\n",  CN.N+CN.N2); 
  ConstrainCorrectPointsOnTwoAxes();
}

void ConstrainFreeTwoAxes() {  free(CN.p); free(CN.p2); }

//--------------------------------------
void ConstrainNullFunc(){};

void ChooseConstrainFuncs()
{
  printf("Choosing the constrain function\n"); 
  if (CN.Type == CONSTRAIN_NULL)
  {
    printf("  No constrain point function is choosen\n");
    CN.Func_create = ConstrainNullFunc;
    CN.Func_points = ConstrainNullFunc;
    CN.Func_forces = ConstrainNullFunc;
    CN.Func_free   = ConstrainNullFunc;
  } 
  else if (CN.Type == CONSTRAIN_AXIS) 
  {
    printf("  Constrain point function is chosen for one axis\n"); 
    CN.Func_create = ConstrainCreateAxis;
    CN.Func_points = ConstrainCorrectPointsOnAxis;
    CN.Func_forces = ConstrainCorrectForcesOnAxis;
    CN.Func_free   = ConstrainFreeAxis;
  }
  else if (CN.Type == CONSTRAIN_TWO_AXES) 
  {
    printf("  Constrain point function is chosen for two axis\n"); 
    CN.Func_create = ConstrainCreateTwoAxes;
    CN.Func_points = ConstrainCorrectPointsOnTwoAxes;
    CN.Func_forces = ConstrainCorrectForcesOnTwoAxes;
    CN.Func_free   = ConstrainFreeTwoAxes;
  }
  else
  {
    printf("  ERROR!! Constrain point function '%s' is undefined\n", CONSTRAIN_CHAR[LE.Type]);
    exit(EXIT_FAILURE);
  }
}


void ConstrainCreate(){ CN.Func_create();}
void ConstrainPoints(){ CN.Func_points();}
void ConstrainForces(){ CN.Func_forces();}
void ConstrainFree(){   CN.Func_free();}