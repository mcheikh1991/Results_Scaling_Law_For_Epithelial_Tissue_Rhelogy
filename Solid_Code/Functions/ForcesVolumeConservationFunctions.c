
inline double CalculatPartialVolumeFromTrianlge(int t);
double CalculateVolumeOfCylinder();
double CalculateTotalVolumeOfCells();
void CalculateVolumeOfCells();
void CalculateCentroidOfCells();

void DefaultVolumeConservationParameters()
{
  VC.Type = FORCE_VC_NULL;
  VC.alpha_cylin = 0.0;
  VC.alpha_cells = 0.0;
  VC.volume = 0.0; VC.volume_0 = 0.0;   VC.volume_cells = 0.0;
  VC.center[0] = 0.0,       VC.center[1] = 0.0,       VC.center[2] = 0.0;
  VC.cylin_circle1_point[0] = 0.0,  VC.cylin_circle1_point[1] = 0.0,  VC.cylin_circle1_point[2] = 0.0;
  VC.cylin_circle1_normal[0] = 0.0, VC.cylin_circle1_normal[1] = 0.0, VC.cylin_circle1_normal[2] = 0.0;
  VC.cylin_circle2_point[0] = 0.0,  VC.cylin_circle2_point[1] = 0.0,  VC.cylin_circle2_point[2] = 0.0;
  VC.cylin_circle2_normal[0] = 0.0, VC.cylin_circle2_normal[1] = 0.0, VC.cylin_circle2_normal[2] = 0.0;
  VC.N_cells = 0; VC.N_cylin_triangles=0; VC.N_cylin_points = 0;

}

void PrintVolumeConservationProperties()
{
  printf("Volume Conservation Properties:\n");
  printf("-----------------------------------\n");
  printf("  VC.Type:\t\t%s\n", FORCE_VC_CHAR[VC.Type]);
  printf("\n");
}

bool ReadVolumeConservationParameters(const char* str_argv1, const char* str_argv2)
{
  if (strcmp(str_argv1,"-vc_type") == 0){    
    VC.Type = (FORCE_VC_TYPE) atol(str_argv2);
    return true;}

  if (strcmp(str_argv1,"-vc_alpha_cylin") == 0){    
    VC.alpha_cylin = atof(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-vc_alpha_cells") == 0){    
    VC.alpha_cells = atof(str_argv2); 
    return true;}

  // center
  if (strcmp(str_argv1,"-vc_center_x") == 0){    
    VC.center[0] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_center_y") == 0){    
    VC.center[1] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_center_z") == 0){    
    VC.center[2] = atof(str_argv2); return true;}

  // First circle face parameters
  if (strcmp(str_argv1,"-vc_cylin_p1_x") == 0){    
    VC.cylin_circle1_point[0] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_cylin_p1_y") == 0){    
    VC.cylin_circle1_point[1] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_cylin_p1_z") == 0){    
    VC.cylin_circle1_point[2] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_cylin_n1_x") == 0){    
    VC.cylin_circle1_normal[0] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_cylin_n1_y") == 0){    
    VC.cylin_circle1_normal[1] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_cylin_n1_z") == 0){    
    VC.cylin_circle1_normal[2] = atof(str_argv2); return true;}

  // Second circle face parameters
  if (strcmp(str_argv1,"-vc_cylin_p2_x") == 0){    
    VC.cylin_circle2_point[0] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_cylin_p2_y") == 0){    
    VC.cylin_circle2_point[1] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_cylin_p2_z") == 0){    
    VC.cylin_circle2_point[2] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_cylin_n2_x") == 0){    
    VC.cylin_circle2_normal[0] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_cylin_n2_y") == 0){    
    VC.cylin_circle2_normal[1] = atof(str_argv2); return true;}

  if (strcmp(str_argv1,"-vc_cylin_n2_z") == 0){    
    VC.cylin_circle2_normal[2] = atof(str_argv2); return true;}

  return false;
}

void ReadVolumeConservationMeshData()
{
  char fileToRead[240];
  if (VC.Type == FORCE_VC_CONST_CYLIN || VC.Type == FORCE_VC_CONST_CYLIN_CELLS)
  {
    // Part 1: Read the index of triangles that will define the triangles to be used in cylinder volume calculation 
    sprintf(fileToRead, "./%s/vc-cylin-triangles.dat", meshInputFolder);
    VC.N_cylin_triangles = GetFileSize(fileToRead);
    printf("Number of triangles for cylinder volume conservation: %d\n", VC.N_cylin_triangles);

    VC.cylin_triangles = (int *)malloc(VC.N_cylin_triangles*sizeof(int));
    sprintf(fileToRead, "./%s/vc-cylin-triangles.dat", meshInputFolder);
    ReadIntFile(fileToRead, VC.N_cylin_triangles, VC.cylin_triangles);

    // Part 2: Fix orientation of the triangular normal
    {
      double Mx, My, Mz, nVx, nVy, nVz;            
      for(int ct=0; ct<VC.N_cylin_triangles; ct++)
      {
        const int t = VC.cylin_triangles[ct];
        CalcTriangleCenter(t, &Mx, &My, &Mz);
        CalcTriangleNormal(t, &nVx, &nVy, &nVz);  
        if (( nVx*(Mx-VC.center[0]) + nVy*(My-VC.center[1]) + nVz*(Mz-VC.center[2]) ) < 0)
        {
          const int swap_index = triangles.I1[t];
          triangles.I1[t] = triangles.I2[t];
          triangles.I2[t] = swap_index;
        }
      }
    }

    // Part 3: Reads the index of points that will define the outer circle of the cylinder
    sprintf(fileToRead, "./%s/vc-cylin-points.dat", meshInputFolder);
    VC.N_cylin_points = GetFileSize(fileToRead);
    printf("Number of points for cylinder volume conservation: %d\n", VC.N_cylin_points);

    VC.cylin_points = (int *)malloc(VC.N_cylin_points*sizeof(int));
    sprintf(fileToRead, "./%s/vc-cylin-points.dat", meshInputFolder);
    ReadIntFile(fileToRead, VC.N_cylin_points, VC.cylin_points);

    // Part 4: Now we want to re-order the points in a counter clockwise direction
    printf("  Re-Orienting the cylinder surface points based on the value of theta to be in counter-clockwise direction\n"); 

    /* Re-Orient a solid boundary based on the value of theta. */
    {
      bool  Found;
      double *Theta, *Theta_old;
      int *old_bound_pts, *new_bound_pts;
      double R;

      if(VC.N_cylin_points<1){ printf("  ERROR!! Volume conservartion function is undefined\n");     exit(EXIT_FAILURE);}

      // Calculating the Theta for each point
      Theta         = (double *)malloc(VC.N_cylin_points*sizeof(double));
      Theta_old     = (double *)malloc(VC.N_cylin_points*sizeof(double));
      new_bound_pts =    (int *)malloc(VC.N_cylin_points*sizeof(int));
      old_bound_pts =    (int *)malloc(VC.N_cylin_points*sizeof(int));
      
      // Calculate the current Theta (i.e Theta_old)
      for(int b=0; b<VC.N_cylin_points; b++) // Looping over all boundary pts
      {
        const int p = VC.cylin_points[b];
        R = sqrt( points.x[p]*points.x[p] + points.y[p]*points.y[p] );     // Calc Radius
        Theta[b] = CalcEllipsoidAngle(points.x[p], points.y[p], R, R);   // Assume its a circle, and calculate Theta
        if (Theta[b] < 0.0) Theta[b] += 2.0*PI;
        Theta_old[b] = Theta[b];
        new_bound_pts[b] = -1; // Initialize
        old_bound_pts[b] = p;
      }
      
      // Now we sort
      quickSortDouble(Theta, 0, VC.N_cylin_points-1);  

      // Now sort the index according to the new Theta
      {
        for (int i=0; i<VC.N_cylin_points; i++)
        {
          int j;
          Found = false;
          for(j=0; j<VC.N_cylin_points; j++)
            if (Theta_old[i] == Theta[j]) {
              Found = true;
              break; }
            
          if (!(Found))              { printf("Error in ReadVolumeConservationMeshData. Theta was not found in Theta_ref");   exit(EXIT_FAILURE);}
          if (new_bound_pts[j] != -1){ printf("Error in ReadVolumeConservationMeshData. Theta was found twice");        exit(EXIT_FAILURE);}

          new_bound_pts[j] = old_bound_pts[i];
        }
      }

      // Check that all points are found and replaced
      {
        for(int b=0; b<VC.N_cylin_points; b++) // Looping over all pts
        {
          if (new_bound_pts[b] == -1) { printf("Error in ReadVolumeConservationMeshData. Not all points where found"); exit(EXIT_FAILURE);} 
          VC.cylin_points[b] = new_bound_pts[b];
        }
      }
      free(Theta);
      free(Theta_old);
      free(new_bound_pts);
      free(old_bound_pts);
    }
  }

  if (VC.Type == FORCE_VC_CONST_CELLS || VC.Type == FORCE_VC_CONST_CYLIN_CELLS)
  {
    sprintf(fileToRead, "./%s/vc-cell-triangles-index.dat", meshInputFolder);
    const int N_triangles = GetFileSize(fileToRead);
    if (N_triangles != triangles.N)  {
      printf("    ERROR!! The cell triangle index file is size %d but it should be %d\n", N_triangles, triangles.N);
      exit(EXIT_FAILURE); }

    VC.cell_traingle_index  = (int *)malloc(triangles.N*sizeof(int));
    ReadIntFile(fileToRead, triangles.N, VC.cell_traingle_index);

    // Check if the min cell index is 0
    int min_cell_index = 1000000;
    for(int t=0; t<triangles.N; t++)
      if (VC.cell_traingle_index[t]<min_cell_index)
        min_cell_index = VC.cell_traingle_index[t];

    if (min_cell_index != 0)  {
      printf("    ERROR!! Cell index should start with 0 not %d\n", min_cell_index);
      exit(EXIT_FAILURE); }
  }
}

//===============================================================================
// Initialization Functions


void VolumeConservationInitializeCells()
{
  int I_max = VC.cell_traingle_index[0], I_min = VC.cell_traingle_index[0];
  for(int t=0; t<triangles.N; t++)  {
    if (VC.cell_traingle_index[t] > I_max) I_max = VC.cell_traingle_index[t];
    if (VC.cell_traingle_index[t] < I_min) I_min = VC.cell_traingle_index[t];  }

  VC.N_cells = I_max - I_min + 1;   // Number of Cells
  VC.cell_V  = (double *)calloc(VC.N_cells, sizeof(double));
  VC.cell_V0 = (double *)calloc(VC.N_cells, sizeof(double));
  VC.cell_center_x = (double *)calloc(VC.N_cells, sizeof(double));
  VC.cell_center_y = (double *)calloc(VC.N_cells, sizeof(double));
  VC.cell_center_z = (double *)calloc(VC.N_cells, sizeof(double));
  VC.cell_num_tri  = (int *)calloc(VC.N_cells, sizeof(int));

  // Check that the orientation the triangle normal is outward
  //---------------------------------------------------------------
  // Find centroid of each cell
  CalculateCentroidOfCells();

  // Fix orientation of the triangular normal
  {
    double Mx, My, Mz, nVx, nVy, nVz;            
    for(int t=0; t<triangles.N; t++)
    {
      const int c = VC.cell_traingle_index[t]; // The cell index of the triangle 't'
      CalcTriangleCenter(t, &Mx, &My, &Mz);
      CalcTriangleNormal(t, &nVx, &nVy, &nVz);  
      if ((nVx*(Mx-VC.cell_center_x[c]) + nVy*(My-VC.cell_center_y[c]) + nVz*(Mz-VC.cell_center_z[c])) < 0)
      {
        const int swap_index = triangles.I1[t];
        triangles.I1[t] = triangles.I2[t];
        triangles.I2[t] = swap_index;
      }
    }
  }
  /* Calculates the initial volume of each cell from triangles. */
  #pragma omp parallel
  {
    int t, t_start, t_end, omp_rank = (int) omp_get_thread_num();

    double *cell_vo_loc = (double *)calloc(VC.N_cells, sizeof(double));

    GetLimitsOMPLoop(triangles.N, omp_nthreads, omp_rank, &t_start, &t_end);
    for(t=t_start; t<t_end; t++)
    {
      const int c = VC.cell_traingle_index[t]; // The cell index of the triangle 't'
      cell_vo_loc[c] += CalculatPartialVolumeFromTrianlge(t); // F.n.A 
    }

    #pragma omp critical
    {
      for(int c=0; c<VC.N_cells; c++)
        VC.cell_V0[c] += cell_vo_loc[c];
    }

    free(cell_vo_loc);
  }

  CalculateVolumeOfCells(); // Same as above but places the volume in VC.cell_V[c]
}


//===========================================================
// Force adding functions

void VolumeConservation3DCylinder()
{ 
  /* Function that conserves the volume of 3D cylinder by adding a volume conservation force */
  VC.volume = CalculateVolumeOfCylinder(); // Calculate the new total volume
  const double alpha = VC.alpha_cylin;
  const double V0_global = VC.volume_0;
  const double V_global = VC.volume;

  // Calculate and new volume force magnitude
  const double Fv = (-alpha*(V_global - V0_global)/V0_global)/3.0;
  #pragma omp parallel
  {
    double *fx_loc = (double *)calloc(forces.N,sizeof(double));
    double *fy_loc = (double *)calloc(forces.N,sizeof(double));
    double *fz_loc = (double *)calloc(forces.N,sizeof(double));

    double nVx, nVy, nVz, Tx, Ty, Tz, A;

    int t, t_start, t_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(triangles.N, omp_nthreads, omp_rank, &t_start, &t_end);

    for(t=t_start; t<t_end; t++)    // Looping over the triangle elements
    {
      const int p0 = triangles.I0[t];
      const int p1 = triangles.I1[t];
      const int p2 = triangles.I2[t];

      A = CalcTriangleArea2(t);       
      CalcTriangleCenter(t,  &Tx,  &Ty,  &Tz);       
      CalcTriangleNormal(t, &nVx, &nVy, &nVz);       

      fx_loc[p0] += Fv*A*nVx;   fy_loc[p0] += Fv*A*nVy;   fz_loc[p0] += Fv*A*nVz;       
      fx_loc[p1] += Fv*A*nVx;   fy_loc[p1] += Fv*A*nVy;   fz_loc[p1] += Fv*A*nVz;       
      fx_loc[p2] += Fv*A*nVx;   fy_loc[p2] += Fv*A*nVy;   fz_loc[p2] += Fv*A*nVz;       
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


void VolumeConservationCells()
{
  /* Function that conserves the volume of a group of cells by adding a volume conservation force */

  //Reset the values to '0'
  memset(VC.cell_V, 0, VC.N_cells* sizeof(double));  
  memset(VC.cell_center_x, 0, VC.N_cells* sizeof(double));  
  memset(VC.cell_center_y, 0, VC.N_cells* sizeof(double));  
  memset(VC.cell_center_z, 0, VC.N_cells* sizeof(double));  
  memset(VC.cell_num_tri,  0, VC.N_cells* sizeof(int));  

  // Calculate new volume of cells
  CalculateVolumeOfCells();

  // Calculate new centroid of each cell 
  // CalculateCentroidOfCells();  // NOT NEEDED NOW, SINCE TRIANGLE ORIENTATION WAS FIXED AT INIT TO CREATE AN OUTWARD NORMAL

  // Calculate the new volume force
  #pragma omp parallel
  {
    const double alpha = VC.alpha_cells;
    int t, t_start, t_end, omp_rank = (int) omp_get_thread_num();
    double Mx, My, Mz;            // Triangle Center
    double nVx, nVy, nVz, A;  // Vector

    double *fx_loc = (double *)calloc(forces.N,sizeof(double));
    double *fy_loc = (double *)calloc(forces.N,sizeof(double));
    double *fz_loc = (double *)calloc(forces.N,sizeof(double));

    GetLimitsOMPLoop(triangles.N, omp_nthreads, omp_rank, &t_start, &t_end);
    for(t=t_start; t<t_end; t++)
    {
      const int p0 = triangles.I0[t];
      const int p1 = triangles.I1[t];
      const int p2 = triangles.I2[t];
      const int c = VC.cell_traingle_index[t]; // The cell index of the triangle 'i'

      const double Fv = (-alpha*(VC.cell_V[c] - VC.cell_V0[c])/VC.cell_V0[c])/3.0;

      A = CalcTriangleArea2(t);
      CalcTriangleCenter(t, &Mx, &My, &Mz);
      CalcTriangleNormal(t, &nVx, &nVy, &nVz);       
      // Check if direction is correct // NOT NEEDED NOW, SINCE TRIANGLE ORIENTATION WAS FIXED AT INIT TO CREATE AN OUTWARD NORMAL
      //if ((nVx*(Mx-VC.cell_center_x[c]) + nVy*(My-VC.cell_center_y[c]) + nVz*(Mz-VC.cell_center_z[c])) < 0)
      //  nVx = -1*nVx, nVy = -1*nVy, nVz = -1*nVz;

      fx_loc[p0] += Fv*A*nVx;   fy_loc[p0] += Fv*A*nVy;   fz_loc[p0] += Fv*A*nVz;       
      fx_loc[p1] += Fv*A*nVx;   fy_loc[p1] += Fv*A*nVy;   fz_loc[p1] += Fv*A*nVz;       
      fx_loc[p2] += Fv*A*nVx;   fy_loc[p2] += Fv*A*nVy;   fz_loc[p2] += Fv*A*nVz; 
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

  // Calculate the new total volume
  VC.volume_cells = CalculateTotalVolumeOfCells();
}

void VolumeConservation3DCylinderAndCells(){
  VolumeConservationCells();
  VolumeConservation3DCylinder();
}

// Final
void ChooseVolumeConservationFunc()
{
  printf("Choosing the volume conservartion force function\n"); 
  if (VC.Type == FORCE_VC_NULL){
    printf("  No volume conservartion force function is choosen\n"); 
    VC.Func = CalculateNullForce;}
  else if (VC.Type == FORCE_VC_CONST_CYLIN){
    printf("  Chosen volume conservartion force function is for cylinder shape domain with constant gammas\n"); 
    VC.volume_0 = CalculateVolumeOfCylinder();              // Initial Volume
    VC.volume = VC.volume_0;
    VC.Func = VolumeConservation3DCylinder;}
  else if (VC.Type == FORCE_VC_CONST_CELLS){
    printf("  Chosen volume conservartion force function is for each cell shape domain with constant gammas\n"); 
    VolumeConservationInitializeCells();
    VC.volume_0 = CalculateTotalVolumeOfCells();            // Initial Total volume
    VC.volume = VC.volume_0;
    VC.Func = VolumeConservationCells;}
  else if (VC.Type == FORCE_VC_CONST_CYLIN_CELLS){
    printf("  Chosen volume conservartion force function is for cylinder shape and each cell shape domain with constant gammas\n"); 
    VolumeConservationInitializeCells();
    VC.volume_0 = CalculateVolumeOfCylinder();              // Initial Volume
    VC.volume = VC.volume_0;
    VC.volume_cells = CalculateTotalVolumeOfCells();
    VC.Func = VolumeConservation3DCylinderAndCells;}
  else{
    printf("  ERROR!! Volume conservartion function is undefined\n");
    exit(EXIT_FAILURE);}
}

void ForcesAddVolumeConservationForce(){ VC.Func();}

// Extra methods
//-----------------------------------------------------------------------
inline double CalculatPartialVolumeFromTrianlge(int t)
{
  /* 
  The volume is calculated based on the divergecne theorem shown below:

    Volume = Int_V dV 
           = Int_V Div.vec(F) dV    where Div.vec(F) = 1
           = Int_S vec(F).hat(n) dA
           = Sum_i vec(F).hat(n) A_i

    for our case we will choose vec(F) = x/3*i + y/3*j + z/3*k so that Div.vec(F) = 1
    where vec(F(x,y,z)) will be calculated at the center of the triangle
    and hat(n) will be the normalized vector of the triangle 

    This function ONLY calculates the partial volume from ONE triangle.
    To get the full volume, loop over all triangles and sum the partial volume
  */
  double Tx, Ty, Tz;  // Local center of triangle
  double A, nVx, nVy, nVz;
  A = CalcTriangleArea2(t);
  CalcTriangleCenter(t,  &Tx,  &Ty,  &Tz);
  CalcTriangleNormal(t, &nVx, &nVy, &nVz);

  return (Tx*nVx + Ty*nVy + Tz*nVz)*A/3.; // F.n.A 
}

double CalculateAreaFromLines()
{
  double A_global = 0.0;
  #pragma omp parallel
  {
    //const double Cx = VC.center[0], Cy=VC.center[1], Cz=VC.center[2]; // Global center of domain
    double A_local = 0; // Local area

    int l, l_start, l_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(lines.N, omp_nthreads, omp_rank, &l_start, &l_end);
    for(l=l_start; l<l_end; l++)
    {
      const int p0 = lines.I0[l];
      const int p1 = lines.I1[l];
      A_local += points.x[p0]*points.y[p1]-points.y[p0]*points.x[p1];
    }

    #pragma omp critical
    {
      A_global += A_local;
    }
  }
  A_global = A_global/2.0;
  return A_global;
}

double CalculateAreaFromPointsOnCylinder()
{
  double A_cylin = 0.0;
  int p0, p1;

  for(int j=0; j<VC.N_cylin_points; j++)
  {
    p0 = VC.cylin_points[j];
    if (j<VC.N_cylin_points-1) p1 = VC.cylin_points[j+1];
    else                       p1 = VC.cylin_points[0];

    A_cylin += points.x[p0]*points.y[p1]-points.y[p0]*points.x[p1];
  }

  A_cylin = A_cylin/2.0;
  return A_cylin;
}


double CalculateVolumeOfCylinder()
{
  /* The function works for a cylinder with triangles defining the outer parameter and lines defining the circles.  */
  double volume = 0.0;
  #pragma omp parallel
  {    
    /*   Uses the CalculatPartialVolumeFromTrianlge function to calculate the volume that the top surface triangles contribute to the volume. */
    double volume_local = 0;       // Local volume

    int ct, ct_start, ct_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(VC.N_cylin_triangles, omp_nthreads, omp_rank, &ct_start, &ct_end);
    for(ct=ct_start; ct<ct_end; ct++)
    {
      const int t = VC.cylin_triangles[ct];  
      volume_local += CalculatPartialVolumeFromTrianlge(t); // F.n.A 
    }

    #pragma omp critical
    {
      volume += volume_local;
    }
  }

  double A_cylinder = CalculateAreaFromPointsOnCylinder();       // Calculate the area of the two circles from points on the surface of cylinder
  volume += (VC.cylin_circle1_point[0]*VC.cylin_circle1_normal[0] + VC.cylin_circle1_point[1]*VC.cylin_circle1_normal[1] + VC.cylin_circle1_point[2]*VC.cylin_circle1_normal[2])*A_cylinder/3.0;
  volume += (VC.cylin_circle2_point[0]*VC.cylin_circle2_normal[0] + VC.cylin_circle2_point[1]*VC.cylin_circle2_normal[1] + VC.cylin_circle2_point[2]*VC.cylin_circle2_normal[2])*A_cylinder/3.0;
  return volume;
}

void CalculateVolumeOfCells()
{
  /* Calculates the volume of each cell from triangles. */
  #pragma omp parallel
  {
    int t, t_start, t_end, omp_rank = (int) omp_get_thread_num();

    double *cell_vo_loc = (double *)calloc(VC.N_cells, sizeof(double));

    GetLimitsOMPLoop(triangles.N, omp_nthreads, omp_rank, &t_start, &t_end);
    for(t=t_start; t<t_end; t++)
    {
      const int c = VC.cell_traingle_index[t]; // The cell index of the triangle 't'
      cell_vo_loc[c] += CalculatPartialVolumeFromTrianlge(t); // F.n.A 
    }

    #pragma omp critical
    {
      for(int c=0; c<VC.N_cells; c++)
        VC.cell_V[c] += cell_vo_loc[c];
    }

    free(cell_vo_loc);
  }
}

double CalculateTotalVolumeOfCells()
{
  double volume = 0.0;
  for(int c=0; c<VC.N_cells; c++) 
    volume += VC.cell_V[c];
  return volume;
}

void CalculateCentroidOfCells()
{
  memset(VC.cell_center_x, 0, VC.N_cells* sizeof(double));  
  memset(VC.cell_center_y, 0, VC.N_cells* sizeof(double));  
  memset(VC.cell_center_z, 0, VC.N_cells* sizeof(double));  
  memset(VC.cell_num_tri,  0, VC.N_cells* sizeof(int));  

  #pragma omp parallel
  {
    double Mx, My, Mz;            // Triangle Center
    int t, t_start, t_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(triangles.N, omp_nthreads, omp_rank, &t_start, &t_end);

    double *cell_center_x_loc = (double *)calloc(VC.N_cells, sizeof(double));
    double *cell_center_y_loc = (double *)calloc(VC.N_cells, sizeof(double));
    double *cell_center_z_loc = (double *)calloc(VC.N_cells, sizeof(double));
    int    *cell_num_tri_loc  = (int *)calloc(VC.N_cells, sizeof(int));

    for(t=t_start; t<t_end; t++)
    {
      const int c = VC.cell_traingle_index[t]; // The cell index of the triangle 't'
      CalcTriangleCenter(t, &Mx, &My, &Mz);
      cell_center_x_loc[c] += Mx;
      cell_center_y_loc[c] += My;
      cell_center_z_loc[c] += Mz;
      cell_num_tri_loc[c]  += 1;
    }

    #pragma omp critical
    {
      for(int c=0; c<VC.N_cells; c++)
      {
        VC.cell_center_x[c] += cell_center_x_loc[c];
        VC.cell_center_y[c] += cell_center_y_loc[c];
        VC.cell_center_z[c] += cell_center_z_loc[c];
        VC.cell_num_tri[c]  += cell_num_tri_loc[c];
      }
    }
    free(cell_center_x_loc);
    free(cell_center_y_loc);
    free(cell_center_z_loc);
    free(cell_num_tri_loc);
  }

  for(int c=0; c<VC.N_cells; c++)
  {
    VC.cell_center_x[c] = VC.cell_center_x[c]/( (double) VC.cell_num_tri[c]);
    VC.cell_center_y[c] = VC.cell_center_y[c]/( (double) VC.cell_num_tri[c]);
    VC.cell_center_z[c] = VC.cell_center_z[c]/( (double) VC.cell_num_tri[c]);
  }
}


