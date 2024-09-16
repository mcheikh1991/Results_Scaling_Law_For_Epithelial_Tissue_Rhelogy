
void DefaultDynamicMeshParameters()
{
  DM.Type = DYNAMIC_MESH_NULL;   // Dynamics Mesh 
  DM.L_max = 1e8;
  DM.L_min = 0.0;
}

void PrintDynamicMeshProperties()
{
  printf("Dynamic Mesh Properties:\n");
  printf("-----------------------------------\n");
  printf("  DM.Type:\t\t%s\n", DYNAMIC_MESH_CHAR[DM.Type]);
  if (DM.Type == DYNAMIC_MESH_ELASTIC_LENGTH)
  {
    printf("  DM.L_max:\t\t%lf\n", DM.L_max);
    printf("  DM.L_min:\t\t%lf\n", DM.L_min);
    if (LE.Type != FORCE_LE_INITIAL_DISP_MULTI){
      printf("    ERROR!! '%s' is only compatible with '%s'\n", DYNAMIC_MESH_CHAR[DM.Type], FORCE_LE_CHAR[FORCE_LE_INITIAL_DISP_MULTI]);
      exit(EXIT_FAILURE);}
    if (BF.Bool_Force_triangle || BF.Bool_Force_square){
      printf("    ERROR!! '%s' is not compatible with squares and triangles yet \n", DYNAMIC_MESH_CHAR[DM.Type]);
      exit(EXIT_FAILURE);}
    if (ST.Type != FORCE_ST_NULL && ST.Type != FORCE_ST_MULTI){
      printf("    ERROR!! '%s' is only compatible with '%s'\n", DYNAMIC_MESH_CHAR[DM.Type], FORCE_ST_CHAR[FORCE_ST_MULTI]);
      exit(EXIT_FAILURE);}
  }
  printf("\n");
}

bool ReadDynamicMeshParameters(const char* str_argv1, const char* str_argv2)
{
  if (strcmp(str_argv1,"-dm_type") == 0){    
    DM.Type = (DYNAMIC_MESH_TYPE) atol(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-dm_l_max") == 0){    
    DM.L_max = atof(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-dm_l_min") == 0){    
    DM.L_min = atof(str_argv2); 
    return true;}
    
  return false;
}


//===============================================================================
int DMAddPoint(const double x, const double y, const double z);
void DMRemovePoint(const int i);
void DMAddLine(const int p0, const int p1, const double L0, const double LE_M);
void DMAddLineWithST(const int p0, const int p1, const double L0, const double LE_M, const double ST_M);
void DMRemoveLine(const int i);

void DyamicMeshDivideLine(const int l)
{
  /* This function will divide line l in half l = l0 + l1
     It will create a new point p and saves it in points and forces
     and replace l with l0 and create a new space for l1.
     Thus:
      - new point:  p
      - new line l0: p0 - p  (replaces old line l)
      - new line l1: p - p1
  */
  //const int p0 = lines.I0[l];
  const int p1    = lines.I1[l];      // p1 on line l
  const double L0 = lines.len[l];     // Initial length
  const double LE_M = LE.M_line[l];     // Spring constant

  // Create new point p and save it
  double px, py, pz;
  CalcMidpoint(l, &px, &py, &pz);
  int p = DMAddPoint(px, py, pz);     // p is the index of the new point

  // Adjust the old line l to become l0
  lines.I1[l] = p;
  lines.len[l] = L0/2.0; // L0
  LE.M_line[l] = 2*LE_M;

  // Add new line, with new L0, and new K=2*K_old
  DMAddLine(p, p1, L0/2.0, 2*LE_M); // index of new line l1 with new L0 and new k

}

void DyamicMeshDivideLineWithST(const int l)
{
  /* This function is similar to 'DyamicMeshDivideLine' but also works on 
      Surface tension multiplier
  */
  //const int p0 = lines.I0[l];
  const int p1    = lines.I1[l];      // p1 on line l
  const double L0 = lines.len[l];     // Initial length
  const double LE_M  = LE.M_line[l];  // Spring constant multiplier
  const double ST_M  = ST.M_line[l];  // Spring constant multiplier

  // Create new point p and save it
  double px, py, pz;
  CalcMidpoint(l, &px, &py, &pz);
  int p = DMAddPoint(px, py, pz);     // p is the index of the new point

  // Adjust the old line l to become l0
  lines.I1[l] = p;
  lines.len[l] = L0/2.0; // L0
  LE.M_line[l] = 2*LE_M;

  // Add new line, with new L0, and new K=2*K_old
  DMAddLineWithST(p, p1, L0/2.0, 2*LE_M, ST_M); // index of new line l1 with new L0 and new k

}


void DyamicMeshMergePointsAndRemoveLine(const int l)
{
  /* This function will merge the points of line l together and remove line l
  */
  const int p0 = lines.I0[l];
  const int p1 = lines.I1[l];
  int p_remove, p_keep;

  // Remove the higher index p
  if (p1 > p0)    
    p_remove = p1, p_keep = p0;
  else            
    p_remove = p0, p_keep = p1;
  DMRemovePoint(p_remove);

  // Remove the line l
  DMRemoveLine(l);

  // Loop over the indices and replace p by p_keep
  for (int i=0; i<lines.N; i++)
  {
    if(lines.I0[i] == p_remove)
      lines.I0[i] = p_keep;
    else if (lines.I1[i] == p_remove)
      lines.I1[i] = p_keep;
  }
}



void DynamicMeshElasticLength()
{
  for (int l=0; l<lines.N;l++)
  {
    const double L = CalcLineLength(l); // Current length
    if(L > DM.L_max)  
    {
      DyamicMeshDivideLine(l);
      //Create Constrain again
      CN.Func_free();
      CN.Func_create();
    }
  }
}

void DynamicMeshElasticLengthAndSurfaceTension()
{
  for (int l=0; l<lines.N;l++)
  {
    const double L = CalcLineLength(l); // Current length
    if(L > DM.L_max)  
    {
      DyamicMeshDivideLineWithST(l);
      //Create Constrain again
      CN.Func_free();
      CN.Func_create();
    }
  }
}

void DynamicMeshNull(){}

void ChooseDynamicMeshFunc()
{
  printf("Choosing the dynamic mesh function\n"); 
  if (DM.Type == DYNAMIC_MESH_NULL){
    printf("  No dynamic mesh function is choosen\n"); 
    DM.Func = DynamicMeshNull;}
  else if (DM.Type == DYNAMIC_MESH_ELASTIC_LENGTH) 
  {
    if (LE.Type != FORCE_LE_INITIAL_DISP_MULTI){
      printf("  ERROR!! Dynamic mesh function needs to have the linear-elasticity option be 'FORCE_LE_INITIAL_DISP_MULTI'\n");
      exit(EXIT_FAILURE);}
    else
    {
      if (ST.Type == FORCE_ST_NULL){
        printf("  Dynamic mesh based on the length is chosen with L_max = %lf and L_min = %lf\n", DM.L_max, DM.L_min); 
        DM.Func = DynamicMeshElasticLength;}
      else if (ST.Type == FORCE_ST_MULTI){
        printf("  Dynamic mesh based on the length with ST option is chosen with L_max = %lf and L_min = %lf\n", DM.L_max, DM.L_min); 
        DM.Func = DynamicMeshElasticLengthAndSurfaceTension;}
    }
  }
  else{
    printf("  ERROR!! Dynamic mesh function is undefined\n");
    exit(EXIT_FAILURE);}
}

void DynamicMesh(){ DM.Func();}


//============================================================
// Private Functions
//============================================================
// Do not use them since they will not resolve all issues


int DMAddPoint(const double x, const double y, const double z)
{
  /* Adds a point at the end of the point arrays (x,y,z)
     and the force arrays.
     Returns the index of the new point */
  const int i = points.N;
  points.N = points.N + 1;
  points.x = realloc(points.x, points.N*sizeof(double));
  points.y = realloc(points.y, points.N*sizeof(double));
  points.z = realloc(points.z, points.N*sizeof(double));

  points.x[i] = x;
  points.y[i] = y;
  points.z[i] = z;

  // Adjusting the forces
  forces.N = points.N;
  forces.x = realloc(forces.x, points.N*sizeof(double));
  forces.y = realloc(forces.y, points.N*sizeof(double));
  forces.z = realloc(forces.z, points.N*sizeof(double));
  forces.x[i] = 0.0;
  forces.y[i] = 0.0;
  forces.z[i] = 0.0;

  return i; // returns the index where the point was added
}

void DMRemovePoint(const int i)
{
  /* Removes a point at the index i of the point arrays (x,y,z)
     and the force arrays.
     NOTE: For this method to be correct one has to also 
           fix the lines, triangles, and sequare arrays if they exist*/

  int j = i;
  points.N = points.N - 1;
  while(j<points.N) { 
    points.x[j] = points.x[j+1];   
    points.y[j] = points.y[j+1];   
    points.z[j] = points.z[j+1];   
    j++; 
  }  
  points.x = realloc(points.x, points.N*sizeof(double));
  points.y = realloc(points.y, points.N*sizeof(double));
  points.z = realloc(points.z, points.N*sizeof(double));

  j = i; // Restart
  forces.N = points.N;
  while(j<forces.N) { 
    forces.x[j] = forces.x[j+1];   
    forces.y[j] = forces.y[j+1];   
    forces.z[j] = forces.z[j+1];   
    j++; 
  }  
  forces.x = realloc(forces.x, forces.N*sizeof(double));
  forces.y = realloc(forces.y, forces.N*sizeof(double));
  forces.z = realloc(forces.z, forces.N*sizeof(double));
}

void DMAddLine(const int p0, const int p1, const double L0, const double LE_M)
{
  const int i = lines.N;
  lines.N = lines.N + 1;
  lines.I0 = realloc(lines.I0, lines.N*sizeof(int));
  lines.I1 = realloc(lines.I1, lines.N*sizeof(int));
  lines.len = realloc(lines.len, lines.N*sizeof(double));
  LE.M_line = realloc(LE.M_line, lines.N*sizeof(double));

  lines.I0[i]  = p0;
  lines.I1[i]  = p1;
  lines.len[i] = L0; 
  LE.M_line[i] = LE_M;
}

void DMAddLineWithST(const int p0, const int p1, const double L0, const double LE_M, const double ST_M)
{
  const int i = lines.N;
  lines.N   = lines.N + 1;
  lines.I0  = realloc(lines.I0,   lines.N*sizeof(int));
  lines.I1  = realloc(lines.I1,   lines.N*sizeof(int));
  lines.len = realloc(lines.len,  lines.N*sizeof(double));
  LE.M_line = realloc(LE.M_line,  lines.N*sizeof(double));
  ST.M_line = realloc(ST.M_line,  lines.N*sizeof(double));

  lines.I0[i]  = p0;
  lines.I1[i]  = p1;
  lines.len[i] = L0; 
  LE.M_line[i] = LE_M;
  ST.M_line[i] = ST_M;
}


void DMRemoveLine(const int i)
{
  int j = i;
  lines.N = lines.N - 1;
  while(j<lines.N) {   
    lines.I0[j] = lines.I0[j+1];
    lines.I1[j] = lines.I1[j+1];
    lines.len[j] = lines.len[j+1];
    LE.M_line[j] = LE.M_line[j+1];
    j++; 
  }  
  lines.I0 = realloc(lines.I0, lines.N*sizeof(int));
  lines.I1 = realloc(lines.I1, lines.N*sizeof(int));
  lines.len = realloc(lines.len, lines.N*sizeof(double));
  LE.M_line = realloc(LE.M_line, lines.N*sizeof(double));
}