
void DefaultBoxParameters()
{
  // Box Properties
  boxes.Max_Points = 100;
  boxes.Max_Lines = 100;
  boxes.Max_Triangles = 100;
  boxes.Max_Squares = 100;
  boxes.M = 1.0;
}


void PrintBoxProperties()
{
  printf("Box Parameters:\n");
  printf("-----------------------------------\n");
  printf("  boxes.Max_Points:\t%d\n", boxes.Max_Points);
  printf("  boxes.Max_Lines:\t%d\n", boxes.Max_Lines);
  printf("  boxes.Max_Triangles:\t%d\n", boxes.Max_Triangles);
  printf("  boxes.Max_Squares:\t%d\n", boxes.Max_Squares);
  printf("  boxes.M:\t\t%lf\n", boxes.M);
  printf("\n");
}

bool ReadBoxParameters(const char* str_argv1, const char* str_argv2)
{
  if (strcmp(str_argv1,"-boxes_max_points") == 0){    
    boxes.Max_Points = atol(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-boxes_max_lines") == 0){    
    boxes.Max_Lines = atol(str_argv2);
    return true;}

  if (strcmp(str_argv1,"-boxes_max_triangles") == 0){    
    boxes.Max_Triangles = atol(str_argv2);
    return true;}

  if (strcmp(str_argv1,"-boxes_max_squares") == 0){    
    boxes.Max_Squares = atol(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-boxes_m") == 0){    
    boxes.M = atof(str_argv2);
    return true;}
  
  return false;
}

//=====================
void BoxesInitialize()
{
  /* Initialized the variables in boxes */
  boxes.N = -1;

  boxes.Lx = -1;
  boxes.Ly = -1;
  boxes.Lz = -1;

  boxes.Nx = -1;
  boxes.Ny = -1;
  boxes.Nz = -1;

  boxes.dx = -1;
  boxes.dy = -1;
  boxes.dz = -1;
 
  boxes.Xmin =  1e10;  
  boxes.Xmax = -1e10;  
  boxes.Ymin =  1e10;  
  boxes.Ymax = -1e10;  
  boxes.Zmin =  1e10;  
  boxes.Zmax = -1e10;  
}



inline int BoxesCalcIndex(const double Px, const double Py, const double Pz, 
  const int Ax, const int Ay, const int Az) 
{
  /* Calculate the index of the box that contains the point (Px, Py, Pz) 
    Ax, Ay, Az if we wanted to get the box to the left(Ax=-1), right(Ax=+1) .. for the main box*/
  int ix = (int) ((Px - boxes.Xmin)/boxes.dx) + Ax; // The index in x direction
  int iy = (int) ((Py - boxes.Ymin)/boxes.dy) + Ay; // The index in y direction
  int iz = (int) ((Pz - boxes.Zmin)/boxes.dz) + Az; // The index in z direction
  ix = MAX(MIN(ix,boxes.Nx-1),0);
  iy = MAX(MIN(iy,boxes.Ny-1),0);
  iz = MAX(MIN(iz,boxes.Nz-1),0);
  return ix + iy*boxes.Nx + iz*boxes.Nx*boxes.Ny*(DIM-2); // The Dim - 2 will make iz = 0 if Dim = 2 
}

bool BoxesCheckIfPointInDomain(const double Px, const double Py, const double Pz) 
{
  /* Checks if the points coordinates (Px, Py, Pz) are in the domain of boxes */
  bool NodeInDomain = true; 
  if      (  Px < boxes.Xmin || Px > boxes.Xmax ) 
    NodeInDomain = false;
  else if (  Py < boxes.Ymin || Py > boxes.Ymax ) 
    NodeInDomain = false;
  else if ( (Pz < boxes.Zmin || Pz > boxes.Zmax ) && (DIM == 3) )
    NodeInDomain = false;
  return NodeInDomain;
}


bool BoxesCheckIfLineInBox(const int B, const int l)
{
  /* Checks if the box already has the line element 'l'. Returns True if it has it False otherwise
       B is the box number, l is the line index */
  bool LineFound = false; 
  for (int k = 0; k < boxes.box_lines[B][0]; k++) // looping over the elements in box
    if ( boxes.box_lines[B][k+1] == l)
    { LineFound = true;
      break;}
  return LineFound;
}

bool BoxesCheckIfPointInBox(const int B, const int p)
{
  /* Checks if the box already has the point 'p'. Returns True if it has it False otherwise
       B is the box number, p is the point index */
  bool PointFound = false; 
  for (int k = 0; k < boxes.box_points[B][0]; k++) // looping over the elements in box
    if ( boxes.box_points[B][k+1] == p)
    { PointFound = true;
      break;}
  return PointFound;
}



int BoxesFindDomainParameters()
{
  /* Find the size of the domain and number of boxes it has */
  printf("  Domain will be divided in to boxes\n");
  // Part 1: Finding the local size
  printf("    Finding the local size of the box domain\n");
  for (int i = 0; i < points.N; i++)
  {
    const double x = points.x[i];
    const double y = points.y[i];
    const double z = points.z[i];

    if ( x > boxes.Xmax) boxes.Xmax = x;
    if ( y > boxes.Ymax) boxes.Ymax = y;
    if ( z > boxes.Zmax) boxes.Zmax = z;

    if ( x < boxes.Xmin) boxes.Xmin = x;
    if ( y < boxes.Ymin) boxes.Ymin = y;
    if ( z < boxes.Zmin) boxes.Zmin = z;
  }

  boxes.Lx = (boxes.Xmax - boxes.Xmin); 
  boxes.Ly = (boxes.Ymax - boxes.Ymin); 
  boxes.Lz = (boxes.Zmax - boxes.Zmin); 

  // Part 2: Extend the length by 5% total (2.5% each direction)
  boxes.Xmin = boxes.Xmin - 0.025*boxes.Lx;
  boxes.Ymin = boxes.Ymin - 0.025*boxes.Ly;
  boxes.Zmin = boxes.Zmin - 0.025*boxes.Lz;
  boxes.Xmax = boxes.Xmax + 0.025*boxes.Lx;
  boxes.Ymax = boxes.Ymax + 0.025*boxes.Ly;
  boxes.Zmax = boxes.Zmax + 0.025*boxes.Lz;

  // Update size now
  boxes.Lx = (boxes.Xmax - boxes.Xmin); 
  boxes.Ly = (boxes.Ymax - boxes.Ymin); 
  boxes.Lz = (boxes.Zmax - boxes.Zmin);

  // Part 3: Finding the size of each box and how many numbers
  double D_avg = 0;
  for (int l=0; l<lines.N; l++)
    D_avg += CalcLineLength(l);

  D_avg = (D_avg/lines.N)*boxes.M; // M times larger than the average length

  // Part 4: Calculate the number of boxes
  // Remove 1 so that we can make the boxe size > element size
  boxes.Nx = (int) fmax( (int) (boxes.Lx/D_avg - 1), 1);
  boxes.Ny = (int) fmax( (int) (boxes.Ly/D_avg - 1), 1);
  boxes.Nz = (int) fmax( (int) (boxes.Lz/D_avg - 1), 1);
  if (DIM == 2) boxes.Nz = 1;

  boxes.N = boxes.Nx*boxes.Ny*boxes.Nz;

  boxes.dx = boxes.Lx / boxes.Nx;
  boxes.dy = boxes.Ly / boxes.Ny;
  boxes.dz = boxes.Lz / boxes.Nz;
  printf("    The boxes created have a size of (%lf, %lf, %lf)\n", boxes.dx, boxes.dy, boxes.dz );
  printf("    The number of boxes is N:%d with Nx:%d, Ny:%d, and Nz:%d \n", boxes.N, boxes.Nx, boxes.Ny, boxes.Nz);

  return 0;
}


int BoxesMallocArrays()
{
  int i,j;
  if (BF.Bool_Force_point)
  {
    printf("  Adding points to boxes\n");
    boxes.box_points = (int **)malloc(boxes.N*sizeof(int *));
    for(i = 0; i < boxes.N; ++i)
    { 
      boxes.box_points[i]  = (int *)malloc(sizeof(int)*(boxes.Max_Points+1));
      boxes.box_points[i][0] = 0;
      for(j = 0; j < boxes.Max_Points; ++j) 
        boxes.box_points[i][j+1] = -1;
    }
  }

  if (BF.Bool_Force_line)
  {
    printf("  Adding lines to boxes\n");
    boxes.box_lines  = (int **)malloc(boxes.N*sizeof(int *));
    for(i = 0; i < boxes.N; ++i)
    { 
      boxes.box_lines[i]   = (int *)malloc(sizeof(int)*(boxes.Max_Lines+1));
      boxes.box_lines[i][0] = 0;
      for(j = 0; j < boxes.Max_Lines; ++j) 
        boxes.box_lines[i][j+1] = -1;
    }
  }

  if (BF.Bool_Force_triangle)
  {
    printf("  Adding triangles to boxes\n");
    boxes.box_triangles = (int **)malloc(boxes.N*sizeof(int *));
    for(i = 0; i < boxes.N; ++i)
    { 
      boxes.box_triangles[i]  = (int *)malloc(sizeof(int)*(boxes.Max_Triangles+1));
      boxes.box_triangles[i][0] = 0;
      for(j = 0; j < boxes.Max_Triangles; ++j) 
        boxes.box_triangles[i][j+1] = -1;
    }
  }

  if (BF.Bool_Force_square)
  {
    printf("  Adding squares to boxes\n");
    boxes.box_squares = (int **)malloc(boxes.N*sizeof(int *));
    for(i = 0; i < boxes.N; ++i)
    { 
      boxes.box_squares[i]  = (int *)malloc(sizeof(int)*(boxes.Max_Squares+1));
      boxes.box_squares[i][0] = 0;
      for(j = 0; j < boxes.Max_Squares; ++j) 
        boxes.box_squares[i][j+1] = -1;
    }
  }

  return 0;
}

//=====================================================================================
void BoxesResetPoints()
{
  #pragma omp parallel
  {
    int i, i_start, i_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(boxes.N, omp_nthreads, omp_rank, &i_start, &i_end);
    for(i = i_start; i < i_end; ++i)
      boxes.box_points[i][0] = 0;
  }
}

void BoxesResetLines()
{  
  #pragma omp parallel
  {
    int i, i_start, i_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(boxes.N, omp_nthreads, omp_rank, &i_start, &i_end);
    for(i = i_start; i < i_end; ++i)
      boxes.box_lines[i][0] = 0;
  }
}

void BoxesResetTriangles()
{
  #pragma omp parallel
  {
    int i, i_start, i_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(boxes.N, omp_nthreads, omp_rank, &i_start, &i_end);
    for(i = i_start; i < i_end; ++i)
      boxes.box_triangles[i][0] = 0;
  }
}


void BoxesResetSquares()
{
  #pragma omp parallel
  {
    int i, i_start, i_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(boxes.N, omp_nthreads, omp_rank, &i_start, &i_end);
    for(i = i_start; i < i_end; ++i)
      boxes.box_squares[i][0] = 0;
  }
}


void BoxesPlacePointsParallel()
{
  BoxesResetPoints();
  #pragma omp parallel
  {
    int B;
    int p, p_start, p_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(points.N, omp_nthreads, omp_rank, &p_start, &p_end);

    int **box_points_local = (int **)malloc(boxes.N*sizeof(int *));
    for(int i = 0; i < boxes.N; ++i) { 
      box_points_local[i]  = (int *)calloc(boxes.Max_Points+1,sizeof(int));}

    for (p = p_start; p < p_end; p++) // Looping over local internal elements
    {
      B = BoxesCalcIndex(points.x[p], points.y[p], points.z[p], 0, 0, 0);
      box_points_local[B][0]++; 
      box_points_local[B][boxes.box_points[B][0]] = p;
    }

    #pragma omp critical
    {
      for(int b=0; b<boxes.N; b++)
        {
          const int b_start = boxes.box_points[b][0];
          for (int j=0; j<box_points_local[b][0]; j++){
            boxes.box_points[b][b_start+j+1] = box_points_local[b][j+1]; }
          boxes.box_points[b][0] += box_points_local[b][0]; 
        }
    }

    for(int i = 0; i < boxes.N; i++){ free(box_points_local[i]);}
      free(box_points_local);
  }
}


void BoxesPlacePoints()
{
  int B;
  BoxesResetPoints();
  for (int p = 0; p < points.N; p++) // Looping over local internal elements
  {
    B = BoxesCalcIndex(points.x[p], points.y[p], points.z[p], 0, 0, 0);
    boxes.box_points[B][0]++; 
    boxes.box_points[B][boxes.box_points[B][0]] = p;
  }
}

void BoxesPlaceLinesPrallel()
{
  BoxesResetLines();
  #pragma omp parallel
  {
    double Mx, My, Mz;
    int B;
    int l, l_start, l_end, omp_rank = (int) omp_get_thread_num();
    GetLimitsOMPLoop(lines.N, omp_nthreads, omp_rank, &l_start, &l_end);

    int **box_lines_local = (int **)malloc(boxes.N*sizeof(int *));
    for(int i = 0; i < boxes.N; ++i) { 
      box_lines_local[i]  = (int *)calloc(boxes.Max_Lines+1,sizeof(int));}

    for (l = l_start; l < l_end; l++) // Looping over local internal elements
    {
      CalcMidpoint(l, &Mx, &My, &Mz);
      B = BoxesCalcIndex(Mx, My, Mz, 0, 0, 0);
      box_lines_local[B][0]++; 
      box_lines_local[B][boxes.box_lines[B][0]] = l;
    }

    #pragma omp critical
    {
      for(int b=0; b<boxes.N; b++)
        {
          const int b_start = boxes.box_lines[b][0];
          for (int j=0; j<box_lines_local[b][0]; j++){
            boxes.box_lines[b][b_start+j+1] = box_lines_local[b][j+1]; }
          boxes.box_lines[b][0] += box_lines_local[b][0]; 
        }
    }

    for(int i = 0; i < boxes.N; i++){ free(box_lines_local[i]);}
      free(box_lines_local);
  }
}


void BoxesPlaceLines()
{
  int B;
  double Mx, My, Mz;
  BoxesResetLines();
  for (int l = 0; l < lines.N; l++) // Looping over local internal elements
  {
    CalcMidpoint(l, &Mx, &My, &Mz);
    B = BoxesCalcIndex(Mx, My, Mz, 0, 0, 0);
    boxes.box_lines[B][0]++; 
    boxes.box_lines[B][boxes.box_lines[B][0]] = l;
    //printf("Box %d: %d, (%lf, %lf, %lf)\n", B, boxes.box_lines[B][0], Mx, My, Mz);
  }
}

void BoxesPlaceTriangles()
{
  int B;
  double Mx, My, Mz;
  BoxesResetTriangles();
  for (int t = 0; t < triangles.N; t++) // Looping over local internal elements
  {
    CalcTriangleCenter(t, &Mx, &My, &Mz);
    B = BoxesCalcIndex(Mx, My, Mz, 0, 0, 0);
    boxes.box_triangles[B][0]++; 
    boxes.box_triangles[B][boxes.box_triangles[B][0]] = t;
    //printf("T in box %d is %d\n",B, boxes.box_triangles[B][0]);
  }
}

void BoxesPlaceSquares()
{
  int B;
  double Mx, My, Mz;
  BoxesResetSquares();
  for (int s = 0; s < squares.N; s++) // Looping over local internal elements
  {
    CalcSquareCenter(s, &Mx, &My, &Mz);
    B = BoxesCalcIndex(Mx, My, Mz, 0, 0, 0);
    boxes.box_squares[B][0]++; 
    boxes.box_squares[B][boxes.box_squares[B][0]] = s;
    //printf("S in box %d is %d\n",B, boxes.box_squares[B][0]);
  }
}

void BoxesNull(){}

void ChooseBoxFunctions()
{  

  printf("Choosing the box function/s\n"); 
  BoxesInitialize();
  if ( (!(BF.Bool_Force_point) && !(BF.Bool_Force_line) && !(BF.Bool_Force_triangle) && !(BF.Bool_Force_square)) || 
       BF.Type == FORCE_BF_SET_POINTS)
  {
    printf("  Domain will not be divided into boxes\n");
    boxes.Func_points = BoxesNull;
    boxes.Func_lines = BoxesNull;
    boxes.Func_triangles = BoxesNull;
    boxes.Func_squares = BoxesNull;
  }
  else
  {
    BoxesFindDomainParameters();
    BoxesMallocArrays();
    if (BF.Bool_Force_point)
    {
      printf("    Box function to place points choosen\n");
      boxes.Func_points = BoxesPlacePoints;
    }
    else                        
    {
      printf("    Box function to place points NOT choosen\n");
      boxes.Func_points = BoxesNull;
    }

    if (BF.Bool_Force_line )    
    {
      printf("    Box function to place lines choosen\n");
      boxes.Func_lines = BoxesPlaceLines;
    }
    else                        
    {
      printf("    Box function to place lines NOT choosen\n");
      boxes.Func_lines = BoxesNull;
    }

    if (BF.Bool_Force_triangle) 
    {
      printf("    Box function to place triangles choosen\n");
      boxes.Func_triangles = BoxesPlaceTriangles;
    }
    else                        
    {
      printf("    Box function to place triangles NOT choosen\n");
      boxes.Func_triangles = BoxesNull;
    }

    if (BF.Bool_Force_square)   
    {
      printf("    Box function to place squares choosen\n");
      boxes.Func_squares = BoxesPlaceSquares;
    }
    else                        
    {
      printf("    Box function to place squares NOT choosen\n");
      boxes.Func_squares = BoxesNull;
    }
  }
}


void BoxesOraganize()
{ 
  /* A function that places points, lines, triangles, squares in boxes.
    Needs 'ChooseBoxFunctions' function to be called first*/
  boxes.Func_points();
  boxes.Func_lines();
  boxes.Func_triangles();
  boxes.Func_squares();
}