void WriteVTKDataFile(const int N)
{
  /* Writes a VTK file with points and lines 
      under the file name "S_%d.vtk, N" */
  char fileName [256];
  sprintf(fileName, "%s/S_%d.vtk", vtkOutputFolder, N);
  fd = fopen(fileName,"w");
  fprintf(fd, "# vtk DataFile Version 2.0\n");
  fprintf(fd, "S, Created by Mohamad Ibrahim Cheikh\n");
  fprintf(fd, "ASCII\n");
  fprintf(fd, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(fd, "POINTS %d double\n", points.N);
  for(int i=0; i<points.N; i++)
    fprintf(fd, "%g %g %g\n", points.x[i], points.y[i], points.z[i]);
  fprintf(fd, "\n");
  fprintf(fd, "CELLS %d %d\n", lines.N, lines.N*3);
  for(int l=0; l < lines.N; l++)
    fprintf(fd, "2 %d %d\n",  lines.I0[l], lines.I1[l]);
  fprintf(fd, "\n");
  fprintf(fd, "CELL_TYPES %d\n", lines.N);
  for(int l=0; l < lines.N; l++) 
    fprintf(fd, "3\n");
  fprintf(fd, "\n");
  // Writing the forces
  fprintf(fd, "POINT_DATA %d\nVECTORS F double\n",forces.N);
  for(int i=0; i<forces.N; i++)
    fprintf(fd, "%g %g %g\n", forces.x[i], forces.y[i], forces.z[i]);
  fprintf(fd, "\n");
  // Write Cell Data
  bool cell_data_written = false;
  if (LD.Type!=L0_DYNAMICS_NULL) 
  {
    if (cell_data_written == false) { fprintf(fd, "CELL_DATA %d\n", lines.N); cell_data_written = true;}
    fprintf(fd, "SCALARS L0 double\nLOOKUP_TABLE default\n");
    for(int i=0; i<lines.N; i++) fprintf(fd, "%g\n", lines.len[i]);
    fprintf(fd, "\n"); 
  }
  if (LD.Type == L0_DYNAMICS_COND_MULTI_LAMBDA)
  {
    if (cell_data_written == false) { fprintf(fd, "CELL_DATA %d\n", lines.N); cell_data_written = true;}
    fprintf(fd, "\nSCALARS L0_M double\nLOOKUP_TABLE default\n");
    for(int i=0; i<lines.N; i++) fprintf(fd, "%g\n", LD.M_lambda[i]);
    fprintf(fd, "\n"); 
  }
  if (LE.Type == FORCE_LE_INITIAL_DISP_MULTI)
  {
    if (cell_data_written == false) { fprintf(fd, "CELL_DATA %d\n", lines.N); cell_data_written = true;}
    fprintf(fd, "\nSCALARS LE_M double\nLOOKUP_TABLE default\n");
    for(int i=0; i<lines.N; i++) fprintf(fd, "%g\n", LE.M_line[i]);
    fprintf(fd, "\n"); 
  }
  fclose(fd); // Close Flie
  printf("  VTK file %s saved\n",fileName);
}

void WriteVTKPointFile(char *Name,const int N, int *I)
{
  /* Writes a VTK file with specified points 
      N:  Number of points to write
      I:  Index of the points to write*/
  char fileName [256];
  sprintf(fileName, "%s/%s.vtk", vtkOutputFolder, Name);
  fd = fopen(fileName,"w");
  fprintf(fd, "# vtk DataFile Version 2.0\n");
  fprintf(fd, "Soft Forces, Created by Mohamad Ibrahim Cheikh\n");
  fprintf(fd, "ASCII\n");
  fprintf(fd, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(fd, "POINTS %d double\n", N);
  for(int i=0; i<N; i++)
    fprintf(fd, "%g %g %g\n", points.x[I[i]], points.y[I[i]], points.z[I[i]]);
  fprintf(fd, "\n");
  fclose(fd); // Close Flie
  printf("  VTK file %s saved\n",fileName);
}

void CreateCSVDataFile()
{
  /* Note that force reported here is the total force from the pervious time step */
  fd = fopen(csvOutputFile,"w");
  if (VC.Type == FORCE_VC_CONST_CYLIN)
    fprintf(fd,"t[ms],Vol[um^3]\n"); 
  else if (VC.Type == FORCE_VC_CONST_CELLS){   
    fprintf(fd,"t[ms],V_All"); 
    for(int c=0; c<VC.N_cells; c++) fprintf(fd,",V%d", c);
    fprintf(fd,"\n");}
  else if (VC.Type == FORCE_VC_CONST_CYLIN_CELLS){   
    fprintf(fd,"t[ms],Vol[um^3],V_All"); 
    for(int c=0; c<VC.N_cells; c++) fprintf(fd,",V%d", c);
    fprintf(fd,"\n");}
  else
    fprintf(fd,"Id,t[ms],x[um],y[um],z[um],fx[nN],fy[nN],fz[nN]\n"); 
  fclose(fd); // Close Flie
}

void AppendToCSVDataFile()
{
  /* Note that force reported here is the total force from the pervious time step */
  fd = fopen(csvOutputFile,"a");
  for (int g=0; g<boundaries.N;g++)
    if (boundaries.groups[g].Type==FORCE_ADD || boundaries.groups[g].Type==FORCE_ADD_ELLIPSE || boundaries.groups[g].Type==FORCE_SPEED_ELLIPSE)
    {
      const int bp = boundaries.groups[g].I[0];
      fprintf(fd,"%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",bp,times.real_time,points.x[bp],points.y[bp],points.z[bp],forces.x[bp],forces.y[bp],forces.z[bp]);
    }
    else if (boundaries.groups[g].Type==NO_SOLID_BC)
    {
      const int bp = boundaries.groups[g].I[0];
      fprintf(fd,"%d,%lf,%lf,%lf,%lf\n",bp,times.real_time,points.x[bp],points.y[bp],points.z[bp]);
    }
  
  if (VC.Type == FORCE_VC_CONST_CYLIN) 
    fprintf(fd,"%lf,%lf\n",times.real_time,VC.volume);
  else if (VC.Type == FORCE_VC_CONST_CELLS)
  {
    fprintf(fd,"%lf,%lf",times.real_time,VC.volume_cells);
    for(int c=0; c<VC.N_cells; c++) fprintf(fd, ",%lf", VC.cell_V[c]);
    fprintf(fd,"\n");
  }
  else if (VC.Type == FORCE_VC_CONST_CYLIN_CELLS)
  {
    fprintf(fd,"%lf,%lf,%lf",times.real_time,VC.volume,VC.volume_cells);
    for(int c=0; c<VC.N_cells; c++) fprintf(fd, ",%lf", VC.cell_V[c]);
    fprintf(fd,"\n");
  }
  fclose(fd); // Close Flie
}

int GetFileSize(const char* fileName)
{
  int N = 0;
  fd   = fopen(fileName,"r");
  if (fd)
  {
    while (1)
    {
      char *ptr = fgets(line, 256, fd);
      if (!ptr) line[0] = 0;

      const int len = strlen(line);
      if ( len==0)        break; // No data was found on that line
      if ( line[0] == 0)  break; // End of file
      N++;
    }
    fclose(fd);
  }
  else
  {
    printf("    ERROR!! File '%s' not found\n", fileName);
    exit(EXIT_FAILURE);
  }
  return N;
}


void ReadDoubleFile(const char* fileName, const int data_size, double *data)
{
  fd   = fopen(fileName,"r");
  if (fd)
  {
    for(int i=0; i<data_size; i++)
    {
      char *ptr = fgets(line, 256, fd);
      if (!ptr) 
        line[0] = 0;

      const int len = strlen(line);
      if ( len==0)        break; // No data was found on that line
      if ( line[0] == 0)  break; // End of file
      data[i] = atof(line);
      //printf("x[%d]=%lf\n",i,x[i]);
    }
    fclose(fd); // Close Flie
    printf("  File '%s' read\n", fileName);
  }
  else
  {
    printf("    ERROR!! File '%s' not found\n", fileName);
    exit(EXIT_FAILURE);
  }
}


void ReadIntFile(const char* fileName, const int data_size, int *data)
{
  fd   = fopen(fileName,"r");
  if (fd)
  {
    for(int i=0; i<data_size; i++)
    {
      char *ptr = fgets(line, 256, fd);
      if (!ptr) 
        line[0] = 0;

      const int len = strlen(line);
      if ( len==0)        break; // No data was found on that line
      if ( line[0] == 0)  break; // End of file
      //data[i] = atoi(line);     // Give problems if string is in sci notation
      data[i] = (int) atof(line); 
    }
    fclose(fd); // Close Flie
    printf("  File '%s' read\n", fileName);
  }
  else
  {
    printf("    ERROR!! File '%s' not found\n", fileName);
    exit(EXIT_FAILURE);
  }
}

void ReadPointsFromMeshFolder(const char meshInputFolder[256])
{
  char fileToRead[240];

  printf("Reading points from mesh folder '%s'\n", meshInputFolder);

  sprintf(fileToRead, "./%s/x.dat", meshInputFolder);

  points.N = GetFileSize(fileToRead);

  printf("  Number of points: %d\n", points.N);

  points.x = (double *)malloc(points.N*sizeof(double));
  points.y = (double *)malloc(points.N*sizeof(double));
  points.z = (double *)malloc(points.N*sizeof(double));

  sprintf(fileToRead, "./%s/x.dat", meshInputFolder);
  ReadDoubleFile(fileToRead, points.N, points.x);

  sprintf(fileToRead, "./%s/y.dat", meshInputFolder);
  ReadDoubleFile(fileToRead, points.N, points.y);
  if (DIM ==3)
  {
    sprintf(fileToRead, "./%s/z.dat", meshInputFolder);
    ReadDoubleFile(fileToRead, points.N, points.z);
  }
  else{
    for(int i=0; i<points.N; i++)   points.z[i] = 0.0;  }
}

void ReloadPointsFromMeshFolder(const char reloadFolder[256])
{
  char fileToRead[240];

  printf("Reloading data from reload folder '%s'\n", reloadFolder);   

  sprintf(fileToRead, "./%s/x_%d.dat", reloadFolder, times.reload_step);
  ReadDoubleFile(fileToRead, points.N, points.x);

  sprintf(fileToRead, "./%s/y_%d.dat", reloadFolder, times.reload_step);
  ReadDoubleFile(fileToRead, points.N, points.y);

  if (DIM ==3){
    sprintf(fileToRead, "./%s/z_%d.dat", reloadFolder, times.reload_step);
    ReadDoubleFile(fileToRead, points.N, points.z);}
}


void ReloadL0FromMeshFolder(const char reloadFolder[256])
{ 
  char fileToRead[240];

  printf("Reloading L0 data from reload folder '%s'\n", reloadFolder);   

  sprintf(fileToRead, "./%s/l0_%d.dat", reloadFolder, times.reload_step);

  const int N = GetFileSize(fileToRead);
  if (N != lines.N) {printf("    ERROR!! ReloadL0FromMeshFolder is reading L0 with size '%d', the size should be '%d'\n", N, lines.N);   exit(EXIT_FAILURE);}
  ReadDoubleFile(fileToRead, lines.N, lines.len);
}


void ReadCellsFromMeshFolder(const char meshInputFolder[256])
{

  char fileToRead[240];

  printf("Reading cells from mesh folder '%s'\n", meshInputFolder);

  sprintf(fileToRead, "./%s/edge1.dat", meshInputFolder);
  lines.N = GetFileSize(fileToRead);

  printf("Number of edges: %d\n", lines.N);

  lines.I0 = (int *)malloc(lines.N*sizeof(int));
  lines.I1 = (int *)malloc(lines.N*sizeof(int));

  sprintf(fileToRead, "./%s/edge1.dat", meshInputFolder);
  ReadIntFile(fileToRead, lines.N, lines.I0);

  sprintf(fileToRead, "./%s/edge2.dat", meshInputFolder);
  ReadIntFile(fileToRead, lines.N, lines.I1);

  lines.len = (double *)malloc(lines.N*sizeof(double));
  for(int l=0; l<lines.N; l++)
    lines.len[l] = CalcLineLength(l);

  if (print_mesh_lines)
    for(int l=0; l<lines.N; l++)
      printf("%d: (%d, %d) L0:%lf\n",l, lines.I0[l], lines.I1[l], lines.len[l]);

  //---------------------------------------------------------------------------------

  if (BF.Bool_Force_triangle || ST.Type == FORCE_ST_CELLS || (VC.Type != FORCE_VC_NULL && DIM == 3) )
  {
    sprintf(fileToRead, "./%s/triangle1.dat", meshInputFolder);
    triangles.N = GetFileSize(fileToRead);

    printf("Number of Triangles: %d\n", triangles.N);

    triangles.I0 = (int *)malloc(triangles.N*sizeof(int));
    triangles.I1 = (int *)malloc(triangles.N*sizeof(int));
    triangles.I2 = (int *)malloc(triangles.N*sizeof(int));

    sprintf(fileToRead, "./%s/triangle1.dat", meshInputFolder);
    ReadIntFile(fileToRead, triangles.N, triangles.I0);

    sprintf(fileToRead, "./%s/triangle2.dat", meshInputFolder);
    ReadIntFile(fileToRead, triangles.N, triangles.I1);

    sprintf(fileToRead, "./%s/triangle3.dat", meshInputFolder);
    ReadIntFile(fileToRead, triangles.N, triangles.I2);

    triangles.area = (double *)malloc(triangles.N*sizeof(double));
    for(int t=0; t<triangles.N; t++) triangles.area[t] = CalcTriangleArea2(t);
  }

  //---------------------------------------------------------------------------------

  if (BF.Bool_Force_square)
  {
    sprintf(fileToRead, "./%s/square1.dat", meshInputFolder);
    squares.N = GetFileSize("./meshfiles/square1.dat");

    printf("Number of Triangles: %d\n", squares.N);

    squares.I0 = (int *)malloc(squares.N*sizeof(int));
    squares.I1 = (int *)malloc(squares.N*sizeof(int));
    squares.I2 = (int *)malloc(squares.N*sizeof(int));
    squares.I3 = (int *)malloc(squares.N*sizeof(int));

    sprintf(fileToRead, "./%s/square1.dat", meshInputFolder);
    ReadIntFile(fileToRead, squares.N, squares.I0);

    sprintf(fileToRead, "./%s/square2.dat", meshInputFolder);
    ReadIntFile(fileToRead, squares.N, squares.I1);

    sprintf(fileToRead, "./%s/square3.dat", meshInputFolder);
    ReadIntFile(fileToRead, squares.N, squares.I2);

    sprintf(fileToRead, "./%s/square4.dat", meshInputFolder);
    ReadIntFile(fileToRead, squares.N, squares.I3);

    squares.area = (double *)malloc(squares.N*sizeof(double));
    for(int t=0; t<squares.N; t++)
      squares.area[t] = CalcSquareArea(t);
  }

}


void ReadPointsFromVTKFile(const char VTK_File[256])
{ 
  /* Reads the 'POINTS' section in the vtk file and creats the coodinate.
      Variables Created or Initialized are:
        points.N
        points.x
        points.y
        points.z
  */
  int k;
  char **data;
  char line[256]="";
  bool flg_POINTS = false;

  fd  = fopen(VTK_File,"r");

  if(fd == NULL) //if file does not exist, create it
  { 
    printf("    ERROR!! File '%s' not found\n", VTK_File);
    exit(EXIT_FAILURE);
  }

  printf("Reading POINTS section in VTK file '%s'\n", VTK_File);
  // Start reading
  while (true)
  {
    ReadFileLine(fd, 256,line); //printf("%s\n",line);
    StrToArray(line,' ', &k, &data); //printf("%s\n",data[0]);
    if (strcmp(data[0],"POINTS")==0) flg_POINTS = true;
    StrToArrayDestroy(k, data);

    if(flg_POINTS)
    { 
      StrToArray(line,' ', &k, &data);
      points.N = atoi(data[1]);
      StrToArrayDestroy(k, data);
      printf("  Number of points: %d\n", points.N);

      points.x = (double *)malloc(points.N*sizeof(double));
      points.y = (double *)malloc(points.N*sizeof(double));
      points.z = (double *)malloc(points.N*sizeof(double));


      // Saving points:
      for(int i=0; i < points.N; i++)
      {
        ReadFileLine(fd, 256,line);    // Get the next line  
        StrToArray(line,' ', &k, &data);
        points.x[i] = atof(data[0]);
        points.y[i] = atof(data[1]);
        points.z[i] = atof(data[2]);
        StrToArrayDestroy(k, data);
      }
      break;
    }

    if (line[0] == 0)  break; // End of file
  }
  // Close VTK Flie
  fclose(fd);

  printf("  Done reading POINTS section in VTK file\n");
  if (flg_POINTS == false) { printf("    ERROR!! No point was read from VTK file\n"); exit(EXIT_FAILURE);  }
}


void ReadCellsFromVTKFile(const char VTK_File[256])
{  
  /* Reads the 'CELL_TYPES', and 'CELLS' section in the vtk file.
      Part 1- Counts the number of line cells, triangle cells
          Variables Initialized are:
            lines.N
            triangles.N

    Part 3-  Inserts the index of the points in the line and triangle cells:
          Variables Initialized are:
            lines.I0
            lines.I1
            triangles.I0
            triangles.I1
            triangles.I2
  */
  char line[256]="";

  //const int VTKPOINTS = 1;
  const int VTKLINES = 3;
  const int VTKTRIANGLES = 5;

  // Boolean variable for each part of the VTK file
  bool flg_CELLS = false, flg_CELL_TYPES = false;

  int  cells_total, cells_types_total;
  int  k;
  char **data;

  printf("Reading Cells from VTK file '%s'\n", VTK_File);
  //------------------------------------------------------------------------------//
  // Part 1: Counting the cell types(lines or triangles)
  //------------------------------------------------------------------------------//
  /* NOTE: One can make the code faster by creating an array of size cells_types_total
           and saving the cell types inside */

  // Opening VTK file on rank 0
  int cell_t;

  lines.N = 0;      // Number of lines
  triangles.N = 0;  // Number of triangles

  fd   = fopen(VTK_File,"r");

  // Start reading
  while (true)
  {
    // Read line
    ReadFileLine(fd, 256,line); 
    if (strncmp(line, "CELL_TYPES", 10) == 0) flg_CELL_TYPES = true;

    if(flg_CELL_TYPES)
    {
      StrToArray(line,' ',&k, &data);
      cells_types_total = atoi(data[1]);
      StrToArrayDestroy(k, data);

      // Saving cell type:
      for(int i=0; i < cells_types_total; i++)
      {
        ReadFileLine(fd, 256,line);                                    // Get the next line  
        StrToArray(line,' ',&k, &data);
        cell_t = atoi(data[0]);
        if (cell_t == VTKLINES)           lines.N++;
        else if (cell_t == VTKTRIANGLES)  triangles.N++;
        StrToArrayDestroy(k, data);
      }
      break;
    }
    if (line[0] == 0)  break; // End of file
  }
  // Close VTK Flie
  fclose(fd);
  printf("  Done reading CELL_TYPES section in VTK File\n");
  printf("    Number of lines = %d\n    Number of triangles = %d\n", lines.N, triangles.N);
  if (!flg_CELL_TYPES) { printf("    ERROR!! No cell type was read from VTK file\n"); exit(EXIT_FAILURE);  }


  //------------------------------------------------------------------------------//
  // Part 2: Saving the points index in the line and triangle cells
  //------------------------------------------------------------------------------//

  lines.I0 = (int *)malloc(lines.N*sizeof(int));
  lines.I1 = (int *)malloc(lines.N*sizeof(int));
  triangles.I0 = (int *)malloc(triangles.N*sizeof(int));
  triangles.I1 = (int *)malloc(triangles.N*sizeof(int));
  triangles.I2 = (int *)malloc(triangles.N*sizeof(int));

  // Opening VTK file on rank 0
  int cell_size, l=0, t=0;

  fd   = fopen(VTK_File,"r");

  // Start reading
  while (true)
  {
    ReadFileLine(fd, 256,line);  // Read line
    if( strncmp(line, "CELLS",  5) == 0) flg_CELLS = true;

    if(flg_CELLS)
    {
      StrToArray(line,' ',&k, &data);
      cells_total = atoi(data[1]);
      StrToArrayDestroy(k, data);

      // Saving Cells
      for(int i=0; i < cells_total; i++)
      {
        ReadFileLine(fd, 256,line);                                    // Get the next line  
        StrToArray(line,' ',&k, &data);
        cell_size = atoi(data[0]);
        if (cell_size == 2) { // Line
          lines.I0[l] = atoi(data[1]);
          lines.I1[l] = atoi(data[2]);
          l++; }
        else if (cell_size == 3) { // Triangles
          triangles.I0[t] = atoi(data[1]);
          triangles.I1[t] = atoi(data[2]);
          triangles.I2[t] = atoi(data[3]);
          t++;}
        StrToArrayDestroy(k, data);
      }
      break;
    }
    if (line[0] == 0)  break; // End of file
  }
  // Close VTK Flie
  fclose(fd);
  printf("  Done reading CELLS section in VTK File\n");
  if (!flg_CELLS) { printf("    ERROR!! No cell was read from VTK file\n"); exit(EXIT_FAILURE);  }


  //------------------------------------------------------------------------------//
  // Part 3: Extra
  //------------------------------------------------------------------------------//

  if (lines.N>0){
    lines.len = (double *)malloc(lines.N*sizeof(double));
    for(int l=0; l<lines.N; l++)  lines.len[l] = CalcLineLength(l);  }

  if (triangles.N>0) {
    triangles.area = (double *)malloc(triangles.N*sizeof(double));
    for(int t=0; t<triangles.N; t++) triangles.area[t] = CalcTriangleArea2(t);  }

}


void ReloadPointsFromVTKFile(const char VTK_File[256])
{ 
  /* Reads the 'POINTS' section in the vtk file and creats the coodinate.  */
  int k;
  char **data;
  char line[256]="";
  bool flg_POINTS = false;

  fd  = fopen(VTK_File,"r");

  if(fd == NULL) //if file does not exist, create it
  { 
    printf("    ERROR!! File '%s' not found\n", VTK_File);
    exit(EXIT_FAILURE);
  }

  printf("Reloading data from POINTS section in VTK file '%s'\n", VTK_File);
  // Start reading
  while (true)
  {
    ReadFileLine(fd, 256,line); //printf("%s\n",line);
    StrToArray(line,' ', &k, &data); //printf("%s\n",data[0]);
    if (strcmp(data[0],"POINTS")==0) flg_POINTS = true;
    StrToArrayDestroy(k, data);

    if(flg_POINTS)
    { 
      StrToArray(line,' ', &k, &data);
      const int N = atoi(data[1]);
      StrToArrayDestroy(k, data);
      if (N != points.N){
        printf("    ERROR!! Number of points in reload file is '%d' not consistant with original file '%d'\n", N, points.N);
        exit(EXIT_FAILURE); }

      // Reloading points:
      for(int i=0; i < points.N; i++)
      {
        ReadFileLine(fd, 256,line);    // Get the next line  
        StrToArray(line,' ', &k, &data);
        points.x[i] = atof(data[0]);
        points.y[i] = atof(data[1]);
        points.z[i] = atof(data[2]);
        StrToArrayDestroy(k, data);
      }
      break;
    }

    if (line[0] == 0)  break; // End of file
  }
  // Close VTK Flie
  fclose(fd);

  printf("  Done reading POINTS section in VTK file\n");
  if (flg_POINTS == false) { printf("    ERROR!! No point was read from VTK file\n"); exit(EXIT_FAILURE);  }
}



void ReloadL0FromVTKFile(const char VTK_File[256])
{ 
  /* Reads the 'CELL_DATA' section in the vtk file and looks for the variable 'Cell_Data_Name'.  */
  int k;
  char **data;
  char line[256]="";
  bool flg_CELL_DATA = false;

  fd  = fopen(VTK_File,"r");

  if(fd == NULL) //if file does not exist, create it
  { 
    printf("    ERROR!! File '%s' not found\n", VTK_File);
    exit(EXIT_FAILURE);
  }

  printf("Reloading data from 'L0' data from the CELL_DATA section in VTK file '%s'\n",  VTK_File);
  // Start reading
  while (true)
  {
    ReadFileLine(fd, 256,line); //printf("%s\n",line);
    StrToArray(line,' ', &k, &data); //printf("%s\n",data[0]);
    if (strcmp(data[0],"CELL_DATA")==0) flg_CELL_DATA = true;
    StrToArrayDestroy(k, data);

    if(flg_CELL_DATA)
    { 
      StrToArray(line,' ', &k, &data);
      const int N = atoi(data[1]); // Number of lines
      StrToArrayDestroy(k, data);
      if (N != lines.N){
        printf("    ERROR!! Number of 'CELL_DATA' variables in reload file is '%d' not consistant with original file '%d'\n", N, lines.N);
        exit(EXIT_FAILURE); }

      ReadFileLine(fd, 256,line);    // Get the next line  
      StrToArray(line,' ', &k, &data);
      if ( strcmp(data[0],"SCALARS")!=0 || strcmp(data[1],"L0")!=0 ) {
        printf("    ERROR!! Line in 'CELL_DATA' section not found. Line looking for is 'SCALARS L0 double' variable found '%s'\n", line);
        exit(EXIT_FAILURE); }
      StrToArrayDestroy(k, data);

      ReadFileLine(fd, 256,line);    // Get the next line  
      StrToArray(line,' ', &k, &data);
      if ( strcmp(data[0],"LOOKUP_TABLE")!=0 || strcmp(data[1],"default")!=0 ) {
        printf(" %d-%d\n", strcmp(data[0],"LOOKUP_TABLE")!=0, strcmp(data[1],"default")!=0 );
        printf("    ERROR!! Line in 'CELL_DATA' section not found. Line looking for is 'LOOKUP_TABLE default' variable found '%s'\n", line);
        exit(EXIT_FAILURE); }
      StrToArrayDestroy(k, data);

      // Reloading points:
      for(int i=0; i < N; i++)
      {
        ReadFileLine(fd, 256,line);    // Get the next line  
        StrToArray(line,' ', &k, &data);
        lines.len[i] = atof(data[0]);
        StrToArrayDestroy(k, data);
      }
      break;
    }

    if (line[0] == 0)  break; // End of file
  }
  // Close VTK Flie
  fclose(fd);

  printf("  Done reading 'L0' portion of the 'CELL_DATA' section in VTK file\n");
  if (flg_CELL_DATA == false) { printf("    ERROR!! No 'L0' portion of the 'CELL_DATA' section was read from VTK file\n"); exit(EXIT_FAILURE);  }
}
