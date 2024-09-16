// Defining a function pointers
typedef void   (*BasicFunctions)(); 
typedef double (*DoubleFunctions)(); 
typedef void   (*BoundaryFunctions)(void *v);
typedef double (*CalculateBarrierForceFunction)(double D, const double D_min, const double C1, const double C2);

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------

typedef struct 
{
  double real_time;       // Actual time in simulation [ms]
  double dt;              // Time Step size [ms]
  double end_time;        // End time [ms]
  int    time_step;       // Time Step counter

  bool use_real_time_to_save;   
  
  double dtime_save_vtk;  // Time Step size to save vtk data [ms]
  double dtime_save_csv;  // Time Step size to save csv data [ms]
  double dtime_print;     // Time Step size to print progress [ms]

  double time_save_vtk;   // Time to save vtk data [ms]
  double time_save_csv;   // Time to save csv data [ms]
  double time_print;      // Time to print progress [ms]

  int ntimeSteps_save_vtk; // number of time steps save vtk data [ms]
  int ntimeSteps_save_csv; // number of time steps save csv data [ms]
  int ntimeSteps_print;    // number of time steps print progress [ms]

  double tstart, ttotal;  // Computational time [s]

  bool   reload;
  double reload_time;
  int    reload_step;
} Time_Struct; 
//---------------------------------------------------------------------------------------------

typedef struct 
{
  int     N;    // Number of total points
  double  *x,  *y,  *z;       // Coordinate of each point 
} Points_Struct;

//---------------------------------------------------------------------------------------------

typedef struct 
{
  SOLID_BC_TYPE Type;
  char    *Name;            // Boundary Name
  int      N;               // Number of total boundary points
  int     *I;               // Index of boundary points 
  double  *x,  *y,  *z;     // Coordinate of each point [if needed]
  double   A,   B;          // Ellipse info [If BP is FORCE_ADD_ELLIPSE]
  int      Dir;             // Direction
  double   Dtheta0;         // Dtheta0 will be used to add an initial jump to create a fake set of solidPoints [FORCE_SPEED_ELLIPSE]
  double   DthetaDt;        // DthetaDt will be used to be the velocity [FORCE_SPEED_ELLIPSE]
  double   K;               // K is the constant for the difference between fake and real points [FORCE_SPEED_ELLIPSE]
  double   *theta;          // Coordinate array to be used by the boundary [FORCE_SPEED_ELLIPSE]
  //  This array will create a fake set of solidPoints infront of solid Points with boundary [FORCE_SPEED_ELLIPSE]
  //  These solid Points will be tracked
  double   *Z_Plane;        // The z-plane on which the fake points reside [FORCE_SPEED_ELLIPSE]
  double   Fx,  Fy,  Fz;    // Forces
  double   Fm;              // Force magnitude

  char    ExternalForceFile [256]; // File to write total external forces on
  double Total_External_Force[3];  // total external forces

  BoundaryFunctions  Force_Func;
  BoundaryFunctions   Move_Func;
  BoundaryFunctions Report_Func;
} Boundary_Points_Struct; 

//---------------------------------------------------------------------------------------------

typedef struct 
{
  int N;                            // Number of total boundary groups
  Boundary_Points_Struct *groups;
} Boundary_Groups_Struct; 

//---------------------------------------------------------------------------------------------

typedef struct 
{
  int        N;   // Number of total lines
  int      *I0;   // Index first points
  int      *I1;   // Index second points
  double  *len;   // Length of each line element (Used to store initial length)
} Lines_Struct; 

//---------------------------------------------------------------------------------------------

typedef struct 
{
  int        N;   // Number of total trianguls
  int      *I0;   // Index first points
  int      *I1;   // Index second points
  int      *I2;   // Index third points
  double  *area;  // Area of each triangluar element (Used to store initial area)
} Triangles_Struct; 

//---------------------------------------------------------------------------------------------
typedef struct 
{
  int        N;   // Number of total squares
  int      *I0;   // Index first points
  int      *I1;   // Index second points
  int      *I2;   // Index thid points
  int      *I3;   // Index thid points
  double  *area;  // Area of each square element (Used to store initial area)
} Squares_Struct; 

//---------------------------------------------------------------------------------------------

typedef struct 
{
  BasicFunctions Func_points;     // Function for points
  BasicFunctions Func_lines;      // Function for lines
  BasicFunctions Func_triangles;  // Function for triangles
  BasicFunctions Func_squares;    // Function for squares

  int     N;              // Number of boxes
  int     Max_Points;     // Number of MAX points per box
  int     Max_Lines;      // Number of MAX lines per box
  int     Max_Triangles;  // Number of MAX triangles per box
  int     Max_Squares;    // Number of MAX squares per box

  // Organization of elements according to different domains
  int **box_points;     // An array that contains a set of points that are found at a specific location
  int **box_lines;      // An array that contains a set of lines that are found at a specific location
  int **box_triangles;  // An array that contains a set of triangles that are found at a specific location
  int **box_squares;    // An array that contains a set of squares that are found at a specific location
    // These array divides the domain into a set of boxes (or squares if 2D) and places 
    // elements in the box depending on the centroid of the element or edges.
    // The row "r" represent the box number which is r = xb + yb*Nx + zb*Nx*Ny
    //    (xb, yb, zb) represent the center of the box
    //    (Nx, Ny, Nz) represent the number of boxes in each direction where
    //     Nx = Lx / dx and Lx is the domain length
    // The first column "c0" represent the number of elements in that box
    // The other columns "c1.." represent the local element index

  double M;          // Multiplier for the average box size
  int    Nx, Ny, Nz; // Number of boxes in x, y and z
  double Lx, Ly, Lz; // Length of total domain in x, y and z
  double dx, dy, dz; // size of boxes in x, y and z
  double Xmin, Xmax; // 
  double Ymin, Ymax; // 
  double Zmin, Zmax; // 

} Boxes_Struct; 

//---------------------------------------------------------------------------------------------

typedef struct 
{
  L0_DYNAMICS_TYPE Type;
  BasicFunctions   Func;

  double     lambda;      // Stretching Constant to activate L-Node Dynamics (i.e. If L > (1+lambda)*L0 ==> L-Node Dynamics activated)
  double     tau;         // Time Constant for L-Node Dynamics (i.e tau*dL0/dt = (L-L0))
  double    *M_lambda;    // Multiplier for lambda
  bool      *Bool_LNode;  // An array of size equal to num of lines (if true ==> this line is in L-Node dynamic mode)
  char       M_lambda_file_name[200];
} L0Dynamic_Struct; 

//---------------------------------------------------------------------------------------------

typedef struct 
{
  FORCE_LE_TYPE  Type;
  BasicFunctions Func;
  double  K;            // Elastic Spring constant for all lines
  double *M_line;       // Elastic Spring constant multiplier for each line (i.e K_line = K*M_line)
  char    M_file_name[200];
} Forces_LE_Struct; 

//---------------------------------------------------------------------------------------------
// Defining a function pointer to calculate the Barrier force

typedef struct  {
  FORCE_BF_TYPE     Type;
  BasicFunctions    Func;

  int               N; // Number of points to apply barrier force (Only FORCE_BF_SET_POINTS)
  int              *I; // Index of points to add barrier force  (Only FORCE_BF_SET_POINTS)

  bool                          Bool_Force_line;
  FUNCTION_TYPE                 Type_Force_line;
  CalculateBarrierForceFunction Calc_Force_line;
  double D_min_line;                 // minimum distance to add barrier force
  double C1_line,  C2_line;          // The Barrier force function constants for lines

  bool                          Bool_Force_point;
  FUNCTION_TYPE                 Type_Force_point;
  CalculateBarrierForceFunction Calc_Force_point;
  double D_min_point;               // minimum distance to add barrier force
  double C1_point, C2_point;        // The Barrier force function constants for points

  bool                          Bool_Force_triangle;
  FUNCTION_TYPE                 Type_Force_triangle;
  CalculateBarrierForceFunction Calc_Force_triangle;
  double D_min_triangle;               // minimum distance to add barrier force
  double C1_triangle, C2_triangle;        // The Barrier force function constants for points

  bool                          Bool_Force_square;
  FUNCTION_TYPE                 Type_Force_square;
  CalculateBarrierForceFunction Calc_Force_square;
  double D_min_square;               // minimum distance to add barrier force
  double C1_square, C2_square;        // The Barrier force function constants for points

} Forces_BF_Struct; // Barrier Forces

//---------------------------------------------------------------------------------------------
typedef struct  {

  FORCE_SF_TYPE     Type;   // Soft Force Type
  BasicFunctions    Func;

  int                  N;    // Number of points to add soft force
  int                 *p;    // Index of points to add soft force
  double         A, B, C;    // The major & minor axis of the single Ellipsoid SF method 
  double          Bd, Bv;    // The minor axis of the two Ellipsoid SF method 
  double              ds;    // Tolarance on the SF method 
  double              sf;    // Magnitude of the soft force

} Forces_SF_Struct; // Soft Forces

//---------------------------------------------------------------------------------------------
typedef struct  {

  FORCE_ST_TYPE     Type;   // Surface Tension Force Type
  BasicFunctions    Func;

  double                 gamma;     // Magnitude of the surface tension

  // FORCES_ST_MULT:
  double                *M_line;

  // FORCES_ST_CELL:
  int                  N_Cells;     // Number of cells
  double                *A,*A0;     // Total Area & Original Total Area of each Cell (Size = N_Cells)
  int              *Cell_Index;     // Cell index for each triangle (Size = Triangles_Struct.N)

} Forces_ST_Struct; // Surface Tension Forces

//---------------------------------------------------------------------------------------------
typedef struct  {

  FORCE_VC_TYPE     Type;     // Volume Conservation Type
  BasicFunctions    Func;

  double        volume_0;     // Initial volume
  double          volume;     // Current volume
  double       center[3];     // Center of the volume of whole domain

  // Parameters for FORCE_VC_CONST_CYLIN
  double           alpha_cylin;       // Magnitude of volume conservation force of the cylinder
  int       N_cylin_triangles;        // Number of triangles on the cylinder surface that will determin cylinder volumes
  int         *cylin_triangles;       // triangles on the cylinder surface used to calculate the cylinder volumes
  int       N_cylin_points;           // Number of points on the cylinder surface that will determin the shape of the two circles
  int         *cylin_points;          // points on the cylinder surface used to calculate area of circles for cylinder volumes
  double    cylin_circle1_point[3];   // Center Point of the first cylinder circle
  double    cylin_circle1_normal[3];  // Normal of the first cylinder circle
  double    cylin_circle2_point[3];   // Center Point of the second cylinder circle
  double    cylin_circle2_normal[3];  // Normal of the second cylinder circle

  // Paremeter for FORCE_VC_CONST_CELLS
  double    volume_cells;                // Total volume of cells
  double    alpha_cells;                 // Magnitude of volume conservation force of the cylinder
  int       N_cells;                     // Number of cells
  double   *cell_V,*cell_V0;             // Volume & Original Volume of each cell (Size = N_Cells)
  int      *cell_traingle_index;         // Cell index for each triangle (Size = Triangles_Struct.N)
  double   *cell_center_x, *cell_center_y, *cell_center_z; // Total Volume & Original Total Volume of each Cell (Size = N_Cells)
  int      *cell_num_tri;                // Number of triangles for each cell


} Forces_VC_Struct; // Volume Conservation Forces

//---------------------------------------------------------------------------------------------

typedef struct {

  DYNAMIC_MESH_TYPE   Type;    // Dynamic Mesh Type
  BasicFunctions   Func;
  double              L_max;  // Maximum Length before dividing lines
  double              L_min;  // Minimum Length before removing lines
  
} Dynamic_Mesh_Struct; // Dynamic Mesh Forces

//---------------------------------------------------------------------------------------------

typedef struct 
{
  CONSTRAIN_TYPE    Type;         // Method of constaint  
  BasicFunctions    Func_create;  // Function to create constrained points
  BasicFunctions    Func_points;  // Function to constrain points
  BasicFunctions    Func_forces;  // Function to fix forces of constrained points
  BasicFunctions    Func_free;    // Function to free constrained points

  int                  N;    // Number of points constrained
  int                 *p;    // Index of points to constrain
  int             axis_i;    // The index of the axis to constrain (0=x,1=y,2=z)
  double          axis_v;    // The constrain value of the axis 

  int                  N2;   // Number of 2nd points constrained
  int                 *p2;   // 2nd Index of points to constrain
  int             axis_i2;   // The 2nd index of the axis to constrain (0=x,1=y,2=z)
  double          axis_v2;   // The 2nd constrain value of the axis 

  int                  N3;   // Number of 3rd points constrained
  int                 *p3;   // 3rd Index of points to constrain
  int             axis_i3;   // The 3rd index of the axis to constrain (0=x,1=y,2=z)
  double          axis_v3;   // The 3rd constrain value of the axis 

  double              ds;    // Tolarance on the Constrain method 

} Constrain_Struct; // Constrains

//--------------------------------------------------------------------------------------------------

typedef struct 
{
  double c;     // Drag coefficient
} Drag_Struct;  // Drag 