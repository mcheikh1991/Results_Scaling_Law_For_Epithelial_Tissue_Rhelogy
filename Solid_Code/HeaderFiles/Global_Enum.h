
// Solid Boundary Condition Type
typedef enum {NO_SOLID_BC=0, FORCE_ADD=1, FORCE_INSERT=2, DIST_FIXED=3, DIST_SPEED=4, FORCE_SPEED=5, SYMMETRY=6 , 
  FORCE_SPEED_ELLIPSE=7, DIST_SPEED_ELLIPSE=8, DIST_SPEED_X_LINEAR_Y=9, AXIS_FIXED=10, FORCE_ADD_ELLIPSE=11} SOLID_BC_TYPE;
const char* SOLID_BC_TYPE_CHAR[] = {"NO_SOLID_BC", "FORCE_ADD", "FORCE_INSERT", "DIST_FIXED", "DIST_SPEED", 
  "FORCE_SPEED", "SYMMETRY", "FORCE_SPEED_ELLIPSE", "DIST_SPEED_ELLIPSE", "DIST_SPEED_X_LINEAR_Y", "AXIS_FIXED", "FORCE_ADD_ELLIPSE"};
/* Info about the SOLID BC:
  0- NO_SOLID_BC: No BC
  1- FORCE_ADD:     Add boundary force at boundary points
  2- FORCE_INSERT:  Insert boundary force at boundary points
  3- DIST_FIXED:  Fix the distance (do not move boundary points)
  4- DIST_SPEED:  Move the points by a prescribed boundary speed
  5- FORCE_SPEED: Create fake points and have them be moved at a prescribed boundary speed. Attach 
                  to each of these fake points a boundary point with an elastic connection. and have it be 
                  dragged behind theses fake point by a large elastic spring. 
  6- SYMMETRY:    Tries to enforce symmetry using one of the axis as a symmetry plane
  7- FORCE_SPEED_ELLIPSE: Create fake point along an ellipse and attach them to boundary points with an elastic spring.
                          Next step move the fake point and transfer the force to the boundary points 
  8- DIST_SPEED_ELLIPSE:  Move boundary points along an ellipse in a constant speed
  9- DIST_SPEED_X_LINEAR_Y: Move boundary points in the x direction with a velocity that varies linearly with y
  10- AXIS_FIXED:           Fix the position of one axis of the boundary point
  11- FORCE_ADD_ELLIPSE: Add a force that is tanget to the surface of an ellipse
*/

// Function Type
typedef enum {FUNC_CONSTANT=0, FUNC_LINEAR=1, FUNC_EXP=2} FUNCTION_TYPE;
const char* FUNCTION_TYPE_CHAR[] = {"Constant", "Linear", "Exponential"};


// L0 Dynamics:
typedef enum {L0_DYNAMICS_NULL=0, L0_DYNAMICS_ALWAYS=1, L0_DYNAMICS_COND=2, L0_DYNAMICS_COND_MULTI_LAMBDA=3} L0_DYNAMICS_TYPE;
const char* L0_DYNAMICS_CHAR[] = {"NO_L0_DYNAMICS", "L0_DYNAMICS_ALWAYS", "L0_DYNAMICS_COND", "L0_DYNAMICS_COND_MULTI_LAMBDA"};
/* Info about the FLUID_VISCOSITY_METHOD:
  0- L0_DYNAMICS_NULL:                No L0 Dynamics
  1- L0_DYNAMICS_ALWAYS:              L0 Dynamics physics are always on (i.e. L0_new = dt/tau*(L-L0_old) + L0_old)
  2- L0_DYNAMICS_COND:                L0 Dynamics physics are activated if L > L0*(1+lambda)
  3- L0_DYNAMICS_COND_MULTI_LAMBDA:   L0 Dynamics physics are activated if L > L0*(1+lambda) with lambda being multiple
  */


// Forces: Barrier Force Models
typedef enum {FORCE_BF_NULL=0, FORCE_BF_ALL_POINTS=1, FORCE_BF_SET_POINTS=2} FORCE_BF_TYPE;
const char* FORCE_BF_CHAR[] = {"NO_BARRIER_FORCE", "BARRIER_FORCE_ALL_POINTS", "BARRIER_FORCE_SET_POINTS"};
/* Info about the SOLID_FORCE_BF_TYPE:
  0- FORCE_BF_NULL:          No Barrier Force
  1- FORCE_BF_ALL_POINTS:    Barrier Force on all points
  2- FORCE_BF_SET_POINTS:    Barrier Force on set of points
  */

// Forces: Linear Elasticity Models
typedef enum {FORCE_LE_NULL=0, FORCE_LE_FIXED_DISP=1, FORCE_LE_INITIAL_DISP=2, FORCE_LE_INITIAL_DISP_MULTI=3} FORCE_LE_TYPE;
const char* FORCE_LE_CHAR[] = {"NO_LINEAR_ELASTIC_FORCE", "LINEAR_ELASTIC_FORCE_FIXED_LENGTH", "LINEAR_ELASTIC_FORCE_INITIAL_LENGTH", "LINEAR_ELASTIC_FORCE_MULTIPLE_INITIAL_LENGTH"};
/* Info about the Solid Internal Forces Type:
  0- FORCE_LE_NULL:                     No Forces between points
  1- FORCE_LE_FIXED_DISP:               Linear Elastic Force = K*(D-D0) where D0 is from input saved in system
  2- FORCE_LE_INITIAL_DISP:             Linear Elastic Force = K*(D-D0) where D0 is calculated at t=0 from initital configutation
  3- FORCE_LE_INITIAL_DISP_MULTI:       Linear Elastic Force = K*(D-D0) where D0 is calculated at t=0 from initital configutation &
                                          K varies at each spring according to a dat file 
                                  */

// Forces: Soft Force Models
typedef enum {FORCE_SF_NULL=0, FORCE_SF_ELLIPSOID=1, FORCE_SF_TWO_ELLIPSOIDS=2, FORCE_SF_ELLIPSOID_REMOVE_NORMAL_FORCE=3, FORCE_SF_PLANAR=4, FORCE_SF_ELLIPSOID_OUTSIDE=5} FORCE_SF_TYPE;
const char* FORCE_SF_CHAR[] = {"NO_SOFT_FORCE", "FORCE_SF_ELLIPSOID", "FORCE_SF_TWO_ELLIPSOIDS" , "FORCE_SF_ELLIPSOID_REMOVE_NORMAL_FORCE","FORCE_SF_PLANAR", "FORCE_SF_ELLIPSOID_OUTSIDE"};
/* Info about the SOLID_FORCE_SF_TYPE:
  0- FORCE_SF_NULL:                           
  1- FORCE_SF_ELLIPSOID:                      Add a soft force to constraint the point on one ellipsoid
  2- FORCE_SF_TWO_ELLIPSOIDS:                 Remove normal force components of the chosen ellipsoid
  3- FORCE_SF_ELLIPSOID_REMOVE_NORMAL_FORCE:  Add a soft force on two ellipsoids to remove normal component of force
  4- FORCE_SF_PLANAR:                         Add a soft force on points on a plane to keep them on the plane
  5- FORCE_SF_ELLIPSOID_OUTSIDE:              Add a soft force to constraint the point on not leave an ellipsoid (they can go in but not outside)
  */

// Forces: Surface Tension Models
typedef enum {FORCE_ST_NULL=0, FORCE_ST_CONST=1, FORCE_ST_MULTI=2, FORCE_ST_CONST_PERC=3, 
  FORCE_ST_MULTI_PERC=4, FORCE_ST_CELLS=5} FORCE_ST_TYPE;
const char* FORCE_ST_CHAR[] = {"NO_SURFACE_TENSION_FORCE", "SURFACE_TENSION_CONST", "SURFACE_TENSION_MULT", "SURFACE_TENSION_PERCENTAGE", "SURFACE_TENSION_MULTI_PERCENTAGE", "SURFACE_TENSION_CELLS"};
/* Info about the SOLID_FORCE_ST_TYPE:
  0- FORCE_ST_NULL:         No Surface Tension
  1- FORCE_ST_CONST:        Const Surface Tension                       (F_st = gamma*nT)
  2- FORCE_ST_MULTI:        Multiple Surface Tension for each edge      (F_st = gamma*nT*M)
  3- FORCE_ST_CONST_PERC:   Const Surface Tension percetange from K     (F_st = gamma_perc*K*D) 
  4- FORCE_ST_MULTI_PERC:   Multiple Surface Tension percetange from K for each edge 
  5- FORCE_ST_CELLS:        Surface Tension from the global change of the area of a cell [Dupin (Modeling the flow of dense suspensions of deformable particles in three dimensions) Eq. 15]
  */

// Forces: Volume Conservation Models
typedef enum {FORCE_VC_NULL=0, FORCE_VC_CONST=1, FORCE_VC_CONST_CYLIN=2, FORCE_VC_CONST_CELLS=3, FORCE_VC_CONST_CYLIN_CELLS=4} FORCE_VC_TYPE;
const char* FORCE_VC_CHAR[] = {"NO_VOLUME_CONSERVATION_FORCE", "VOLUME_CONSERVATION_CONST", "VOLUME_CONSERVATION_CONST_CYLINDER", "VOLUME_CONSERVATION_CONST_CELLS",
"VOLUME_CONSERVATION_CONST_CYLINDER_AND_CELLS"};
/* Info about the SOLID_FORCE_VC_TYPE:
  0- FORCE_VC_NULL:               No Volume Conservation
  1- FORCE_VC_CONST:              Const Volume Conservation force is added 
  2- FORCE_VC_CONST_CYLIN:        Const Volume Conservation force is added but the method of volume calculation is for cylinders:
                                    It is based upon the total area of triangles and the area from specific lines that make up the circles
  3- FORCE_VC_CELLS:              Const Volume Conservation force is added to each cell
  4- FORCE_VC_CONST_CYLIN_CELLS:  Both FORCE_VC_CONST_CYLIN + FORCE_VC_CELLS are applied 
  */  

// Solid: Constrain
typedef enum {CONSTRAIN_NULL=0, CONSTRAIN_ELLIPSOID=1, CONSTRAIN_TWO_ELLIPSOIDS=2, CONSTRAIN_AXIS=3, CONSTRAIN_TWO_AXES=4, CONSTRAIN_THREE_AXES=5} CONSTRAIN_TYPE;
const char* CONSTRAIN_CHAR[] = {"NO_CONSTRAINT", "CONSTRAIN_ELLIPSOID", "CONSTRAIN_TWO_ELLIPSOIDS", "CONSTRAIN_AXIS", "CONSTRAIN_TWO_AXES", "CONSTRAIN_THREE_AXES"};
/* Info about the Constrain:
  0- CONSTRAIN_NULL:                      No Geometric constrain on the model
  1- CONSTRAIN_ELLIPSOID:                 Constrain the point on one ellipsoid
  2- CONSTRAIN_TWO_ELLIPSOIDS:            Constrain the point on two ellipsoids [An ellipsoid (A, Bd, C) for y>0 and another ellipsoid (A, Bv, C) for y<=0]
  3- CONSTRAIN_AXIS:                      Constrain on one axis
  4- CONSTRAIN_TWO_AXES:                  Constrain on two axes
  5- CONSTRAIN_THREE_AXES:                Constrain on three axes
  */

// Dynamic Mesh Models
typedef enum {DYNAMIC_MESH_NULL=0, DYNAMIC_MESH_ELASTIC_LENGTH=1} DYNAMIC_MESH_TYPE;
const char* DYNAMIC_MESH_CHAR[] = {"NO_DYNAMIC_MESHING", "DYNAMIC_MESH_ELASTIC_LENGTH"};
/* Info about the SOLID_FORCE_SF_TYPE:
  0- DYNAMIC_MESH_NULL:                           
  1- DYNAMIC_MESH_ELASTIC_LENGTH:    Divide lines based on the length of L 
                                    (i.e. if L > L_divid ==> Divide an edge into two edges)
                                          or L < L_merge ==> Merge edges to become one)
  */

