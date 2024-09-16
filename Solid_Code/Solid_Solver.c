#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
//#include <sys/resource.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>
#include <ctype.h>
#include <unistd.h>

// Global Variable declarations:
//------------------------------
#include "HeaderFiles/Global_Variables.h"
#include "HeaderFiles/Global_Enum.h"
#include "HeaderFiles/Structs.h"

FILE *fd; char line[256]="";

int DIM;
int omp_nthreads;

Time_Struct               times;        // All times
Points_Struct             points;       // Coordinate of points
Points_Struct             forces;       // Forces on each point
Lines_Struct              lines;        // Edges
Triangles_Struct          triangles;    
Squares_Struct            squares; 
Boundary_Groups_Struct    boundaries;
Boxes_Struct              boxes;
Drag_Struct               drag;

char meshInputFolder[200] = "meshfiles/";
char csvOutputFile[256]   = "out.csv";
char reloadFolder[200];
char vtkOutputFolder[200];
char vtkInputFile[256] = "";
char reloadVTKFile[256] = "";

// Material Properties:
//------------------------------
Forces_LE_Struct        LE; // Linear Elasticity Parameters
Forces_BF_Struct        BF; // Barrier Parameters
Forces_SF_Struct        SF; // Soft Force Parameters
Forces_ST_Struct        ST; // Surface Tension Parameters
Forces_VC_Struct        VC; // Volume Conservation Parameters
L0Dynamic_Struct        LD; // L0 Dynamics
Constrain_Struct        CN; // Constrain Parameters
Dynamic_Mesh_Struct     DM; // Dynamics Mesh

// Functions:
//------------------------------
#include "Functions/Misc.c"
#include "Functions/TimeFunctions.c"
#include "Functions/MathFunctions.c"
#include "Functions/NoiseFunctions.c"
#include "Functions/OMPFunctions.c"
#include "Functions/InputOutputFunctions.c"
#include "Functions/BoundaryFunctions.c"
#include "Functions/BoxesFunctions.c"
#include "Functions/ForcesFunctions.c"
#include "Functions/MoveFunctions.c"
#include "Functions/L0DynamicsFunctions.c"
#include "Functions/ConstrainFunctions.c"
#include "Functions/DynamicMeshFunctions.c"
#include "Functions/Misc2.c"
#include "Functions/ReloadData.c"


int main(int argc, char *argv[]) 
{
  printf("\n");
  printf("====================================\n");
  printf("Solid Solver v8.31\n");
  printf("By: Mohamad Ibrahim Cheikh\n");
  printf("====================================\n");

  // Default Parameters:
  //=====================================
  DIM = 2;
  omp_nthreads = 1;
  if(!getcwd(vtkOutputFolder, sizeof(vtkOutputFolder))) exit(EXIT_FAILURE);

  drag.c = 1.0;                           //  Drag Default Parameter
  DefaultTimeParameters();                //  Time Default Parameters
  DefaultBoundaryParameters();            //  Boundary Default Parameters
  DefaultBoxParameters();                 //  Boxes Default Parameters
  DefaultLinearElasticParameters();       //  Linear Elastic Default Parameters
  DefaultBarrierForceParameters();        //  Barrier Force Default Parameters
  DefaultVolumeConservationParameters();  //  Volume Conservation Default Parameters
  DefaultL0DynamicParameters();           //  L0 Dynamics Default Parameters
  DefaultSoftForceParameters();           //  Soft Force Parameters
  DefaultSurfaceTensionParameters();      //  Surface Tension Parameters
  DefaultConstrainParameters();           //  Constrain Parameters
  DefaultDynamicMeshParameters();         //  Dynamic Mesh Parameters

  // Input Parameters
  //=====================================
  #include "HeaderFiles/ReadInputData.h"

  // Reading Mesh
  //=====================================
  #include "HeaderFiles/ReadMeshFiles.h"

  // Defining Boundary
  //=====================================
  ChooseBoundaryGroupFunc();

  // Creating Boxes
  //=====================================
  ChooseBoxFunctions();

  // Choosing Force Functions
  //=====================================
  ChooseLinearElasticFunc();
  ChooseBarrierForcesFunc();
  ChooseSoftForcesFunc(print_sf_points);
  ChooseL0DynamicsFunction();
  ChooseSurfaceTensionFunc();
  ChooseVolumeConservationFunc();

  // Choosing Constrain Functions
  //=====================================
  ChooseConstrainFuncs();
  ConstrainCreate();

  // Choosing Dynamic Mesh Function
  //=====================================
  ChooseDynamicMeshFunc();

  // Intialize/Reload Data
  //=====================================
  InitailizeOrReloadData();

  // Enter The time step loop
  //=====================================
  printf("\nEntering the time step loop:\n");  times.tstart = myclock(); // Global Clock
  for(; times.real_time<=times.end_time; times.real_time+=times.dt)
  {
    ReportData(); 

    BoxesOraganize();

    Zero3ArrayDouble(forces.N, forces.x, forces.y, forces.z);
    
    ForcesAddLinearElasticForce();
    
    ForcesAddBarrierForces();
    
    ForcesAddSoftForces();
    
    ForcesAddSurfaceTensionForces();
    
    ForcesAddVolumeConservationForce();
    
    ForcesAddBoundary(); 
    
    ConstrainForces();
    
    MovePoints();
    
    MoveBoundary();
    
    ApplyL0Dynamics();
    
    ConstrainPoints();
    
    DynamicMesh();
  
    times.time_step++;
  }

  printf("Finished Open-MP Simulation of Threads:%d\n", omp_nthreads); 

  return 0;
}