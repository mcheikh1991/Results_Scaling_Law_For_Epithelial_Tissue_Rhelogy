

void DefaultBoundaryParameters()
{  
  boundaries.N = 0;        // Default no boundaries
}

// NO_SOLID_BC = 0
//=============================================================

void BoundaryPointsNull(void *v){}

// FORCE_ADD = 1
//=============================================================

void BoundaryPointsAddForce(void *v)
{
  Boundary_Points_Struct* _bp = (Boundary_Points_Struct *)v;
  for(int i=0; i<_bp->N; i++)
  {
    const int p  = _bp->I[i];
    forces.x[p] += _bp->Fx;
    forces.y[p] += _bp->Fy;
    forces.z[p] += _bp->Fz;
  }
}


// DIST_FIXED = 3
//=============================================================
void InitializePointsDistFixed(Boundary_Points_Struct *_bp)
{
  printf("  Total number of boundary points for DIST_FIXED boundary '%s' is %d\n", _bp->Name ,_bp->N);

  _bp->x = (double *)malloc(_bp->N*sizeof(double));
  _bp->y = (double *)malloc(_bp->N*sizeof(double));
  _bp->z = (double *)malloc(_bp->N*sizeof(double));

  // Save original position
  for(int i=0; i<_bp->N; i++)
  {
    const int p = _bp->I[i];
    _bp->x[i] = points.x[p]; 
    _bp->y[i] = points.y[p]; 
    _bp->z[i] = points.z[p]; 
  }
}


void BoundaryPointsDistFixed(void *v)
{
  Boundary_Points_Struct* _bp = (Boundary_Points_Struct *)v;
  for(int i=0; i<_bp->N; i++)
  {
    const int p = _bp->I[i];
    points.x[p] = _bp->x[i];
    points.y[p] = _bp->y[i];
    points.z[p] = _bp->z[i];
  }
}


// FORCE_ADD_ELLIPSE = 11
//=============================================================
void InitializePointsAddForceEllipse(Boundary_Points_Struct *_bp)
{
  printf("  Total number of boundary points for Add Force Ellipse boundary labeled '%s' is %d\n", _bp->Name ,_bp->N);
  _bp->z = (double *)malloc(_bp->N*sizeof(double));

  // Save original position
  for(int i=0; i<_bp->N; i++)
  {
    const int p = _bp->I[i];
    _bp->z[i] = points.z[p]; 
  }
}

void BoundaryPointsAddForceEllipse(void *v)
{
  Boundary_Points_Struct* _bp = (Boundary_Points_Struct *)v;
  const double A = _bp->A, B = _bp->B;
  const int Dir  = _bp->Dir;
  const double Mag =  _bp->Fm; // Force Mag

  // Boundary Points
  double x3,y3;
  double nVx, nVy, Vm;
  for (int i=0; i<_bp->N; i++)
  {
    const int p = _bp->I[i]; 
    const double x0 = points.x[p];
    const double y0 = points.y[p];

    // Find intersection p'(x3,y3) of point p(x0,y0) with Ellipse (A,B)
    const double x1 =  A*B/sqrt(A*A*y0*y0 + B*B*x0*x0)*x0;
    const double y1 =  A*B/sqrt(A*A*y0*y0 + B*B*x0*x0)*y0;
    const double x2 = -A*B/sqrt(A*A*y0*y0 + B*B*x0*x0)*x0;
    const double y2 = -A*B/sqrt(A*A*y0*y0 + B*B*x0*x0)*y0;
    if ((x0>0.0) && (x1>0.0)){ x3=x1; }
    else                     { x3=x2; }
    if ((y0>0.0) && (y1>0.0)){ y3=y1; }
    else                     { y3=y2; }
    const double Theta3 = CalcEllipsoidAngle(x3, y3, A, B);
    const double Theta4 = Theta3 + Dir*1e-3; // Theta of the next point after p'

    // Normal Vector from the Gradient of Ellipse (A,V) 
    const double Nx = 2.0*x3/(A*A);
    const double Ny = 2.0*y3/(B*B);

    const double T1x = -Ny;
    const double T1y = Nx;
    const double ThetaT1x = CalcEllipsoidAngle(x3 + T1x, y3 +T1y, A, B);

    const double T2x = Ny;
    const double T2y = -Nx;
    const double ThetaT2x = CalcEllipsoidAngle(x3 + T2x, y3 +T2y, A, B);

    if (fabs(ThetaT1x - Theta4) < fabs(ThetaT2x - Theta4)) {
      Vm = sqrt(T1x*T1x + T1y*T1y);    nVx = T1x/Vm;       nVy = T1y/Vm;}
    else {
      Vm = sqrt(T2x*T2x + T2y*T2y);    nVx = T2x/Vm;       nVy = T2y/Vm;}

    forces.x[p] += Mag*nVx;   // [N]
    forces.y[p] += Mag*nVy;   // [N]
    forces.z[p] += 0.0;   // [N]
    points.z[p] = _bp->z[i];
    //printf("p(%lf, %lf, %lf) | nV(%lf, %lf, 0.0) | (%lf, %lf) | Theta3=%lf | Theta4=%lf\n",points.x[p],points.y[p],points.z[p], nVx, nVy, x3, y3, Theta3, Theta4);
  }
  //printf("In\n");
}


// FORCE_SPEED_ELLIPSE = 7
//=============================================================

void InitailizeForceEllipseSpeedPoints(Boundary_Points_Struct *_bp)
{
  printf("  Total number of boundary points for Force Speed Ellipse points is %d\n", _bp->N);

  _bp->theta = (double *)malloc(_bp->N*sizeof(double));
  _bp->Z_Plane = (double *)malloc(_bp->N*sizeof(double));
  const double A = _bp->A;
  const double B = _bp->B;
  const double Dtheta0 = _bp->Dtheta0;
  // Save original position
  for(int i=0; i<_bp->N; i++)
  {
    const int p = _bp->I[i];
    _bp->theta[i] = CalcEllipsoidAngle(points.x[p], points.y[p], A, B) + Dtheta0;
    _bp->Z_Plane[i] = points.z[p];
  }

  // Create CVS File for total external force reporting
  strcpy(_bp->ExternalForceFile,vtkOutputFolder);
  strcat(_bp->ExternalForceFile,"/ExternalForces.csv");
  fd = fopen(_bp->ExternalForceFile,"w");
  fprintf(fd,"tims[ms],fx[nN],fy[nN],fz[nN]\n");
  fclose(fd); 

  // Initialize Total External Force
  _bp->Total_External_Force[0] = 0.0;
  _bp->Total_External_Force[1] = 0.0;
  _bp->Total_External_Force[2] = 0.0;
}


void BoundaryPointsAddForceSpeedEllipse(void *v)
{
  Boundary_Points_Struct* _bp = (Boundary_Points_Struct *)v;
  const double A = _bp->A;
  const double B = _bp->B;
  const double K =  _bp->K;             // P controler

  // Initialize Total External Force
  _bp->Total_External_Force[0] = 0.0;
  _bp->Total_External_Force[1] = 0.0;
  _bp->Total_External_Force[2] = 0.0;

  // Boundary Points
  for (int i=0; i<_bp->N; i++)
  {
    const int p = _bp->I[i];
    const double theta = _bp->theta[i]; 
    const double Z_Plane = _bp->Z_Plane[i];

    // Calculate the tangential vector from real p --> imag p
    const double Vx = A*cos(theta) - points.x[p];
    const double Vy = B*sin(theta) - points.y[p];
    const double Vz = Z_Plane - points.z[p];

    // Extran_Force = Sum(Elastic Forces) + Fictious_Force
    _bp->Total_External_Force[0] += -1*forces.x[p] + K*Vx; // D-D0 is in Vx
    _bp->Total_External_Force[1] += -1*forces.y[p] + K*Vy; // D-D0 is in Vy
    _bp->Total_External_Force[2] += -1*forces.z[p] + K*Vz; // D-D0 is in Vz

    // Do not use here += since we are assuming these points only have BC force applies
    forces.x[p] = K*Vx; // [N]
    forces.y[p] = K*Vy; // [N]
    forces.z[p] = K*Vz; // [N]
  }
}

void BoundaryPointsReportTotalExternalForceSpeedEllipse(void *v)
{
  Boundary_Points_Struct* _bp = (Boundary_Points_Struct *)v;
  fd = fopen(_bp->ExternalForceFile,"a");
  fprintf(fd,"%lf,%lf,%lf,%lf\n", times.real_time, _bp->Total_External_Force[0], _bp->Total_External_Force[1], _bp->Total_External_Force[2]);
  fclose(fd); // Close Flie
}

void BoundaryPointsMoveForceSpeedEllipse(void *v)
{
  Boundary_Points_Struct* _bp = (Boundary_Points_Struct *)v;

  const double DthetaDt = _bp->DthetaDt;
  // Move the ghost point for BC "FORCE_SPEED_ELLIPSE"
  for(int i=0; i<_bp->N; i++)
    _bp->theta[i] += times.dt*DthetaDt;
}


//===============================================================================
//==============================================================================

// Compile all the above

void ChooseBoundaryFunc(void *v)
{
  Boundary_Points_Struct* _bp = (Boundary_Points_Struct *)v;
  printf("  Choosing the boundary function for boundary labeled '%s':\n", _bp->Name); 
  if ( _bp->Type == NO_SOLID_BC)
  {
    printf("    No boundary function choosen\n"); 
    _bp->Force_Func  = BoundaryPointsNull;
    _bp->Move_Func   = BoundaryPointsNull;
    _bp->Report_Func = BoundaryPointsNull;
  }
  else if ( _bp->Type == DIST_FIXED)        // Imployed in function 'MoveBoundary'
  {
    printf("    A fixed distance boundary function is choosen for labeled '%s'.\n", _bp->Name); 
    InitializePointsDistFixed(_bp);
    _bp->Force_Func  = BoundaryPointsNull;
    _bp->Move_Func   = BoundaryPointsDistFixed;
    _bp->Report_Func = BoundaryPointsNull;
  }
  else if ( _bp->Type == FORCE_ADD)         // Imployed in function 'ForcesAddBoundary'
  {
    printf("    An additional force will be added to boundary points labeled '%s'.\n", _bp->Name); 
    _bp->Force_Func  = BoundaryPointsAddForce;
    _bp->Move_Func   = BoundaryPointsNull;
    _bp->Report_Func = BoundaryPointsNull;
  }
  else if ( _bp->Type == FORCE_ADD_ELLIPSE) // Imployed in function 'ForcesAddBoundary'
  {
    printf("    An additional force tangent to an ellipse will be added to boundary points labeled '%s'.\n", _bp->Name);
    InitializePointsAddForceEllipse(_bp);
    _bp->Force_Func  = BoundaryPointsAddForceEllipse;
    _bp->Move_Func   = BoundaryPointsNull;
    _bp->Report_Func = BoundaryPointsNull;
  }
  else if ( _bp->Type == FORCE_SPEED_ELLIPSE)
  {
    printf("    An additional force controlloing the speed around an ellipse will be added to boundary points labeled '%s'.\n", _bp->Name);
    InitailizeForceEllipseSpeedPoints(_bp);
    _bp->Force_Func  = BoundaryPointsAddForceSpeedEllipse;
    _bp->Move_Func   = BoundaryPointsMoveForceSpeedEllipse;
    _bp->Report_Func = BoundaryPointsReportTotalExternalForceSpeedEllipse;
  }
  else{
    printf("  ERROR!! Unknown boundary function for type '%s'\n", SOLID_BC_TYPE_CHAR[_bp->Type]);
    exit(EXIT_FAILURE);}
}

void ChooseBoundaryGroupFunc(){
  printf("Number of total boundaries is %d:\n",boundaries.N);
  for(int g=0; g<boundaries.N; g++)
    ChooseBoundaryFunc(&boundaries.groups[g]);
}

void ForcesAddBoundary(){ 
  for(int g=0; g<boundaries.N; g++)
    boundaries.groups[g].Force_Func(&boundaries.groups[g]);
}

void MoveBoundary(){
  for(int g=0; g<boundaries.N; g++)
    boundaries.groups[g].Move_Func(&boundaries.groups[g]);
}

void ReportBoundary(){
  for(int g=0; g<boundaries.N; g++)
    boundaries.groups[g].Report_Func(&boundaries.groups[g]);
}
