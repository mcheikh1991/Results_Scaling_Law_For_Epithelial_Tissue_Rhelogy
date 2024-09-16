void ZeroArrayDouble(const int ArraySize, double *_a){
  for(int i=0; i<ArraySize; i++)
    _a[i] = 0.0;
}

void Zero3ArrayDouble(const int ArraySize, double *_a1, double *_a2, double *_a3){
  for(int i=0; i<ArraySize; i++)
    _a1[i] = 0.0, _a2[i] = 0.0, _a3[i] = 0.0;
}

void ZeroArrayInt(const int ArraySize, int *_a){
  for(int i=0; i<ArraySize; i++)
    _a[i] = 0;
}

double CalcDistance(const int p0, const int p1) 
{
 double Dx = points.x[p1]-points.x[p0];
 double Dy = points.y[p1]-points.y[p0];
 double Dz = points.z[p1]-points.z[p0];
 return sqrt( Dx*Dx + Dy*Dy + Dz*Dz);
}


double CalcDistance2(const double p0x, const double p0y, const double p0z, 
  const double p1x, const double p1y, const double p1z) 
{
  /* Calculates the distance between two points p0 (p0x, p0y, p0z) and p1 (p1x, p1y, p1z)
     x, y, z are a 1D array for all the coordinates */
 double Dx = p1x-p0x;
 double Dy = p1y-p0y;
 double Dz = p1z-p0z;
 return sqrt( Dx*Dx + Dy*Dy + Dz*Dz );
}


double CalcLineLength(const int l){
 return CalcDistance(lines.I0[l], lines.I1[l]); } 

void CalcNormalVector(const int p0, const int p1, double *nX, double *nY, double *nZ, double *N)
{
  /* Calculates the normalized vector from two points p0 and p1*/
  *nX = points.x[p1] - points.x[p0];
  *nY = points.y[p1] - points.y[p0];
  *nZ = points.z[p1] - points.z[p0];
  *N = sqrt( (*nX)*(*nX) + (*nY)*(*nY) + (*nZ)*(*nZ));
  *nX = (*nX) / (*N);
  *nY = (*nY) / (*N);
  *nZ = (*nZ) / (*N);
}

void CalcNormalVector2( const double p0x, const double p0y, const double p0z,
  const double p1x, const double p1y, const double p1z, double *nX, double *nY, double *nZ, double *N)
{
  /* Calculates the normalized vector from two points p0 and p1*/
  *nX = p1x - p0x;
  *nY = p1y - p0y;
  *nZ = p1z - p0z;
  *N = sqrt( (*nX)*(*nX) + (*nY)*(*nY) + (*nZ)*(*nZ));
  *nX = (*nX) / (*N);
  *nY = (*nY) / (*N);
  *nZ = (*nZ) / (*N);
}


void CalcVector(const int p0, const int p1, double *vX, double *vY, double *vZ)
{
  /* Calculates the normalized vector from two points p0 and p1*/
  *vX = points.x[p1] - points.x[p0];
  *vY = points.y[p1] - points.y[p0];
  *vZ = points.z[p1] - points.z[p0];
}

double CalcDotProduct( const double V1x, const double V1y, const double V1z,
                       const double V2x, const double V2y, const double V2z){
  /* Do dot porduct of V1.V2 */
  return V1x*V2x + V1y*V2y + V1z*V2z;
}


void CalcCrossProduct( const double V1x, const double V1y, const double V1z,
                       const double V2x, const double V2y, const double V2z,
                             double *Vx,       double *Vy,       double *Vz){
  /* Do cross porduct of V1xV2 */
  *Vx = V1y*V2z - V1z*V2y;
  *Vy = V1z*V2x - V1x*V2z;
  *Vz = V1x*V2y - V1y*V2x;
}

void CalcCrossProductNormalzied( const double V1x, const double V1y, const double V1z,
                                 const double V2x, const double V2y, const double V2z,
                                       double *nVx,      double *nVy,      double *nVz){
  /* Do cross porduct of V1xV2 */
  CalcCrossProduct(V1x, V1y, V1z, V2x, V2y, V2z, nVx, nVy, nVz);
  double Vm = sqrt( (*nVx)*(*nVx) + (*nVy)*(*nVy) + (*nVz)*(*nVz));
  *nVx = (*nVx) / (Vm);
  *nVy = (*nVy) / (Vm);
  *nVz = (*nVz) / (Vm);  
}

double CalcTriangleArea( const double p0x, const double p0y, const double p0z, 
                         const double p1x, const double p1y, const double p1z,
                         const double p2x, const double p2y, const double p2z)
{
  double Vx, Vy, Vz;
  CalcCrossProduct( p1x-p0x, p1y-p0y, p1z-p0z, 
                    p2x-p0x, p2y-p0y, p2z-p0z,
                        &Vx,     &Vy,     &Vz);
  return 0.5*sqrt( (Vx)*(Vx) + (Vy)*(Vy) + (Vz)*(Vz));
}

double CalcTriangleArea2(const int t)
{
  const int p0 = triangles.I0[t];
  const int p1 = triangles.I1[t];
  const int p2 = triangles.I2[t];
  return CalcTriangleArea(points.x[p0], points.y[p0], points.z[p0],
                          points.x[p1], points.y[p1], points.z[p1],
                          points.x[p2], points.y[p2], points.z[p2]);
}

void CalcTriangleNormal(const int t, double *nVx, double *nVy, double *nVz)
{
  const int p0 = triangles.I0[t];
  const int p1 = triangles.I1[t];
  const int p2 = triangles.I2[t];
  CalcCrossProductNormalzied( points.x[p1]-points.x[p0], points.y[p1]-points.y[p0], points.z[p1]-points.z[p0],
                              points.x[p2]-points.x[p0], points.y[p2]-points.y[p0], points.z[p2]-points.z[p0],
                                                    nVx,                       nVy,                       nVz);
}

double CalcSquareArea(const int s)
{
  const int p0 = squares.I0[s];
  const int p1 = squares.I1[s];
  const int p2 = squares.I2[s];
  const int p3 = squares.I3[s];
  return  CalcTriangleArea( points.x[p0], points.y[p0], points.z[p0],
                            points.x[p1], points.y[p1], points.z[p1],
                            points.x[p2], points.y[p2], points.z[p2]) +
          CalcTriangleArea( points.x[p0], points.y[p0], points.z[p0],
                            points.x[p2], points.y[p2], points.z[p2],
                            points.x[p3], points.y[p3], points.z[p3]);
}

double CalcEllipsoidAngle(const double x, const double y, const double A, const double B)
{
  /* Calcualte the Azimuthal Theta angle on an ellipsoid A, B */
  const double   PI = 3.14159265358979323846; // PI
  double Theta = 0;
  if (x==0) {
    if      (y<0)   Theta = -PI/2;
    else if (y>0)   Theta =  PI/2;
    else if (y==0)  Theta =   0.0; }
  else {
    Theta = atan((y*A)/(x*B));
    if (x<0) {
      if      (y<0)  Theta = Theta - PI;  
      else if (y>0)  Theta = Theta + PI;
      else if (y==0) Theta = PI;  }
  }
  return Theta;
}

void CalcMidpoint(const int l, double *Mx, double *My, double *Mz) 
{
  const int p0 = lines.I0[l];
  const int p1 = lines.I1[l];
  *Mx = (points.x[p0] + points.x[p1])/2.0;
  *My = (points.y[p0] + points.y[p1])/2.0;
  *Mz = (points.z[p0] + points.z[p1])/2.0;
}

void CalcTriangleCenter(const int t, double *Mx, double *My, double *Mz)
{
  const int p0 = triangles.I0[t];
  const int p1 = triangles.I1[t];
  const int p2 = triangles.I2[t];
  *Mx = (points.x[p0] + points.x[p1]+ points.x[p2])/3.0;
  *My = (points.y[p0] + points.y[p1]+ points.y[p2])/3.0;
  *Mz = (points.z[p0] + points.z[p1]+ points.z[p2])/3.0;
}


void CalcSquareCenter(const int s, double *Mx, double *My, double *Mz)
{
  const int p0 = squares.I0[s];
  const int p1 = squares.I1[s];
  const int p2 = squares.I2[s];
  const int p3 = squares.I3[s];
  *Mx = (points.x[p0] + points.x[p1]+ points.x[p2]+ points.x[p3])/4.0;
  *My = (points.y[p0] + points.y[p1]+ points.y[p2]+ points.y[p3])/4.0;
  *Mz = (points.z[p0] + points.z[p1]+ points.z[p2]+ points.z[p3])/4.0;
}


bool checkIfPointIn2DTriangle( const double px, const double py, 
  const double p0x, const double p0y,  const double p1x, const double p1y, const double p2x, const double p2y)
{
  // px,  py  is the point coordinate
  // p0x, p0y are the for the first node of the triangle
  // p1x, p1y are the for the second node of the triangle
  // p2x, p2y are the for the thrid node of the triangle
  // We will solve the sytem of Equations below to find s and t
  // Eq1 = p0x-px + (p1x-p0x)*s + (p2x-p0x)*t
  // Eq2 = p0y-py + (p1y-p0y)*s + (p2y-p0y)*t
  // where:
  //   s =  (-(p0x - p2x)*(p0y - py) + (p0x - px)*(p0y - p2y))/((p0x - p1x)*(p0y - p2y) - (p0x - p2x)*(p0y - p1y));
  //   t =  ( (p0x - p1x)*(p0y - py) - (p0x - px)*(p0y - p1y))/((p0x - p1x)*(p0y - p2y) - (p0x - p2x)*(p0y - p1y));
  // The point p is inside the triangle if 0 <= s <= 1 and 0 <= t <= 1 and s + t <= 1.

  //const PetscScalar tol = 1e-12;

  double s = (-(p0x - p2x)*(p0y - py) + (p0x - px)*(p0y - p2y))/((p0x - p1x)*(p0y - p2y) - (p0x - p2x)*(p0y - p1y));
  double t = ( (p0x - p1x)*(p0y - py) - (p0x - px)*(p0y - p1y))/((p0x - p1x)*(p0y - p2y) - (p0x - p2x)*(p0y - p1y));

  return (((s>=0) && (s<=1) && (t>=0) && (t<=1) && ((s+t)<=1))? true: false);
}

bool checkIfPointIn3DTriangle( const double px,  const double py,  const double pz,
                               const double p0x, const double p0y, const double p0z, 
                               const double p1x, const double p1y, const double p1z,
                               const double p2x, const double p2y, const double p2z  )
{
  // px,  py  is the point coordinate
  // p0x, p0y are the for the first node of the triangle
  // p1x, p1y are the for the second node of the triangle
  // p2x, p2y are the for the thrid node of the triangle

  double PP0x, PP0y, PP0z;
  double PP1x, PP1y, PP1z;
  double PP2x, PP2y, PP2z;
  double Ux, Uy, Uz;
  double Vx, Vy, Vz;
  double Wx, Wy, Wz;

  PP0x = px - p0x;
  PP0y = py - p0y;
  PP0z = pz - p0z;

  PP1x = px - p1x;
  PP1y = py - p1y;
  PP1z = pz - p1z;

  PP2x = px - p2x;
  PP2y = py - p2y;
  PP2z = pz - p2z;

  CalcCrossProduct( PP0x, PP0y, PP0z, PP1x, PP1y, PP1z, &Ux, &Uy, &Uz);
  CalcCrossProduct( PP1x, PP1y, PP1z, PP2x, PP2y, PP2z, &Vx, &Vy, &Vz);
  CalcCrossProduct( PP2x, PP2y, PP2z, PP0x, PP0y, PP0z, &Wx, &Wy, &Wz);

  if ( CalcDotProduct( Ux, Uy, Uz, Vx, Vy, Vz) < 0)
    return false;

  if ( CalcDotProduct( Ux, Uy, Uz, Wx, Wy, Wz) < 0)
    return false;

  return true;
}

bool checkIfPointIn3DSquare( const double px,  const double py,  const double pz,
                             const double p0x, const double p0y, const double p0z, 
                             const double p1x, const double p1y, const double p1z,
                             const double p2x, const double p2y, const double p2z,
                             const double p3x, const double p3y, const double p3z  )
{
  // px,  py  is the point coordinate
  // p0x, p0y are the for the first node of the triangle
  // p1x, p1y are the for the second node of the triangle
  // p2x, p2y are the for the thrid node of the triangle

  double PP0x, PP0y, PP0z;
  double PP1x, PP1y, PP1z;
  double PP2x, PP2y, PP2z;
  double PP3x, PP3y, PP3z;
  double Ux, Uy, Uz;
  double Vx, Vy, Vz;
  double Wx, Wy, Wz;
  double Tx, Ty, Tz;

  PP0x = px - p0x;
  PP0y = py - p0y;
  PP0z = pz - p0z;

  PP1x = px - p1x;
  PP1y = py - p1y;
  PP1z = pz - p1z;

  PP2x = px - p2x;
  PP2y = py - p2y;
  PP2z = pz - p2z;

  PP3x = px - p3x;
  PP3y = py - p3y;
  PP3z = pz - p3z;

  CalcCrossProduct( PP0x, PP0y, PP0z, PP1x, PP1y, PP1z, &Ux, &Uy, &Uz);
  CalcCrossProduct( PP1x, PP1y, PP1z, PP2x, PP2y, PP2z, &Vx, &Vy, &Vz);
  CalcCrossProduct( PP2x, PP2y, PP2z, PP3x, PP3y, PP3z, &Wx, &Wy, &Wz);
  CalcCrossProduct( PP3x, PP3y, PP3z, PP0x, PP0y, PP0z, &Tx, &Ty, &Tz);

  if ( CalcDotProduct( Ux, Uy, Uz, Vx, Vy, Vz) < 0)
    return false;

  if ( CalcDotProduct( Vx, Vy, Vz, Wx, Wy, Wz) < 0)
    return false;

  if ( CalcDotProduct( Wx, Wy, Wz, Tx, Ty, Tz) < 0)
    return false;

  if ( CalcDotProduct( Tx, Ty, Tz, Ux, Uy, Uz) < 0)
    return false;
  
  return true;
}
