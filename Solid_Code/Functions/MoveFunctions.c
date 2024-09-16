void MovePoints()
{  
  for(int p=0; p<points.N; p++)
  {
    points.x[p] += forces.x[p]*times.dt/drag.c;
    points.y[p] += forces.y[p]*times.dt/drag.c;
    points.z[p] += forces.z[p]*times.dt*(DIM-2)/drag.c;
  }
}