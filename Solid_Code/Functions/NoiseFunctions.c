

void AddNoise()
{
  for(int i=0; i<points.N; i++)
  {
    double rnd = ((double) rand())/((double) RAND_MAX);
    points.x[i] += 0.01*(2*rnd-1);
    rnd = ((double) rand())/((double) RAND_MAX);
    points.y[i] += 0.01*(2*rnd-1);
    rnd = ((double) rand())/((double) RAND_MAX);
    points.z[i] += 0.01*(2*rnd-1)*(DIM-2);
  }
}

