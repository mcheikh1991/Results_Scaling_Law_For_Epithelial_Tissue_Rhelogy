

void GetLimitsOMPLoop(const int N, const int OMP_N_threads, const int OMP_Rank, int *i_start, int *i_end){
  int di = (int) N/OMP_N_threads;
  *i_start = OMP_Rank*di;
  *i_end = (OMP_Rank+1)*di;
  if (OMP_N_threads-1 == OMP_Rank)
    *i_end = N;
}
