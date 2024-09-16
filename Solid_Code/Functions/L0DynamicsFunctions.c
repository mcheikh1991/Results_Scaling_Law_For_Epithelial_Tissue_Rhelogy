
void DefaultL0DynamicParameters()
{
  LD.Type = L0_DYNAMICS_NULL;
  LD.lambda = 0.2;
  LD.tau = 2000.0;
}

bool ReadL0DynamicsParameters(const char* str_argv1, const char* str_argv2)
{
  if (strcmp(str_argv1,"-ld_type") == 0){    
    LD.Type = (L0_DYNAMICS_TYPE) atol(str_argv2);  
    return true;}

  if (strcmp(str_argv1,"-ld_lambda") == 0){    
    LD.lambda = atof(str_argv2); 
    return true;}

  if (strcmp(str_argv1,"-ld_tau") == 0){    
    LD.tau = atof(str_argv2);  
    return true;}

  if (strcmp(str_argv1,"-ld_lambda_multiplier_file_loc") == 0) { 
    strcpy(LD.M_lambda_file_name, str_argv2); 
    return true;}

  return false;
}


void PrintL0DynamicsProperties()
{
  printf("L0-Dynamics Properties:\n");
  printf("-----------------------------------\n");
  printf("  LD.Type:\t\t%s\n", L0_DYNAMICS_CHAR[LD.Type]);
  if (LD.Type != L0_DYNAMICS_NULL)
  {
    printf("  LD.tau:\t\t%lf\n", LD.tau);
    if (LD.Type == L0_DYNAMICS_COND)
      printf("  LD.lambda:\t\t%lf\n", LD.lambda);
  }
  printf("\n");
}


void ReadL0DynamicsMeshData()
{
  char fileToRead[240];
  if (LD.Type == L0_DYNAMICS_COND_MULTI_LAMBDA)
  {
    printf("  Reading the LD-lambda-M for each line\n");
    LD.M_lambda = (double *)malloc(lines.N*sizeof(double));

    sprintf(fileToRead, "./%s/ld-lambda-m.dat", LD.M_lambda_file_name);
    ReadDoubleFile(fileToRead, lines.N, LD.M_lambda);
  }
}

//==============================================================================

void L0DynamicsAlways()
{
  #pragma omp parallel
  {
    int l, l_start, l_end, omp_rank = (int) omp_get_thread_num();
    double L, L0;
    const double tau = LD.tau;

    GetLimitsOMPLoop(lines.N, omp_nthreads, omp_rank, &l_start, &l_end);
    for(l=l_start; l<l_end; l++)
    {
      L = CalcLineLength(l);
      L0 = lines.len[l];
      lines.len[l] = (times.dt/tau)*(L-L0) + L0;
    }
  }
}

void L0DynamicsConditional()
{
  #pragma omp parallel
  {
    int l, l_start, l_end, omp_rank = (int) omp_get_thread_num();
    double L, L0;
    const double lambda = LD.lambda;
    const double tau = LD.tau;

    GetLimitsOMPLoop(lines.N, omp_nthreads, omp_rank, &l_start, &l_end);
    for(l=l_start; l<l_end; l++)
    {
      L = CalcLineLength(l);
      L0 = lines.len[l];
      if( L > (1+lambda)*L0)      
        LD.Bool_LNode[l] = true;

      if (LD.Bool_LNode[l])   
        lines.len[l] = (times.dt/tau)*(L-L0) + L0;
    }
    
  }
}


void L0DynamicsConditionalMultiLambda()
{
  #pragma omp parallel
  {
    int l, l_start, l_end, omp_rank = (int) omp_get_thread_num();
    double L, L0;
    const double lambda = LD.lambda;
    const double tau = LD.tau;

    GetLimitsOMPLoop(lines.N, omp_nthreads, omp_rank, &l_start, &l_end);
    for(l=l_start; l<l_end; l++)
    {
      const double M = LD.M_lambda[l]; // lambda multiplier
      L = CalcLineLength(l);
      L0 = lines.len[l];
      if( L > (1+lambda*M)*L0)      
        LD.Bool_LNode[l] = true;

      if (LD.Bool_LNode[l])   
        lines.len[l] = (times.dt/tau)*(L-L0) + L0;
    }
    
  }
}



void L0DynamicsNull()
{}

void ChooseL0DynamicsFunction()
{
  printf("Choosing the correct L0-Dynamic function:\n"); 
  if (LD.Type == L0_DYNAMICS_ALWAYS){
    printf("  L0-Dynamic is 'Always' on, with tau = %lf\n", LD.tau); 
    LD.Func=L0DynamicsAlways;}
  else if (LD.Type == L0_DYNAMICS_COND){
    printf("  L0-Dynamic is 'Conditional' for L>L0*%lf, with tau = %lf\n", (1+LD.lambda), LD.tau); 
    LD.Bool_LNode = (bool *)calloc(lines.N,sizeof(bool));
    LD.Func=L0DynamicsConditional;}
  else if (LD.Type == L0_DYNAMICS_COND_MULTI_LAMBDA){
    printf("  L0-Dynamic is 'Conditional' for L>L0*(1+lambda) with multiple lambdas and tau = %lf\n", LD.tau); 
    LD.Bool_LNode = (bool *)calloc(lines.N,sizeof(bool));
    LD.Func=L0DynamicsConditionalMultiLambda;}
  else {
    printf("  No L0-Dynamic function chosen\n"); 
    LD.Func=L0DynamicsNull;}

}

void ApplyL0Dynamics(){ LD.Func();}