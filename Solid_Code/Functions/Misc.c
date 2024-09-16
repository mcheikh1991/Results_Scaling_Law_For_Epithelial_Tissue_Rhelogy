#include <errno.h>

int ReadFileLine(FILE *fp,size_t len,char string[])
{
  char *ptr = fgets(string, len, fp);
  if (!ptr) {
    string[0] = 0;
    if (!feof(fp)) { 
      printf("    Error reading from file %d", errno);
    exit(EXIT_FAILURE);}
  }
  return(0);
}

int StrToArray(const char s[], char sp, int *argc, char ***args)
{
  /*@C
      StrToArray - Separates a string by a character (for example ' ' or '\n') and creates an array of strings

      Not Collective

      Input Parameters:
   +  s - pointer to string
   -  sp - separator character

      Output Parameters:
   +   argc - the number of entries in the array
   -   args - an array of the entries with a null at the end
  @*/
  int  i,j,n,*lens,cnt = 0;
  bool flg = false;

  if (!s) n = 0;
  else    n = strlen(s);
  *argc = 0;
  *args = NULL;
  for (; n>0; n--) {   /* remove separator chars at the end - and will empty the string if all chars are separator chars */
   if (s[n-1] != sp) break;
  }
  if (!n) return 0;
  for (i=0; i<n; i++) {
   if (s[i] != sp) break;
  }
  for (;i<n+1; i++) {
   if ((s[i] == sp || s[i] == 0) && !flg) {flg = true; (*argc)++;}
   else if (s[i] != sp) {flg = false;}
  }

  (*args) = (char**) malloc(((*argc)+1)*sizeof(char*)); if (!*args) {printf("    Error in StrToArray\n"); exit(EXIT_FAILURE);}
  lens    = (int*) malloc((*argc)*sizeof(int)); if (!lens) {printf("    Error in StrToArray\n"); exit(EXIT_FAILURE);}
  for (i=0; i<*argc; i++) lens[i] = 0;

  *argc = 0;
  for (i=0; i<n; i++) {
   if (s[i] != sp) break;
  }
  for (;i<n+1; i++) {
   if ((s[i] == sp || s[i] == 0) && !flg) {flg = true; (*argc)++;}
  else if (s[i] != sp) {lens[*argc]++;flg = false;}
  }

  for (i=0; i<*argc; i++) {
    (*args)[i] = (char*) malloc((lens[i]+1)*sizeof(char));
    if (!(*args)[i]) {
      free(lens);
      for (j=0; j<i; j++) free((*args)[j]);
      free(*args);
      printf("    Error in StrToArray\n"); exit(EXIT_FAILURE);
    }
  }
  free(lens);
  (*args)[*argc] = NULL;

  *argc = 0;
  for (i=0; i<n; i++) {
    if (s[i] != sp) break;
  }
  for (;i<n+1; i++) {
    if ((s[i] == sp || s[i] == 0) && !flg) {flg = true; (*args)[*argc][cnt++] = 0; (*argc)++; cnt = 0;}
   else if (s[i] != sp && s[i] != 0) {(*args)[*argc][cnt++] = s[i]; flg = false;}
  }

  // remove the new line
  int len = strlen((*args)[*argc-1]);
  if( (*args)[*argc-1][len-1] == '\n' || (*args)[*argc-1][len-1] == '\r' )
      (*args)[*argc-1][len-1] = 0;

  return 0;
}

int StrToArrayDestroy(int argc, char **args)
{

 /*@C
    StrToArrayDestroy - Frees array created with PetscStrToArray().

    Not Collective

     Output Parameters:
  +  argc - the number of arguments
  -  args - the array of arguments
 */
  
  for (int i = 0; i < argc; ++i) free(args[i]);
  if (args) free(args);
  return 0;
}


// A utility function to swap two elements
void swap(int* a, int* b)
{
  int t = *a;
  *a = *b;
  *b = t;
}

void swapDouble(double* a, double* b)
{
  double t = *a;
  *a = *b;
  *b = t;
}
 
/* This function takes last element as pivot, places the pivot element at its correct position in sorted
array, and places all smaller (smaller than pivot) to left of pivot and all greater elements to right
of pivot */
int partition (int arr[], int low, int high)
{
  int pivot = arr[high]; // pivot
  int i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far

  for (int j = low; j <= high - 1; j++)
  {
      // If current element is smaller than the pivot
      if (arr[j] < pivot)
      {
          i++; // increment index of smaller element
          swap(&arr[i], &arr[j]);
      }
  }
  swap(&arr[i + 1], &arr[high]);
  return (i + 1);
}

/* The main function that implements QuickSort arr[] --> Array to be sorted, low --> Starting index,
high --> Ending index */
void quickSort(int arr[], int low, int high)
{
  if (low < high)
  {
      /* pi is partitioning index, arr[p] is now
      at right place */
      int pi = partition(arr, low, high);

      // Separately sort elements before
      // partition and after partition
      quickSort(arr, low, pi - 1);
      quickSort(arr, pi + 1, high);
  }
}


int partitionDouble (double arr[], int low, int high)
{
  /* Same as partition but for Double */
  double pivot = arr[high]; // pivot
  int i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far

  for (int j = low; j <= high - 1; j++)
  {
      // If current element is smaller than the pivot
      if (arr[j] < pivot)
      {
          i++; // increment index of smaller element
          swapDouble(&arr[i], &arr[j]);
      }
  }
  swapDouble(&arr[i + 1], &arr[high]);
  return (i + 1);
}
 

void quickSortDouble(double arr[], int low, int high)
{
  /* Same as quickSort but for Double */
  if (low < high)
  {
      /* pi is partitioning index, arr[p] is now
      at right place */
      int pi = partitionDouble(arr, low, high);

      // Separately sort elements before
      // partition and after partition
      quickSortDouble(arr, low, pi - 1);
      quickSortDouble(arr, pi + 1, high);
  }
}

