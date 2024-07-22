void ScalarMult(double x[], double y[], int n);
void ScalarMultParallel(double x[], double y[], int n, int  nthreads);
void MatrixMult(double A[], double B[], double C[], int n);
void MatrixMultParallel(double A[], double B[], double C[], int n, int  nthreads);
void SlayFind(double A[], double U[], double x[], double y[], double b[], int n);
void SlayFindParallel(double A[], double U[], double x[], double y[], double b[], int n, int  nthreads);
void MatrixVectorMult(double A[], double x[], double y[], int n);
void MatrixVectorMultParallel(double A[], double x[], double y[], int n, int  nthreads);