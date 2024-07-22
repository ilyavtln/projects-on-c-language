#include "pch.h"
#include "source.h"

int main()
{
    setlocale(LC_ALL, "");

    //Задание 1
    int n1 = 1000;
    double* x = NULL;
    double* y = NULL;
    if (x == NULL) x = new double[n1]; 
    if (y == NULL) y = new double[n1]; 
    //ScalarMult( x, y, n1);
    //ScalarMultParallel(x, y, n1, 2);
    if (x) { delete[] x; x = NULL; }
    if (y) { delete[] y; y = NULL; }

    //Задание 2
    int n2 = 2500;
    double* A = NULL;
    double* B = NULL;
    double* C = NULL;
    if (A == NULL) A = new double[n2 * n2];
    if (B == NULL) B = new double[n2 * n2];
    if (C == NULL) C = new double[n2 * n2];
    //MatrixMult(A, B, C, n2);
    //MatrixMultParallel(A, B, C, n2, 2);
    if (A) { delete[] A; A = NULL; }
    if (B) { delete[] B; B = NULL; }
    if (C) { delete[] C; C = NULL; }

    //Задание 3
    int n3 = 5000;
    double* U = NULL;
    double* b = NULL;
    if (A == NULL) A = new double[n3 * n3];
    if (x == NULL) x = new double[n3];
    if (y == NULL) y = new double[n3];
    if (U == NULL) U = new double[n3 * n3];
    if (b == NULL) b = new double[n3];
    SlayFind(A, U, x, y, b, n3);
    //SlayFindParallel(A, U, x, y, b, n3, 16);
    if (A) { delete[] A; A = NULL; }
    if (U) { delete[] U; U = NULL; }
    if (y) { delete[] y; y = NULL; }
    if (b) { delete[] b; b = NULL; }
    if (x) { delete[] x; x = NULL; }

    //Задание 6
    int n6 = 3;
    if (A == NULL) A = new double[n6 * n6];
    if (x == NULL) x = new double[n6];
    if (y == NULL) y = new double[n6];
    //MatrixVectorMult(A, x, y, n6);
    //MatrixVectorMultParallel(A, x, y, n6, );
    if (A) { delete[] A; A = NULL; }
    if (x) { delete[] x; x = NULL; }
    if (y) { delete[] y; y = NULL; }
}
