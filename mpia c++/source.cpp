#include "pch.h"

void ScalarMult(double x[], double y[], int n)
{
    ofstream f("ScalarMult.txt");
    double skal = 0.0;
    for (int i = 0; i < n; i++)
    {
        x[i] = i % 11;
        y[i] = i % 11;
    }
    cout << endl;
    auto start = steady_clock::now();
    for (int i = 0; i < n; i++)
    {
        skal += x[i] * y[i];
    }
    auto end = steady_clock::now();
    cout << duration<double, milli>(end - start).count() << " ms" << endl;
    f << skal << endl;
    f.close();
}

void ScalarMultParallel(double x[], double y[], int n, int  nthreads)
{
    ofstream f("ScalarMultParallel.txt");
    double skal = 0;
    for (int i = 0; i < n; i++)
    {
        x[i] = i % 11;
        y[i] = i % 11;
    }
    auto start = steady_clock::now();
    #pragma omp parallel num_threads(nthreads)
    {
        #pragma omp for reduction(+:skal)
        for (int i = 0; i < n; i++)
        {
            skal += x[i] * y[i];
        }
    }
    auto end = steady_clock::now();
    cout << duration<double, milli>(end - start).count() << " ms" << endl;
    f << skal << endl;
    f.close();
}

void MatrixMult(double A[], double B[], double C[], int n)
{
    ofstream f("MatrixMult.txt");
    double norm = 0.0;
    for (int i = 0; i < n * n; i++)
    {
        A[i] = i % 11;
        B[i] = i % 11;
    }
    auto start = steady_clock::now();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            C[i * n + j] = 0.0;
            for (int k = 0; k < n; k++)
            {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
    auto end = steady_clock::now();
    cout << duration<double, milli>(end - start).count() << " ms" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            f << C[i * n + j] << " ";
            norm += pow(C[i * n + j], 2);
        }
        f << endl;
    }
    norm = sqrt(norm);
    cout << "Норма матрицы" << " " << setprecision(15) << norm << endl;
    f.close();
}

void MatrixMultParallel(double A[], double B[], double C[], int n, int  nthreads)
{
    ofstream f("MatrixMultParallel.txt");
    double norm = 0.0;
    int k = 0;
    int j = 0;
    for (int i = 0; i < n * n; i++)
    {
        A[i] = i % 11;
        B[i] = i % 11;
    }
    auto start = steady_clock::now();
    #pragma omp parallel num_threads(nthreads)
    {
        #pragma omp for private(j, k)
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                C[i * n + j] = 0.0;
                for (int k = 0; k < n; k++)
                {
                    C[i * n + j] += A[i * n + k] * B[k * n + j];
                }
            }
        }
    }
    auto end = steady_clock::now();
    cout << duration<double, milli>(end - start).count() << " ms" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            f << C[i * n + j] << " ";
            norm += pow(C[i * n + j], 2);
        }
        f << endl;
    }
    norm = sqrt(norm);
    cout << "Норма матрицы" << " " <<  setprecision(15) << norm << endl;
    f.close();
}

double VectorNorm(double y[], int n)
{
    double norm = 0.0;
    for (int i = 0; i < n; i++)
    {
        norm += pow(y[i], 2);
    }
    return sqrt(norm);
}

void SlayFind(double A[], double U[], double x[], double y[], double b[], int n)
{
    ofstream fx("SlayFindX.txt");
    ofstream fy("SlayFindY.txt");
    double norm = 0.0;
    for (int i = 0; i < n * n; i++)
    {
        U[i] = 0;
    }
    for (int i = 0, k = n; i < n; i++, k--)
    {
        U[i + n * i] = 10;
        x[i] = i % 11;
        for (int j = i + n * i + 1; j < i + n * i + k; j++)
        {
            U[j] = 1;
        }
    }
    for (int i = 0; i < n; i++)
    {
        b[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            b[i] += U[i * n + j] * x[j];
        }
    }
    auto start = steady_clock::now();
    for (int i = n - 1; i >= 0; i--)
    {
        double s = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            s += U[i * n + j] * y[j];
        }
        y[i] = (b[i] - s) / (U[i * n + i]);
    }
    auto end = steady_clock::now();
    cout << duration<double, milli>(end - start).count() << " ms" << endl;
    for (int i = 0; i < n; i++)
    {
        fx << x[i] << " ";
    }
    for (int i = 0; i < n; i++)
    {
        fy << y[i] << " ";
    }
    fy.close();
    cout << "Норма вектора x: " << setprecision(15) << VectorNorm(x, n) << endl;
    cout << "Норма вектора y: " << setprecision(15) << VectorNorm(y, n) << endl;
    fx.close();
    fy.close();
}

void SlayFindParallel(double A[], double U[], double x[], double y[], double b[], int n, int  nthreads)
{
    ofstream fx("SlayFindParallelX.txt");
    ofstream fy("SlayFindParallelY.txt");
    double norm = 0.0;
    int j = 0;
    for (int i = 0; i < n * n; i++)
    {
        U[i] = 0;
    }
    for (int i = 0, k = n; i < n; i++, k--)
    {
        U[i + n * i] = 10;
        x[i] = i % 11;
        for (int j = i + n * i + 1; j < i + n * i + k; j++)
        {
            U[j] = 1;
        }
    }
    for (int i = 0; i < n; i++)
    {
        b[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            b[i] += U[i * n + j] * x[j];
        }
    }
    auto start = steady_clock::now();
    for (int i = n - 1; i >= 0; i--)
    {
        y[i] = 0.0;
        double s = 0.0;
        #pragma omp parallel for reduction(+:s) num_threads(nthreads)
        for (int j = i + 1; j < n; j++)
        {
            s += U[i * n + j] * y[j];
        }
        y[i] = (b[i] - s) / (U[i * n + i]);
    }
    auto end = steady_clock::now();
    cout << duration<double, milli>(end - start).count() << " ms" << endl;
    for (int i = 0; i < n; i++)
    {
        fx << x[i] << endl;
    }
    for (int i = 0; i < n; i++)
    {
        fy << y[i] << endl;
    }
    fy.close();
    cout << "Норма вектора x: " << setprecision(15) << VectorNorm(x, n) << endl;
    cout << "Норма вектора y: " << setprecision(15) << VectorNorm(y, n) << endl;
    fx.close();
    fy.close();
}

void MatrixVectorMult(double A[], double x[], double y[], int n)
{
    ofstream f("MatrixVectorMult.txt");
    double norm = 0.0;
    for (int i = 0; i < n; i++)
    {
        x[i] = 1;
    }
    for (int i = 0; i < n * n; i++)
    {
        A[i] = i % 11;
    }
    auto start = steady_clock::now();
    for (int i = 0; i < n; i++)
    {
        double temp = 0;
        y[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            temp += A[i * n + j] * x[j];
        }
        y[i] = temp;
    }
    auto end = steady_clock::now();
    cout << duration<double, milli>(end - start).count() << " ms" << endl;
    for (int i = 0; i < n; i++)
    {
        f << y[i] << endl;
    }
    norm = VectorNorm(y, n);
    cout << setprecision(15) << norm << endl;
    f.close();
}

void MatrixVectorMultParallel(double A[], double x[], double y[], int n, int  nthreads)
{
    ofstream f("MatrixVectorMultParallel.txt");
    double norm = 0.0;
    for (int i = 0; i < n; i++)
    {
        x[i] = 1;
    }
    for (int i = 0; i < n * n; i++)
    {
        A[i] = i % 11;
    }
    auto start = steady_clock::now();
    #pragma omp parallel num_threads(nthreads)
    {
        #pragma omp for
        for (int i = 0; i < n; i++)
        {
            double temp = 0;
            y[i] = 0.0;
            for (int j = 0; j < n; j++)
            {
                temp += A[i * n + j] * x[j];
            }
            y[i] = temp;
        }
    }
    auto end = steady_clock::now();
    cout << duration<double, milli>(end - start).count() << " ms" << endl;
    for (int i = 0; i < n; i++)
    {
        f << y[i] << " ";
    }
    norm = VectorNorm(y, n);
    cout << setprecision(15) << norm << endl;
    f.close();
}
