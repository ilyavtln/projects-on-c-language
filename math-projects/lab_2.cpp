#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>

#define STRIN "%lf"
#define STROUT "%.15lf "
#define EE "%.5e\t"

#pragma warning(disable:4996)

typedef double type;
typedef double summa;

using namespace std;

int n = 0; // размер матрицы
int m = 0; // количество нулевых диагоналей
type w = 0; // параметр релаксации
type e = 0; // точность
int maxiter = 0; // максимальное количество итераций
vector<vector<type>> di;// не нулевые диагонали
vector<type> f;// Вектор
vector<type> xk;
vector<type> x0;// начальное приближение
vector<type> Ax;// для подсчета невязки
vector<type> f_Ax;// для подсчета невязки


void Diagonal(ifstream& matrix, int size, int num) // файл, размер диагонали, номер диагонали
{
    di.resize(n);
    di[num].resize(size);
    for (int j = 0; j < size; j++)
        matrix >> di[num][j];
}

void Input(ifstream& matrix, ifstream& vector, ifstream& x_0)
{
    matrix >> n;
    if (n <= 0) return;

    matrix >> m;
    if (m <= 0) return;

    matrix >> w >> e >> maxiter;
    if (maxiter <= 0) return;

    Diagonal(matrix, n - m - 2, 0);
    Diagonal(matrix, n - m - 1, 1);
    Diagonal(matrix, n - m, 2);
    Diagonal(matrix, n - 1, 3);
    Diagonal(matrix, n, 4);
    Diagonal(matrix, n - 1, 5);
    Diagonal(matrix, n - m, 6);
    Diagonal(matrix, n - m - 1, 7);
    Diagonal(matrix, n - m - 2, 8);

    x0.resize(n);
    xk.resize(n);
    Ax.resize(n);
    f_Ax.resize(n);
    for (int i = 0; i < n; i++)
    {
        x_0 >> x0[i];
        xk[i] = x0[i];
    }

    f.resize(n);
    for (int i = 0; i < n; i++)
        vector >> f[i];
}

void Iterations(vector<type>& xk, vector<type>& x0)
{
    for (int i = 0; i < n; i++)
    {
        type sum = 0, a = 0;
        for (int j = 0; j < n; j++)
        {
            if (i == j) a = di[4][i]; // главная диагональ
            else
            {
                if (j - i == -m - 2) a = di[0][j]; // нижние диагонали
                else if (j - i == -m - 1) a = di[1][j];
                else if (j - i == -m) a = di[2][j];
                else if (j - i == -1) a = di[3][j];
                else if (j - i == 1) a = di[5][i]; // верхние диагонали
                else if (j - i == m) a = di[6][i];
                else if (j - i == m + 1) a = di[7][i];
                else if (j - i == m + 2) a = di[8][i];
                else  a = 0;
            }
            sum += a * x0[j];
        }
        Ax[i] = sum;
        xk[i] = x0[i] + w / di[4][i] * (f[i] - sum);
    }
}

type Norm(vector<type>& y)
{
    type norma = 0;
    for (int i = 0; i < n; i++)
        norma += y[i] * y[i];
    return sqrt(norma);
}

void Jacobi()
{
    if (w <= 0 || w > 1) return;
    for (int i = 1; i < maxiter; i++)
    {
        Iterations(xk, x0);
        for (int i = 0; i < n; i++)
        {
            f_Ax[i] = f[i] - Ax[i];
            x0[i] = xk[i];
        }

        type nevyazka;
        nevyazka = Norm(f_Ax) / Norm(f);
        if (nevyazka < e)
            return;
    }

}

void Gauss()
{
    if (w <= 0 || w >= 2) return;
    for (int i = 1; i < maxiter; i++)
    {
        Iterations(x0, x0);
        for (int i = 0; i < n; i++)
            f_Ax[i] = f[i] - Ax[i];

        type nevyazka;
        nevyazka = Norm(f_Ax) / Norm(f);
        cout << nevyazka << " " << i << endl;
        if (nevyazka < e)
            return;
    }
}

void Output()
{
    FILE* out;
    out = fopen("out.txt", "w");
    for (int i = 0; i < n; i++)
        fprintf(out, STROUT, x0[i]);
}

int main()
{
    ifstream matrix("matrix.txt");
    ifstream vector("vector.txt");
    ifstream x_0("x_0.txt");

    Input(matrix, vector, x_0);
    //Jacobi();
    clock_t start = clock();
    Gauss();
    clock_t end = clock();
    double cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("\n");
    printf("for loop took %f seconds to execute \n", cpu_time_used);
    Output();

    matrix.close();
    vector.close();
    x_0.close();
}

