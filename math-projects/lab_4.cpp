#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

double h = 1e-8;
int n = 0; // количество переменных
int m = 0; // количество уравнений
int maxiter; //	максимальное количество итераций
double e1 = 0; // точность
double e2 = 0; // точность

vector<double> x0; // начальное приближение
vector<double> xk; // решение на k итерации
vector<double> delta_xk; // решение на k итерации
vector<double> F; // значение вектор-функции
vector<vector<double>> A; // матрица
vector<double> xkv; // для подсчета нормы

double Function(int i, double x, double y) // задаются сами уравнения системы
{
    switch (i)
    {
        case 0: return (x - 5) * (x - 5) + (y - 3) * (y - 3) - 9;
        case 1: return (x - 11) * (x - 11) + (y - 3) * (y - 3) - 9;
        /*case 0: return y+x-5;
        case 1: return y-4;
        case 2: return y-x+5;
        case 0: return y - sin(x);
        case 1: return y - x + 5;*/
    }
    return 0;
}

void Input(ifstream& kuslau, ifstream& x_0)
{
    kuslau >> n >> m >> maxiter >> e1 >> e2;
    x0.resize(n);
    xk.resize(n);
    xkv.resize(n);
    delta_xk.resize(n);
    A.resize(m);
    for (int i = 0; i < m; i++)
        A[i].resize(n);
    F.resize(m);
    kuslau >> x0[0] >> x0[1];
    for (int i = 0; i < n; i++)
    {
        x_0 >> x0[i];
        xk[i] = x0[i];
    }
}

double norma_f(vector<double> xkk)
{
    double norm = 0;
    for (int i = 0; i < m; i++)
        norm += Function(i, xkk[0], xkk[1]) * Function(i, xkk[0], xkk[1]);
    return sqrt(norm);
}

double Diff(int i, int j, double x, double y) // производные для аналитической матрицы Якоби
{
    switch (j) // по какой переменной производная
    {
    case 0: // производная по х
        return (Function(i, x + h, y) - Function(i, x - h, y)) / (2 * h);
    case 1: // производная по y
        return (Function(i, x, y + h) - Function(i, x, y - h)) / (2 * h);
    }
    return 0;
}

void Jacobi_analytical() // аналитический подсчет матрицы Якоби
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
        {
            A[i][j] = Diff(i, j, xk[0], xk[1]);
        }

}

double Diff_function(int i, int j, double x, double y) // производные для численной матрицы Якоби
{
    switch (j) // по какой переменной производная
    {
    case 0:// производная по х
        switch (i) // какая функция
        {
            case 0: return 2 * (x - 5);
            case 1: return 2 * (x - 11);
            /*case 0: return 1;
            case 1: return 0;
            case 0: return -cos(x);
            case 1: return -1;*/

        }
    case 1:// производная по y
        switch (i) // какая функция
        {
            case 0: return 2 * (y - 3);
            case 1: return 2 * (y - 3);
            /*case 0: return 1;
            case 1: return 1;
            case 0: return 1;
            case 1: return 1;*/
        }
    }
    return 0;
}


void Jacobi_numerical() // численный подсчет матрицы Якоби
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
        {
            A[i][j] = Diff_function(i, j, xk[0], xk[1]);

        }

}

void Gauss() // с выбором ведущего элемента
{
    for (int k = 0; k < n - 1; k++)
    {
        double max = abs(A[k][k]);
        int ind = k;
        for (int j = k + 1; j < n; j++) // поиск наибольшего по модулю в столбце
        {
            if (abs(A[j][k]) > max)
            {
                max = abs(A[j][k]);
                ind = j;
            }
        }
        if (max > 0)
        {
            if (ind != k) // замена строк, если ведущий элемент не на диагонали
            {
                swap(A[k], A[ind]);
                swap(F[k], F[ind]);
            }
        }
        if (A[k][k] != 0) // проверка что на диагонали не 0
        {
            for (int j = k + 1; j < n; j++)
            {
                double b = A[j][k] / A[k][k];
                A[j][k] = 0;
                F[j] -= b * F[k];
                for (int i = k + 1; i < n; i++)
                    A[j][i] -= b * A[k][i];
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) // поиск вектора delta_xk
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++)
            sum += A[i][j] * delta_xk[j];
        delta_xk[i] = (F[i] - sum) / A[i][i];
    }
}

void Sort_matrix()
{
    for (int i = 1; i < m; i++)
    {
        int j = i - 1;
        while (j >= 0 && abs(F[j]) < abs(F[j + 1]))
        {
            swap(F[j], F[j + 1]);
            swap(A[j], A[j + 1]);
            j--;
        }
    }
}

void Excluded_expressions() //суммирование исключаемых выражений в квадрате
{
    int i = n - 1;
    for (int k = i; k < m; k++)
    {
        for (int j = 0; j < n; j++)
        {
            A[k][j] *= 2 * Function(k, xk[0], xk[1]);
            F[k] = -pow(F[k], 2);
        }
    }

    for (int k = i + 1; k < m; k++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] += A[k][j];
        }
        F[i] -= F[k];
    }
}

void SNU_v7()
{
    int exit = 1;
    for (int k_iter = 1; k_iter <= maxiter && exit; k_iter++)
    {

        for (int i = 0; i < n; i++)
            F[i] = -Function(i, xk[0], xk[1]);

        Jacobi_numerical();
        Sort_matrix();
        Excluded_expressions();
        Gauss();

        double beta = 1;

        for (int i = 0; i < n; i++)
            xkv[i] = xk[i] + beta * delta_xk[i];

        while (norma_f(xkv) > norma_f(xk) && exit)
        {
            beta /= 2;

            for (int i = 0; i < n; i++)
                xkv[i] = xk[i] + beta * delta_xk[i];

            if (beta < e1) exit = 0;
        }

        if (beta < e1 || norma_f(xk) / norma_f(x0) < e2) exit = 0;

        for (int i = 0; i < n; i++)
            xk[i] = xkv[i];

        cout << k_iter << " " << beta << " " << setprecision(15) << xk[0] << " " << setprecision(15) << xk[1] << " " << norma_f(xk) << endl;
    }
}

void SNU_v3()
{
    int exit = 1;
    for (int k_iter = 1; k_iter <= maxiter && exit; k_iter++)
    {

        for (int i = 0; i < n; i++)
            F[i] = -Function(i, xk[0], xk[1]);

        Jacobi_analytical();
        Sort_matrix();
        Excluded_expressions();
        Gauss();

        double beta = 1;

        for (int i = 0; i < n; i++)
            xkv[i] = xk[i] + beta * delta_xk[i];

        while (norma_f(xkv) > norma_f(xk) && exit)
        {
            beta /= 2;

            for (int i = 0; i < n; i++)
                xkv[i] = xk[i] + beta * delta_xk[i];

            if (beta < e1) exit = 0;
        }

        if (beta < e1 || norma_f(xk) / norma_f(x0) < e2) exit = 0;

        for (int i = 0; i < n; i++)
            xk[i] = xkv[i];

        cout << k_iter << " " << beta << " " << setprecision(15) << xk[0] << " " << setprecision(15) << xk[1] << " " << norma_f(xk) << endl;
    }
}

void SNU_v2() // аналитический якоби, исключение уравнений с наименьшим абсолютным значением
{
    int exit = 1;
    for (int k_iter = 1; k_iter <= maxiter && exit; k_iter++)
    {
        for (int i = 0; i < n; i++)
            F[i] = -Function(i, xk[0], xk[1]);

        Jacobi_analytical();
        Sort_matrix();
        Gauss();

        double beta = 1;

        for (int i = 0; i < n; i++)
            xkv[i] = xk[i] + beta * delta_xk[i];

        while (norma_f(xkv) > norma_f(xk))
        {
            beta /= 2;

            for (int i = 0; i < n; i++)
                xkv[i] = xk[i] + beta * delta_xk[i];
        }

        if (beta < e1 || norma_f(xk) / norma_f(x0) < e2) exit = 0;

        for (int i = 0; i < n; i++)
            xk[i] = xkv[i];

        cout << k_iter << " " << beta << " " << setprecision(15) << xk[0] << " " << setprecision(15) << xk[1] << " " << norma_f(xk) << endl;
    }
}


int main()
{
    ifstream kuslau("kuslau.txt");
    ifstream x_0("x0.txt");

    Input(kuslau, x_0);
    SNU_v3();
    //SNU_v2();
    //SNU_v7();
}
