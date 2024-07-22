#include <math.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <chrono>

typedef double type;
using namespace std;
using namespace chrono;

#define STRIN "%lf\t"
#define STROUT "%.15lf\t"
#pragma warning(disable:4996)

int n = 0; // размер матрицы
int maxiter; // максимальное количество итераций
int k_iter; // фактическое количество итераций
type e = 0; // величина требуемой относительной невязки
vector<int> ig; //указатели начала строк
vector<int> jg; //номера столбцов внедиагональных элементов
vector<type> di; //диагональные элементы матрицы
vector<type> gu; //внедиагональные элементы верхнего треугольника матрицы
vector<type> gl; //внедиагональные элементы нижнего треугольника матрицы
vector<type> pr; //вектор правой части
vector<type> x0; // начальное приближение
vector<type> ggl; //нижний треугольник
vector<type> ggu; //верхний треугольник
vector<type> diu; //диагональ u
vector<type> dil; //диагональ l
vector<type> r; // для реализации функций
vector<type> rk; // для реализации функций
vector<type> Ar; // для реализации функций
vector<type> z; // для реализации функций
vector<type> p; // для реализации функций
vector<type> y; // для реализации функций
vector<type> q; // для реализации функций

void Input(ifstream& kuslau, ifstream& iginput, ifstream& jginput, ifstream& diinput, ifstream& glinput, ifstream& guinput, ifstream& prinput, ifstream& x0input)
{
    kuslau >> n >> maxiter >> e;
    ig.resize(n + 1);
    for (int i = 0; i < n + 1; i++)
        iginput >> ig[i];

    int k = ig[n] - ig[0]; // количество внедиагональных
    jg.resize(k);
    for (int i = 0; i < k; i++)
        jginput >> jg[i];

    di.resize(n);
    for (int i = 0; i < n; i++)
        diinput >> di[i];

    gl.resize(k);
    for (int i = 0; i < k; i++)
        glinput >> gl[i];

    gu.resize(k);
    for (int i = 0; i < k; i++)
        guinput >> gu[i];

    pr.resize(n);
    for (int i = 0; i < n; i++)
        prinput >> pr[i];

    x0.resize(n);
    for (int i = 0; i < n; i++)
        x0input >> x0[i];

    if (ig[0] == 1)
    {
        for (int i = 0; i < n + 1; i++)
            ig[i] -= 1;

        for (int i = 0; i < k; i++)
            jg[i] -= 1;
        k--; // ig[n] - 1;
    }

    ggl.resize(k);
    ggu.resize(k);
    diu.resize(n);
    dil.resize(n);
    r.resize(n);
    rk.resize(n);
    Ar.resize(n);
    z.resize(n);
    p.resize(n);
    y.resize(n);
    q.resize(n);
}

type Norm(vector<type>& y)
{
    type norma = 0;
    for (int i = 0; i < n; i++)
        norma += y[i] * y[i];
    return sqrt(norma);
}

type vector_multiplication(vector<type>& a, vector<type>& b)
{
    type s = 0;
    for (int i = 0; i < n; i++)
        s += a[i] * b[i];
    return s;
}
// Aa=b
void multiplication_matrix_on_vector(vector<type>& x, vector<type>& b)
{
    for (int i = 0; i < n; i++)
    {
        b[i] = di[i] * x[i];
        int i0 = ig[i];
        int i1 = ig[i + 1];
        for (int k = i0; k < i1; k++)
        {
            int j = jg[k];
            b[i] += gl[k] * x[j]; // нижний треугольник
            b[j] += gu[k] * x[i]; // верхний треугольник
        }
    }
}

// d = a + b
void summ(vector<type>& a, vector<type>& b, vector<type>& d)
{
    for (int i = 0; i < n; i++)
        d[i] = a[i] + b[i];
}

// d = a * b
void multCoeff(vector<type>& a, type b, vector<type>& d)
{
    for (int i = 0; i < n; i++)
        d[i] = a[i] * b;
}

void calcLU()
{
    for (int i = 0; i < n; i++)
    {
        type sumDi = 0;
        int i0 = ig[i];
        int i1 = ig[i + 1];
        for (int k = i0; k < i1; k++)
        {
            double suml = 0, sumu = 0;
            int j = jg[k];
            int j0 = ig[j], j1 = ig[j + 1];
            for (int ik = i0, kj = j0; ik < i1 && kj < j1; )
            {
                if (jg[ik] > jg[kj]) kj++;
                else if (jg[ik] < jg[kj]) ik++;
                else
                {
                    suml += ggl[ik] * ggu[kj];
                    sumu += ggl[kj] * ggu[ik];
                    ik++;
                    kj++;
                }
            }
            ggl[k] = (gl[k] - suml) / diu[j];
            ggu[k] = (gu[k] - sumu) / diu[j];
            sumDi += ggl[k] * ggu[k];
        }
        diu[i] = sqrt(di[i] - sumDi);
        dil[i] = diu[i];
    }
}

void CalcY(vector<double> &b, vector<double>& y)
{
    for (int i = 0; i < n; i++)
    {
        type sum = 0;
        int i0 = ig[i];
        int i1 = ig[i + 1];
        for (int k = i0; k < i1; k++)
        {
            int j = jg[k];
            sum += ggl[k] * y[j];
        }
        y[i] = (b[i] - sum) / dil[i];
    }
}

void CalcX(vector<double>& y, vector<double>& x)
{
    vector <double> v = y;

    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = v[i] / diu[i];
        for (int k = ig[i]; k < ig[i + 1]; k++)
        {
            int j = jg[k];
            v[j] -= x[i] * ggu[k];
        }
    }
}

double calcDiscrepancy(vector <double>& x, vector <double>& b, double& normb)
{
    for (int i = 0; i < n; i++)
    {
        normb += b[i] * b[i];
        r[i] = b[i] - di[i] * x[i];
        for (int k = ig[i]; k < ig[i + 1]; k++)
        {
            int j = jg[k];
            r[i] -= gl[k] * x[j];
            r[j] -= gu[k] * x[i];
        }
    }
    return sqrt(vector_multiplication(r, r) / normb);
}


void LOS_LU()
{
    double normb = 0;
    calcDiscrepancy(x0, pr, normb);
    CalcY(r, r);
    CalcX(r, z);
    multiplication_matrix_on_vector(z, p);
    CalcY(p, p);
    double scalarr = vector_multiplication(r, r);
    double discrepancy = sqrt(scalarr / normb);
    for (int k = 1; k < maxiter && discrepancy > e; k++)
    {
        double scalarp = vector_multiplication(p, p);
        double alpha = vector_multiplication(p, r) / scalarp;
        multCoeff(z, alpha, y);
        summ(x0, y, x0);
        multCoeff(p, -alpha, y);
        summ(r, y, r);

        CalcX(r, rk);
        multiplication_matrix_on_vector(rk, Ar);
        CalcY(Ar, q);
        double betta = -vector_multiplication(p, q) / scalarp;
        multCoeff(z, betta, y);
        summ(rk, y, z);
        multCoeff(p, betta, y);
        summ(q, y, p);
        discrepancy = sqrt(vector_multiplication(r, r) / scalarr);
        cout << k << " " << discrepancy << endl;
        k_iter = k;
    }
    normb = 0;
    calcDiscrepancy(x0, pr, normb);
    discrepancy = sqrt(vector_multiplication(r, r) / normb);
    cout << "Final discrepancy: " << discrepancy << endl;
}

//Гильберт
void MakeGilbert(int k, int iter, type nev)
{
    //Построение матрицы Гильберта
    int k1 = (int)((k * k - k) / 2);
    vector <double> di(k), ggl(k1), ggu(k1), x(k), b(k);
    vector <int> ig(k + 1), jg(k1);
    ig[0] = 0;
    for (int i = 0; i < k; i++)
    {
        di[i] = 1.0 / (2 * i + 1);
        x[i] = i + 1;
        ig[i + 1] = ig[i] + i;
        for (int k = ig[i], j = 0; k < ig[i + 1]; k++, j++)
        {
            jg[k] = j;
            double elem = 1.0 / (i + j + 1);
            ggl[k] = elem;
            ggu[k] = elem;
        }
    }

    for (int i = 0; i < k; i++)
    {
        b[i] = di[i] * x[i];
        for (int k = ig[i]; k < ig[i + 1]; k++)
        {
            int j = jg[k];
            b[i] += ggl[k] * x[j];
            b[j] += ggu[k] * x[i];
        }
    }

    FILE* kuslau = fopen("kuslau.txt", "w");
    fprintf(kuslau, "%d %d %.15lf", k, iter, nev);
    fclose(kuslau);

    FILE* dia = fopen("di.txt", "w");
    for (int i = 0; i < k; i++)
        fprintf(dia, "%.15lf\n", di[i]);
    fclose(dia);

    FILE* gl = fopen("gl.txt", "w");
    for (int i = 0; i < k1; i++)
        fprintf_s(gl, "%.15lf\n", ggl[i]);
    fclose(gl);

    FILE* gu = fopen("gu.txt", "w");
    for (int i = 0; i < k1; i++)
        fprintf_s(gl, "%.15lf\n", ggu[i]);
    fclose(gu);

    FILE* igf = fopen("ig.txt", "w");
    for (int i = 0; i < k + 1; i++)
        fprintf_s(igf, "%d\n", ig[i]);
    fclose(igf);

    FILE* jgf = fopen("jg.txt", "w");
    for (int i = 0; i < k1; i++)
        fprintf_s(jgf, "%d\n", jg[i]);
    fclose(jgf);

    FILE* pr = fopen("pr.txt", "w");
    for (int i = 0; i < k; i++)
        fprintf_s(pr, "%.15lf\n", b[i]);
    fclose(pr);

    FILE* x0 = fopen("x0.txt", "w");
    for (int i = 0; i < k; i++)
        fprintf_s(x0, "%.15lf\n", 0.0);
    fclose(x0);
}

void Output()
{
    FILE* out;
    out = fopen("out.txt", "w");
    for (int i = 0; i < n; i++)
        fprintf(out, "%.15lf\n", x0[i]);
    fclose(out);
}

int main()
{
    ifstream kuslau("kuslau.txt");
    ifstream iginput("ig.txt");
    ifstream jginput("jg.txt");
    ifstream diinput("di.txt");
    ifstream glinput("gl.txt");
    ifstream guinput("gu.txt");
    ifstream prinput("pr.txt");
    ifstream x0input("x0.txt");

    //MakeGilbert(10, 10000, 1e-15);
    Input(kuslau, iginput, jginput, diinput, glinput, guinput, prinput, x0input);
    calcLU();
    auto start = steady_clock::now();
    LOS_LU();
    auto end = steady_clock::now();
    double time_used = duration<double>(end - start).count();
    cout << endl;
    cout << "Time: " << time_used << endl;
    Output();
}
