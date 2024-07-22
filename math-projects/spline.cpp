#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

int n = 0;
double ht = 0;
double ki = 0;
double x = 0;
double r = 0;

vector<double> xi;
vector<double> fi;
vector<double> hi;
vector<double> fd;

void Input(ifstream& in)
{
    //Количество узлов
    in >> n;
    fi.resize(n);
    xi.resize(n);
    hi.resize(n - 1);
    fd.resize(n);
    //Значения xi
    for (int i = 0; i < n; i++)
        in >> xi[i];
    //Значения f(xi)
    for (int i = 0; i < n; i++)
        in >> fi[i];
}

//Шаги
void add_hi()
{
    for (int i = 0; i < n - 1; i++)
        hi[i] = xi[i + 1] - xi[i];
}

double Psi1(double x, int k)
{
    ki = (x - xi[k]) / hi[k];
    return 1 - 3 * pow(ki, 2) + 2 * pow(ki, 3);
}

double Psi2(double x, int k)
{
    ki = (x - xi[k]) / hi[k];
    return hi[k] * (ki - 2 * pow(ki, 2) + pow(ki, 3));
}

double  Psi3(double x, int k)
{
    ki = (x - xi[k]) / hi[k];
    return 3 * pow(ki, 2) - 2 * pow(ki, 3);
}

double Psi4(double x, int k)
{
    ki = (x - xi[k]) / hi[k];
    return hi[k] * (-pow(ki, 2) + pow(ki, 3));
}

void find_fd()
{
    for (int i = 1; i < n - 1; i++)
    {
        fd[i] = -fi[i - 1] * hi[i] / (hi[i - 1] * (hi[i - 1] + hi[i])) + (hi[i] - hi[i - 1])\
            * fi[i] / (hi[i - 1] * hi[i]) + hi[i - 1] * fi[i + 1] / (hi[i] * (hi[i - 1] + hi[i]));
    }

    fd[0] = -fi[0] * (2 * hi[0] + hi[1]) / (hi[0] * (hi[0] + hi[1]))\
        + fi[1] * (hi[0] + hi[1]) / (hi[0] * hi[1]) - fi[2] * hi[0] / (hi[1] * (hi[0] + hi[1]));

    fd[n - 1] = hi[n - 2] * fi[n - 3] / ((hi[n - 3] + hi[n - 2]) * hi[n - 3])\
        - (hi[n - 3] + hi[n - 2]) * fi[n - 2] / (hi[n - 2] * hi[n - 3])\
        + (2 * hi[n - 2] + hi[n - 3]) * fi[n - 1] / (hi[n - 2] * (hi[n - 3] + hi[n - 2]));
}

void Func()
{
    add_hi();
    find_fd();
    for (int i = 0; i < n - 1; i++)
    {
        ht = hi[i] / 100;
        for (int j = 0; j < 101; j++)
        {
            x = xi[i] + j * ht;
            r = fi[i] * Psi1(x, i) + fi[i + 1] * Psi3(x, i)\
                + fd[i] * Psi2(x, i) + fd[i + 1] * Psi4(x, i);
            cout << setprecision(15) << x << " " << setprecision(15) << r << endl;
        }
    }
}

int main()
{
    ifstream in("input.txt");
    Input(in);
    Func();
}
