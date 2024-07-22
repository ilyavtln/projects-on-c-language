#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;

// Структура для хранения КЭ
struct rectangle
{
    //Номер элемента
    int number;
    //Координаты для х
    int x_left, x_right;
    //Координаты для y
    int y_down, y_up;
};

// Структура для хранения шагов
struct steps
{
    double deltaT = 0;
    double deltaT0 = 0, deltaT1 = 0, deltaT2 = 0;
};

// Структура для Гаусса
struct QuadratureNode
{
    double x, y;
    double weight;
};

vector<double> X_w; //Массив X без разбиений
vector<double> Y_w; //Массив Y без разбиений
vector<double> T_w; //Массив T без разбиений

int Nx_w, Ny_w, Nt_w; //Длины массивов Xw, Yw и Tw

int L; //Количество подобластей
vector<rectangle> rectangles; //Массив подобластей

vector<double> X; //Массив X c разбиениями
vector<double> Y; //Массив Y с разбиениями
vector<double> T; //Массив T с разбиениями

int Nx, Ny, Nt; //Размеры массивов X, Y и T

vector<int> IX; //Массив индекса элементов в Х, которые есть в X_w
vector<int> IY; //Массив индекса элементов в Y, которые есть в Y_w

int countKraev; //Количество краевых условий
vector<vector<int>> MatrixKraev; //Матрица краевых условий

//Глобальная матрица
vector<double> di;
vector<int> ig;
vector<int> jg;
vector<double> ggl;
vector<double> ggu;
int countEl; //Количество элементов в ggl или ggu
int globalN; //Размерность глобальной матрицы (Nx * Ny)

//Глобальный вектор
vector<double> pr;

//ЛОС 
int maxiter = 10000; // максимальное количество итераций
int k_iter = 0; // фактическое количество итераций
double e = 1e-32; // величина требуемой относительной невязки

vector<double> z;
vector<double> r;
vector<double> p;
vector<double> y; // для реализации функций
vector<double> t; // для реализации функций

//Значения весов на временных слоях от j-3 до j-го
vector<double> q_j3; //Вектор - ответ
vector<double> q_j2; //Вектор - ответ
vector<double> q_j1; //Вектор - ответ
vector<double> q_j; //Вектор - ответ

// Искомая функция
double functionU(double x, double y, double t, int i)
{
    switch (i)
    {
    case 1: return 5 * x * t + 2 * t;
    default: return 0;
    }
}

// Функция правой части
double func(double x, double y, double t, int i)
{
    switch (i)
    {
    case 1: return 10 * x + 4;
    default: cout << "Такой функции не существет" << endl;
    }
    return -1;
}

double lambda(double x, double y, double t, int i)
{
    switch (i)
    {
    case 1: return 1;
    default: cout << "Такой лямбды не существет" << endl;
    }
    return -1;
}

double sigma(double x, double y, double t, int i)
{
    switch (i)
    {
    case 1: return 2;
    default: cout << "Такой сигмы не существет" << endl;
    }
    return -1;
}

//1 краевое
double UgKr(int i, double x, double y, double t)
{
    switch (i)
    {
    case 1: return 5 * x * t + 2 * t;
    case 2: return 50 * t + 2 * t;
    default: cout << "Такого 1го краевого не существет" << endl;
    }
    return -1;
}

//2 краевое
double TettaKr(int i, double x, double y, double t)
{
    switch (i)
    {
    case 1: return 0;
    default: cout << "Такой тетты не существет" << endl;
    }
    return -1;
}

//3 краевое
double bettaKr(int i)
{
    switch (i)
    {
    case 1: return 1;
    default: cout << "Такой бетты не существет" << endl; ;
    }
    return -1;
}

//3 краевое
double UbettaKr(int i, double x, double y, double t)
{
    switch (i)
    {
    case 1: return -3 * t;
    default: cout << "Такой Убетты не существет" << endl;
    }
    return -1;
}

// Функция u при t в начальный момент времени
double u0(double x, double y, int i)
{
    return functionU(x, y, T[0], i);
}

// Функция u при t в момент времени tvalue
double uAtTime(double x, double y, double tValue, int i)
{
    return functionU(x, y, tValue, i);
}

double lambdaMiddle(vector<vector<double>> nodes, double curT, int number)
{
    double valueSum = 0;
    for (int i = 0; i < 4; i++)
    {
        valueSum += lambda(nodes[i][0], nodes[i][1], curT, number);
    }
    return valueSum / 4.0;
}

double sigmaMiddle(vector<vector<double>> nodes, double curT, int number)
{
    double valueSum = 0;
    for (int i = 0; i < 4; i++)
    {
        valueSum += sigma(nodes[i][0], nodes[i][1], curT, number);
    }
    return valueSum / 4.0;
}

vector<vector<double>> locaMatrixM(double gammaAvr, double hx, double hy)
{
    vector<vector<double>> M(4, vector<double>(4));
    double generalCoeff = (gammaAvr * hx * hy) / 36;
    double valueDi = 4 * generalCoeff;
    double value_0_1 = 2 * generalCoeff;
    double value_0_2 = 2 * generalCoeff;
    double value_0_3 = generalCoeff;
    double value_1_2 = 1 * generalCoeff;
    double value_1_3 = 2 * generalCoeff;
    double value_2_3 = 2 * generalCoeff;
    //Диагональ
    M[0][0] = valueDi;
    M[1][1] = valueDi;
    M[2][2] = valueDi;
    M[3][3] = valueDi;
    //1 строка and 1 столбец 
    M[0][1] = value_0_1;
    M[1][0] = value_0_1;
    M[0][2] = value_0_2;
    M[2][0] = value_0_2;
    M[0][3] = value_0_3;
    M[3][0] = value_0_3;
    //2 строка and 2 столбец 
    M[1][2] = value_1_2;
    M[2][1] = value_1_2;
    M[1][3] = value_1_3;
    M[3][1] = value_1_3;
    //3 строка and 3 столбец
    M[2][3] = value_2_3;
    M[3][2] = value_2_3;

    return M;
}

vector<vector<double>> locaMatrixG(double lambdaAvr, double hx, double hy)
{
    double lambdaDiv6 = lambdaAvr / 6;
    double coeff1 = lambdaDiv6 * (hy / hx);
    double coeff2 = lambdaDiv6 * (hx / hy);
    vector<vector<double>> G(4, vector<double>(4));
    double valueDi = 2 * (coeff1 + coeff2);

    //Диагональ
    G[0][0] = valueDi;
    G[1][1] = valueDi;
    G[2][2] = valueDi;
    G[3][3] = valueDi;
    //1 строка and 1 столбец 
    G[0][1] = -2 * coeff1 + coeff2;
    G[1][0] = -2 * coeff1 + coeff2;
    G[0][2] = coeff1 + -2 * coeff2;
    G[2][0] = coeff1 + -2 * coeff2;
    G[0][3] = -coeff1 - coeff2;
    G[3][0] = -coeff1 - coeff2;
    //2 строка and 2 столбец 
    G[1][2] = -coeff1 - coeff2;
    G[2][1] = -coeff1 - coeff2;
    G[1][3] = coeff1 - 2 * coeff2;
    G[3][1] = coeff1 - 2 * coeff2;
    //3 строка and 3 столбец
    G[2][3] = -2 * coeff1 + coeff2;
    G[3][2] = -2 * coeff1 + coeff2;

    return G;
}

vector<double> locaVectorB(double hx, double hy, double f1, double f2, double f3, double f4)
{
    double generalCoeff = (hx * hy) / 36;
    vector<double> b(4);
    b[0] = generalCoeff * (4 * f1 + 2 * f2 + 2 * f3 + f4);
    b[1] = generalCoeff * (2 * f1 + 4 * f2 + f3 + 2 * f4);
    b[2] = generalCoeff * (2 * f1 + f2 + 4 * f3 + 2 * f4);
    b[3] = generalCoeff * (f1 + 2 * f2 + 2 * f3 + 4 * f4);

    return b;
}

void inputGrid()
{
    ifstream in("grid.txt");
    if (!in) exit(-1);
    //Чтение данных об области по x
    in >> Nx_w;
    X_w.resize(Nx_w);
    for (int i = 0; i < Nx_w; i++)
    {
        in >> X_w[i];
    }
    //Чтение данных об области по y
    in >> Ny_w;
    Y_w.resize(Ny_w);
    for (int i = 0; i < Ny_w; i++)
    {
        in >> Y_w[i];
    }
    //Чтение информации о подобластях
    in >> L;
    rectangles.resize(L);
    for (int i = 0; i < L; i++)
    {
        in >> rectangles[i].number 
           >> rectangles[i].x_left
           >> rectangles[i].x_right
           >> rectangles[i].y_down
           >> rectangles[i].y_up;
    }
    in.close();
}

void inputPartitions()
{
    ifstream in("partitions.txt");
    if (!in) exit(-1);
    IX.resize(Nx_w);
    IY.resize(Ny_w);

    //Построение разбиения по х
    int nXk = 0;
    X.resize(1, X_w[0]);
    for (int i = 0, j = 1; i < Nx_w - 1; i++, j++)
    {
        int countInterval = 0;
        double coeff = 0, step = 0;
        in >> countInterval >> coeff;
        nXk += countInterval;
        X.resize(nXk + 1);
        if (coeff != 1)
        {
            double sumProgression = (pow(coeff, countInterval) - 1.0) / (coeff - 1.0);
            step = (X_w[i + 1] - X_w[i]) / sumProgression;
            int jk = 1;
            for (j; j < nXk; j++, jk++)
            {
                X[j] = X_w[i] + step * (pow(coeff, jk) - 1.0) / (coeff - 1.0);
            }
        }
        else
        {
            step = (X_w[i + 1] - X_w[i]) / countInterval;
            int jk = 1;
            for (j; j < nXk; j++, jk++)
            {
                X[j] = X_w[i] + step * jk;
            }
        }
        X[j] = X_w[i + 1];
        IX[i + 1] = j;
        Nx = nXk + 1;
    }

    //Построение разбиения по y
    int nYk = 0;
    Y.resize(1, Y_w[0]);
    for (int i = 0, j = 1; i < Ny_w - 1; i++, j++)
    {
        int countInterval = 0;
        double coeff = 0, step = 0;
        in >> countInterval >> coeff;
        nYk += countInterval;
        Y.resize(nYk + 1);
        if (coeff != 1)
        {
            double sumProgression = (pow(coeff, countInterval) - 1.0) / (coeff - 1.0);
            step = (Y_w[i + 1] - Y_w[i]) / sumProgression;
            int jk = 1;
            for (j; j < nYk; j++, jk++)
            {
                Y[j] = Y_w[i] + step * (pow(coeff, jk) - 1.0) / (coeff - 1.0);
            }
        }
        else
        {
            step = (Y_w[i + 1] - Y_w[i]) / countInterval;
            int jk = 1;
            for (j; j < nYk; j++, jk++)
            {
                Y[j] = Y_w[i] + step * jk;
            }
        }
        Y[j] = Y_w[i + 1];
        IY[i + 1] = j;
        Ny = nYk + 1;
    }
}

void inputTimeLayers()
{
    ifstream in("time_grid.txt");
    if (!in) exit(-1);

    //Чтение данных о временных слоях
    in >> Nt_w;
    T_w.resize(Nt_w);
    for (int i = 0; i < Nt_w; i++)
    {
        in >> T_w[i];
    }

    //Построение разбиения по t
    int nTk = 0;
    T.resize(1, T_w[0]);
    for (int i = 0, j = 1; i < Nt_w - 1; i++, j++)
    {
        int countInterval = 0;
        double coeff = 0, step = 0;
        in >> countInterval >> coeff;
        nTk += countInterval;
        T.resize(nTk + 1);
        if (coeff != 1)
        {
            double sumProgression = (pow(coeff, countInterval) - 1.0) / (coeff - 1.0);
            step = (T_w[i + 1] - T_w[i]) / sumProgression;
            int jk = 1;
            for (j; j < nTk; j++, jk++)
            {
                T[j] = T_w[i] + step * (pow(coeff, jk) - 1.0) / (coeff - 1.0);
            }
        }
        else
        {
            step = (T_w[i + 1] - T_w[i]) / countInterval;
            int jk = 1;
            for (j; j < nTk; j++, jk++)
            {
                T[j] = T_w[i] + step * jk;
            }
        }
        T[j] = T_w[i + 1];
        Nt = nTk + 1;
    }
}

void inputKraev()
{
    ifstream in("kraevie.txt");
    in >> countKraev;
    MatrixKraev.resize(countKraev);
    for (int i = 0; i < countKraev; i++)
    {
        MatrixKraev[i].resize(6);
        for (int j = 0; j < 6; j++)
        {
            in >> MatrixKraev[i][j];
        }
    }
}

int numberRectangle(int s, int p)
{
    for (int i = 0; i < L; i++)
    {
        int x0 = IX[rectangles[i].x_left - 1], x1 = IX[rectangles[i].x_right - 1];
        int y0 = IY[rectangles[i].y_down - 1], y1 = IY[rectangles[i].y_up - 1];
        if ((p >= x0 && p <= x1) && (s >= y0 && s <= y1))
            if ((p + 1 >= x0 && p + 1 <= x1) && (s >= y0 && s <= y1))
                if ((p >= x0 && p <= x1) && (s + 1 >= y0 && s + 1 <= y1))
                    if ((p + 1 >= x0 && p + 1 <= x1) && (s + 1 >= y0 && s + 1 <= y1))
                        return rectangles[i].number;
    }
    return -1;
}

int numberAreaForPoint(double x, double y)
{
    for (int i = 0; i < L; i++)
    {
        double x0 = X[IX[rectangles[i].x_left - 1]], x1 = X[IX[rectangles[i].x_right - 1]];
        double y0 = Y[IY[rectangles[i].y_down - 1]], y1 = Y[IY[rectangles[i].y_up - 1]];
        if ((x >= x0 && x <= x1) && (y >= y0 && y <= y1))
            return rectangles[i].number;
    }
    return -1;
}

bool IsFictionNode(int s, int p)
{
    for (int i = 0; i < L; i++)
    {
        int x0 = IX[rectangles[i].x_left - 1], x1 = IX[rectangles[i].x_right - 1];
        int y0 = IY[rectangles[i].y_down - 1], y1 = IY[rectangles[i].y_up - 1];
        if (((p >= x0) && (p <= x1) && (s >= y0) && (s <= y1)))
            return false;
    }
    return true;
}

void generatePortrait()
{
    globalN = Nx * Ny;
    vector<set<int>> list;
    list.resize(globalN);
    di.resize(globalN); 

    for (int s = 0; s < Ny - 1; s++)
    {
        for (int p = 0; p < Nx - 1; p++)
        {
            vector<int> globalNumber{ (s + 1) * Nx + p + 1, (s + 1) * Nx + p,
                                       s * Nx + p + 1, s * Nx + p };

            for (int i = 0; i < 4; i++)
            {
                int ind1 = globalNumber[i];
                for (int j = i + 1; j < 4; j++)
                {
                    int ind2 = globalNumber[j];
                    list[ind1].insert(ind2);
                }
            }
        }
    }

    ig.resize(globalN + 1);
    ig[0] = 0;
    ig[1] = 0;
    for (int i = 0; i < globalN; i++)
    {
        ig[i + 1] = ig[i] + list[i].size();
    }
    countEl = ig[globalN];
    jg.resize(countEl);
    ggl.resize(countEl);
    ggu.resize(countEl);
    for (int i = 0, k = 0; i < globalN; i++)
    {
        for (int j : list[i])
        {
            jg[k] = j;
            k++;
        }
    }

    pr.resize(globalN);

    q_j.resize(globalN);
    q_j1.resize(globalN);
    q_j2.resize(globalN);
    q_j3.resize(globalN);
}

void addLocalElement(double elem, int i, int j)
{
    if (i == j)
    {
        di[i] += elem;
    }
    else
    {
        if (i > j)
        {
            for (int ind = ig[i]; ind < ig[i + 1]; ind++)
            {
                if (jg[ind] == j)
                {
                    ggl[ind] += elem;
                }
            }
        }
        else
        {
            for (int ind = ig[j]; ind < ig[j + 1]; ind++)
            {
                if (jg[ind] == i)
                {
                    ggu[ind] += elem;
                }
            }
        }
    }
}

double valueFuncAtPoint(double x, double y, vector<double>& q)
{
    if (x < X[0] || x > X[Nx - 1] || y < Y[0] || y > Y[Ny - 1])
    {
        cout << "Точка не принадлежит области" << endl;
        return -1;
    }

    int begX = 0, endX = Nx - 1, begY = 0, endY = Ny - 1;
    while (!(X[begX] <= x && x <= X[begX + 1]))
    {
        int indX = (begX + endX) / 2.0;
        if (X[indX] < x)
            begX = indX;
        else
            endX = indX;
    }
    while (!(Y[begY] <= y && y <= Y[begY + 1]))
    {
        int indY = (begY + endY) / 2.0;
        if (Y[indY] < y)
            begY = indY;
        else
            endY = indY;
    }
    vector <int> globalNumbers = { begY * Nx + begX, begY * Nx + begX + 1, (begY + 1) * Nx + begX, (begY + 1) * Nx + begX + 1 };
    double x0 = X[begX], x1 = X[begX + 1], y0 = Y[begY], y1 = Y[begY + 1];
    double hx = x1 - x0, hy = y1 - y0;
    double valueFunc = q[globalNumbers[0]] * ((x1 - x) / hx) * ((y1 - y) / hy) +
        q[globalNumbers[1]] * ((x - x0) / hx) * ((y1 - y) / hy) +
        q[globalNumbers[2]] * ((x1 - x) / hx) * ((y - y0) / hy) +
        q[globalNumbers[3]] * ((x - x0) / hx) * ((y - y0) / hy);
    return valueFunc;
}

void getQuadratures2D(vector<QuadratureNode>& quadratures)
{
    double p1 = sqrt((3. - 2. * sqrt(1.2)) / 7.);
    double p2 = sqrt((3. + 2. * sqrt(1.2)) / 7.);

    double w1 = (18. + sqrt(30)) / 36.;
    double w2 = (18. - sqrt(30)) / 36.;

    vector<double> quadratures1D = { -p1, -p2, p1, p2 };
    vector<double> weights1D = { w1, w2, w1, w2 };


    QuadratureNode node{ };

    for (int s = 0; s < 4; s++)
    {
        for (int p = 0; p < 4; p++)
        {
            node.x = quadratures1D[p];
            node.y = quadratures1D[s];
            node.weight = weights1D[p] * weights1D[s];

            quadratures.push_back(node);
        }
    }
}

double uNumMinusUReal(double x, double y, double t, int i, vector<double> q)
{
    double u = valueFuncAtPoint(x, y, q);
    double uR = functionU(x, y, t, i);

    return (u - uR) * (u - uR);
}

double normL2(double time, vector<double>& q)
{

    vector<QuadratureNode> nodes{ };
    getQuadratures2D(nodes);

    double normL2 = 0;

    for (int s = 0; s < Ny - 1; s++)
    {
        for (int p = 0; p < Nx - 1; p++)
        {
            int areaNumber = numberRectangle(s, p);

            if (areaNumber != -1) 
            {
                double hx = X[p + 1] - X[p];
                double hy = Y[s + 1] - Y[s];

                double integralValue = 0;

                for (auto& node : nodes)
                {
                    double x = hx / 2. * node.x + (X[p] + X[p + 1]) / 2.;
                    double y = hy / 2. * node.y + (Y[s] + Y[s + 1]) / 2.;

                    integralValue += node.weight * uNumMinusUReal(x, y, time, areaNumber, q);
                }
                normL2 += hx * hy * integralValue / 8.;
            };
        }
    }
    return sqrt(normL2);
}

void clearVector(vector<double> &curQ)
{
    for (int i = 0; i < globalN; i++)
    {
        curQ[i] = 0;
    }
}

void clearMatrix()
{
    for (int i = 0; i < globalN; i++)
    {
        di[i] = 0;
        pr[i] = 0;
    }

    for (int i = 0; i < countEl; i++)
    {
        ggl[i] = 0;
        ggu[i] = 0;
    }
}

vector<double> mult_mv_with_coeff(vector<vector<double>>& a, vector<double>& b, double coeff)
{
    vector<double> answer(4, 0);
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            answer[i] += coeff * a[i][j] * b[j];
        }
    }

    return answer;
}

// a = b
void vector_eq_vector(vector<double>& a, vector<double>& b)
{
    for (int i = 0; i < 4; i++)
    {
        a[i] = b[i];
    }
}

vector<double> summ_two_vectors(vector<double> a, vector<double> b)
{
    vector<double> answer(4, 0);
    for (int i = 0; i < 4; i++)
    {
        answer[i] = a[i] + b[i];
    }

    return answer;
}

vector<double> summ_three_vectors(vector<double> a, vector<double> b, vector<double> c)
{
    vector<double> answer(4, 0);
    for (int i = 0; i < 4; i++)
    {
        answer[i] = a[i] + b[i] + c[i];
    }

    return answer;
}

double setCoeffForMatrix(int typeScheme, steps schemeSteps)
{
    double coeffForMA = 1;
    switch (typeScheme)
    {
    case 2:
    {
        double deltaT = schemeSteps.deltaT;
        coeffForMA = 1.0 / deltaT;
        break;
    }
    case 3:
    {
        double deltaT = schemeSteps.deltaT;
        double deltaT0 = schemeSteps.deltaT0;
        coeffForMA = (deltaT + deltaT0) / (deltaT * deltaT0);
        break;
    }
    case 4:
    {
        double deltaT = schemeSteps.deltaT;
        double deltaT0 = schemeSteps.deltaT0;
        double deltaT1 = schemeSteps.deltaT1;
        double deltaT2 = schemeSteps.deltaT2;
        coeffForMA = (deltaT0 * (deltaT0 + deltaT1) + deltaT0 * deltaT + deltaT * (deltaT0 + deltaT1)) / (deltaT * deltaT0 * (deltaT0 + deltaT1));
        break;
    }
    default:
        cout << "Такой схемы не существует" << endl;
        break;
    }

    return coeffForMA;
}

double findAdditiveB(int typeScheme, steps schemeSteps, double elemM, int index)
{
    double incideB = 0;
    switch (typeScheme)
    {
    case 2:
    {
        double deltaT = schemeSteps.deltaT;

        double n1 = 1.0 / deltaT;

        incideB = elemM * n1 * q_j3[index];
        break;
    }
    case 3:
    {
        double deltaT = schemeSteps.deltaT;
        double deltaT0 = schemeSteps.deltaT0;
        double deltaT1 = schemeSteps.deltaT1;

        double n1 = -deltaT0 / (deltaT * deltaT1);
        double n2 = deltaT / (deltaT1 * deltaT0);

        incideB = elemM * (n1 * q_j3[index] + n2 * q_j2[index]);
        break;
    }
    case 4:
    {
        double deltaT = schemeSteps.deltaT;
        double deltaT0 = schemeSteps.deltaT0;
        double deltaT1 = schemeSteps.deltaT1;
        double deltaT2 = schemeSteps.deltaT2;

        double n1 = -1.0 * (((deltaT0 + deltaT1) * deltaT) / ((deltaT1 + deltaT2) * deltaT1 * deltaT0));
        double n2 = ((deltaT0 * deltaT) / (deltaT2 * deltaT1 * (deltaT0 + deltaT1)));
        double n3 = -1.0 * ((deltaT0 * (deltaT0 + deltaT1)) / (deltaT * deltaT2 * (deltaT1 + deltaT2)));

        incideB = -elemM * (n3 * q_j3[index] + n2 * q_j2[index] + n1 * q_j1[index]);
        break;
    }
    default:
        cout << "Такой схемы не существует" << endl;
        break;
    }

    return incideB;
}

void makeGlobalMatrixAndVector(double curT, int schemeType, steps schemeSteps)
{
    for (int s = 0; s < Ny - 1; s++)
    {
        for (int p = 0; p < Nx - 1; p++)
        {
            int areaNumber = numberRectangle(s, p);

            vector<int> globalNum = { Nx * s + p, Nx * s + p + 1, Nx * (s + 1) + p, Nx * (s + 1) + p + 1 };

            if (areaNumber != -1)
            {
                //Находим значения x и y на КЭ
                double
                    x0 = X[p],
                    x1 = X[p + 1],
                    y0 = Y[s],
                    y1 = Y[s + 1];
                //Вычисляем шаги
                double hx = x1 - x0, hy = y1 - y0;
                //Матрица узлов, описывающая КЭ
                vector<vector<double>> nodes =
                {
                    {x0, y0},
                    {x1, y0},
                    {x0, y1},
                    {x1, y1},
                };
                //Находим значения лямбды и гаммы
                double avrLambda = lambdaMiddle(nodes, curT, areaNumber);
                double avrSigma = sigmaMiddle(nodes, curT, areaNumber);
                //Считаем функцию
                double f1 = func(x0, y0, curT, areaNumber),
                       f2 = func(x1, y0, curT, areaNumber),
                       f3 = func(x0, y1, curT, areaNumber),
                       f4 = func(x1, y1, curT, areaNumber);
                //Находим локальную матрицу жесткости
                vector<vector<double>> localG = locaMatrixG(avrLambda, hx, hy);
                //Находим локальную матрицу массы
                vector<vector<double>> localM = locaMatrixM(avrSigma, hx, hy);
                //Находим локальный вектор правой части
                vector<double> localB = locaVectorB(hx, hy, f1, f2, f3, f4);
                // Находим коэффициент для матрицы жесткости
                double coeffForMA = setCoeffForMatrix(schemeType, schemeSteps);

                //Добавление в глобальную матрицу и вектор
                for (int i = 0; i < 4; i++)
                {
                    double sumBi = 0;
                    for (int j = 0; j < 4; j++)
                    {
                        double elemM = localM[i][j];
                        double elemAij = localG[i][j] + coeffForMA * elemM;
                        addLocalElement(elemAij, globalNum[i], globalNum[j]);
                        sumBi += findAdditiveB(schemeType, schemeSteps, elemM, globalNum[j]);
                    }
                    pr[globalNum[i]] += localB[i] + sumBi;
                }
            }
            else
            {
                vector<vector<int>> localNum = { {p, s}, {p + 1, s}, {p, s + 1}, {p + 1, s + 1} };
                for (int i = 0; i < 4; i++)
                {
                    if(IsFictionNode(localNum[i][1], localNum[i][0]))
                    {
                        di[globalNum[i]] = 1;
                    }
                }
            }
        }
    }
}

bool compareRows(const vector<int>& row1, const vector<int>& row2)
{
    return row1[0] > row2[0];  //Сортируем строки, чтобы краевые 1 го рода оказались в конце
}

void kraevoe_3(vector<int> &kr_3, double curT)
{
    //Вертикаль
    if (kr_3[2] == kr_3[3])
    {
        for (int k = IY[kr_3[4] - 1]; k < IY[kr_3[5] - 1]; k++)
        {
            double hy = Y[k + 1] - Y[k]; //Шагы
            double betta = bettaKr(kr_3[1]); //номер формулы для бетты
            double u_betta_1 = UbettaKr(kr_3[1], X[IX[kr_3[2] - 1]], Y[k], curT); //Ubetta1 
            double u_betta_2 = UbettaKr(kr_3[1], X[IX[kr_3[2] - 1]], Y[k + 1], curT); //Ubetta2
            double coeff = (betta * hy) / 6.0;
            vector<vector<double>> localA =
            {
                {2 * coeff, coeff },
                {coeff, 2 * coeff }
            };

            vector<double> localB =
            {
                {coeff * (2 * u_betta_1 + u_betta_2)},
                {coeff * (u_betta_1 + 2 * u_betta_2)}
            };

            vector<int> globalNum = { Nx * k + IX[kr_3[2] - 1], Nx * (k + 1) + IX[kr_3[2] - 1]};

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    double elemAij = localA[i][j];
                    addLocalElement(elemAij, globalNum[i], globalNum[j]);
                }
                pr[globalNum[i]] += localB[i];
            }
        }
    }
    //Горизонталь
    else
    {
        for (int k = IX[kr_3[2] - 1]; k < IX[kr_3[3] - 1]; k++)
        {
            double hx = X[k + 1] - X[k]; //Шаг
            double betta = bettaKr(kr_3[1]); //номер формулы для бетты
            double u_betta_1 = UbettaKr(kr_3[1], X[k], Y[IY[kr_3[4] - 1]], curT), u_betta_2 = UbettaKr(kr_3[1], X[k + 1], Y[IY[kr_3[4] - 1]], curT);
            double coeff = (betta * hx) / 6.0;
            vector<vector<double>> localA =
            {
                {2 * coeff, coeff },
                {coeff, 2 * coeff }
            };

            vector<double> localB =
            {
                {coeff * (2 * u_betta_1 + u_betta_2)},
                {coeff * (u_betta_1 + 2 * u_betta_2)}
            };

            vector<int> globalNum = { Nx * (IY[kr_3[4] - 1]) + k, Nx * (IY[kr_3[4] - 1]) + k + 1};

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    addLocalElement(localA[i][j], globalNum[i], globalNum[j]);
                }
                pr[globalNum[i]] += localB[i];
            }
        }
    }
}

void kraevoe_2(vector<int>& kr_2, double curT)
{
    if (kr_2[2] == kr_2[3])
    {
        for (int k = IY[kr_2[4] - 1]; k < IY[kr_2[5] - 1]; k++)
        {
            double hy = Y[k + 1] - Y[k]; //Шаг
            double u_tetta_1 = TettaKr(kr_2[1], X[IX[kr_2[2] - 1]], Y[k], curT), u_tetta_2 = TettaKr(kr_2[1], X[IX[kr_2[2] - 1]], Y[k + 1], curT);
            double coeff = hy / 6.0;

            vector<double> localB =
            {
                {coeff * (2 * u_tetta_1 + u_tetta_2)},
                {coeff * (u_tetta_1 + 2 *  u_tetta_2)}
            };

            vector<int> globalNum = { Nx * k + IX[kr_2[2] - 1], Nx * (k + 1) + IX[kr_2[2] - 1] };

            for (int i = 0; i < 2; i++)
            {
                pr[globalNum[i]] += localB[i];
            }
        }
    }
    else
    {
        for (int k = IX[kr_2[2] - 1]; k < IX[kr_2[3] - 1]; k++)
        {
            double hx = X[k + 1] - X[k]; //Шаг
            double u_tetta_1 = TettaKr(kr_2[1], X[k], Y[IY[kr_2[4] - 1]], curT), u_tetta_2 = TettaKr(kr_2[1], X[k + 1], Y[IY[kr_2[4] - 1]], curT);
            double coeff = hx / 6.0;

            vector<double> localB =
            {
                {coeff * (2 * u_tetta_1 + u_tetta_2)},
                {coeff * (u_tetta_1 + 2 * u_tetta_2)}
            };

            vector<int> globalNum = { Nx * (IY[kr_2[4] - 1]) + k, Nx * (IY[kr_2[4] - 1]) + k + 1 };

            for (int i = 0; i < 2; i++)
            {
                pr[globalNum[i]] += localB[i];
            }
        }
    }
}

void kraevoe_1(vector<int>& kr_1, double curT)
{
    //Вертикаль
    if (kr_1[2] == kr_1[3])
    {
        for (int k = IY[kr_1[4] - 1]; k <= IY[kr_1[5] - 1]; k++)
        {
            double u_g = UgKr(kr_1[1], X[IX[kr_1[2] - 1]], Y[k], curT);

            int globalNum = Nx * k + IX[kr_1[2] - 1];
            //Зануляем строку в нижнем треугольнике
            for (int i = ig[globalNum]; i < ig[globalNum + 1]; i++)
            {
                ggl[i] = 0;
            }

            //Зануляем столбец в верхнем треугольнике
            for (int j = ig[globalNum + 1]; j < countEl; j++)
            {
                if (jg[j] == globalNum)
                {
                    ggu[j] = 0;
                }
            }
            pr[globalNum] = u_g;
            di[globalNum] = 1;
        }
    }
    //Горизонталь
    else
    {
        for (int k = IX[kr_1[2] - 1]; k <= IX[kr_1[3] - 1]; k++)
        {
            double u_g = UgKr(kr_1[1], X[k], Y[IY[kr_1[4] - 1]], curT);

            int globalNum = Nx * IY[kr_1[4] - 1] + k;
            //Зануляем строку в нижнем треугольнике
            for (int i = ig[globalNum]; i < ig[globalNum + 1]; i++)
            {
                ggl[i] = 0;
            }

            //Зануляем столбец в верхнем треугольнике
            for (int j = ig[globalNum + 1]; j < countEl; j++)
            {
                if (jg[j] == globalNum)
                {
                    ggu[j] = 0;
                }
            }
            pr[globalNum] = u_g;
            di[globalNum] = 1;
        }
    }
}

void boundaryConditions(double curT)
{
    //Сортируем в порядке убывания
    sort(MatrixKraev.begin(), MatrixKraev.end(), compareRows);
    for (int i = 0; i < countKraev; i++) 
    {
        if (MatrixKraev[i][0] == 3)
        {
            kraevoe_3(MatrixKraev[i], curT);
        }

        else if (MatrixKraev[i][0] == 2)
        {
            kraevoe_2(MatrixKraev[i], curT);
        }

        else
        {
            kraevoe_1(MatrixKraev[i], curT);
        }
    }
}

double Norm(vector<double>& y)
{
    double norma = 0;
    for (int i = 0; i < globalN; i++)
        norma += y[i] * y[i];
    return sqrt(norma);
}

double vector_multiplication(vector<double>& a, vector<double>& b)
{
    double s = 0;
    for (int i = 0; i < globalN; i++)
        s += a[i] * b[i];
    return s;
}

// Aa=b
vector<double> multiplication_matrix_on_vector(vector<double>& a, vector<double>& b)
{
    for (int i = 0; i < globalN; i++)
        b[i] = di[i] * a[i];

    for (int i = 1; i < globalN; i++)
    {
        int i0 = ig[i];
        int i1 = ig[i + 1];
        for (int j = 0; j < (i1 - i0); j++)
        {
            b[i] += ggl[i0 + j] * a[jg[i0 + j]]; // нижний треугольник
            b[jg[i0 + j]] += ggu[i0 + j] * a[i]; // верхний треугольник
        }
    }
    return b;
}

// d = a + b * c
vector<double> summ(vector<double>& a, double b, vector<double>& c, vector<double>& d)
{
    for (int i = 0; i < globalN; i++)
        d[i] = a[i] + b * c[i];
    return d;
}

// ЛОС с предобусловливанием
void LOS_D(vector<double> &q)
{
    z.resize(globalN);
    r.resize(globalN);
    p.resize(globalN);
    y.resize(globalN);
    t.resize(globalN);

    k_iter = 0;

    double alpha, betta;
    double pk_1_rk_1, pk_1_pk_1;
    vector<double> r0;
    r0.resize(globalN);
    multiplication_matrix_on_vector(q, y); //y = Aq
    summ(pr, -1, y, r);	//	r0 = f - Aq
    for (int i = 0; i < globalN; i++)	// r0 = (f - Aq)/sqrt(D)
    {
        r[i] /= sqrt(di[i]);
        r0[i] = r[i];
    }

    for (int i = 0; i < globalN; i++)	// z0 = r0/sqrt(D)
        z[i] = r[i] / sqrt(di[i]);
    multiplication_matrix_on_vector(z, p); // p = Az0
    for (int i = 0; i < globalN; i++)	// z0 = r0/sqrt(D)
        p[i] /= sqrt(di[i]);
    for (int k = 1; k < maxiter; k++)
    {
        pk_1_rk_1 = vector_multiplication(p, r); // (p_(k-1),r_(k-1))
        pk_1_pk_1 = vector_multiplication(p, p); // (p_(k-1),p_(k-1))
        alpha = pk_1_rk_1 / pk_1_pk_1; //	alpha_k = (p_(k-1),r_(k-1)) / (p_(k-1),p_(k-1))
        summ(q, alpha, z, q); // q_k = q_(k-1) + alpha_k * z_(k-1)
        summ(r, -alpha, p, r); // r_k = r_(k-1) - alpha_k * p_(k-1)
        for (int i = 0; i < globalN; i++)	// y = r_k/sqrt(D)
        {
            y[i] = r[i] / sqrt(di[i]);
        }
        multiplication_matrix_on_vector(y, t); // t = A * r_k / sqrt(D)
        for (int i = 0; i < globalN; i++)	// t = (A * r_k / sqrt(D) ) / sqrt(D)
        {
            t[i] /= sqrt(di[i]);
        }
        betta = -vector_multiplication(p, t) / pk_1_pk_1; // betta_k = (p_k-1,D_1*A*D_1*r_k) / (p_(k-1),p_(k-1))
        summ(y, betta, z, z);//	z_k = r_k / sqrt D + betta_k * z_(k-1)
        summ(t, betta, p, p);//p_k = D_1*A*D_1* r_k + betta_k * p_(k-1)

        if (vector_multiplication(r, r) < e) // (r_k, r_k) < e
        {
            k_iter = k;
            return;
        }
    }
    k_iter = maxiter;
}

void output(double x0, double x1, double hx, double y0, double y1, double hy, double curT, vector<double>& q)
{
    for(double y = y0; y <= y1; y += hy)
    {
        for (double x = x0; x <= x1; x += hx)
        {
            int areaNum = numberAreaForPoint(x, y);
            if (areaNum != -1)
            {
                double valueFind = valueFuncAtPoint(x, y, q);
                double valueReal = functionU(x, y, curT, areaNum);
                printf_s("x:%f y:%f t:%f %.15lf %.15lf %.15lf\n", x, y, curT, valueReal, valueFind, abs(valueFind - valueReal));
                //printf_s("%.15lf\n", valueFind);
                //printf_s("%.15lf\n", valueReal);
            }
            else
            {
                printf_s("Фиктивная область\n");
            }
        }
    }
}

void findQ(steps schemeSteps, double curT)
{
    makeGlobalMatrixAndVector(curT, 4, schemeSteps);
    boundaryConditions(curT);
    LOS_D(q_j);
    clearMatrix();
}

// Нахождение q(j-1) через 3 слойную
void findQ1(steps schemeSteps, double curT)
{
    makeGlobalMatrixAndVector(curT, 3, schemeSteps);
    boundaryConditions(curT);
    LOS_D(q_j1);
    clearMatrix();
}

// Нахождение q(j-2) через 2 слойную
void findQ2(steps schemeSteps, double curT)
{
    makeGlobalMatrixAndVector(curT, 2, schemeSteps);
    boundaryConditions(curT);
    LOS_D(q_j2);
    clearMatrix();
}

// Нахождение q(j-3) из начального приближения
void findQ3()
{
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            int globalNumber = i * Nx + j;
            q_j3[globalNumber] = u0(X[j], Y[i], numberAreaForPoint(X[j], Y[i]));
        }
    }
}

void initApproximation()
{
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            int globalNumber = i * Nx + j;
            int areaNum = numberAreaForPoint(X[j], Y[i]);
            q_j3[globalNumber] = uAtTime(X[j], Y[i], T[0], areaNum);
            q_j2[globalNumber] = uAtTime(X[j], Y[i], T[1], areaNum);
            q_j1[globalNumber] = uAtTime(X[j], Y[i], T[2], areaNum);
        }
    }
}

void calcSolution()
{
    // Структура для хранения шагов
    steps schemeSteps;
    // Для тестирования только 4 слойной схемы
    bool isTest = 0;

    if (!isTest)
    {
        for (int ti = 0; ti < Nt && ti < 3; ti++)
        {
            // Нахождение q(j-3) из начального приближения
            if (ti == 0)
            {
                findQ3();
                printf_s("\nt:%d = %.5lf, k_iter = %d\n", ti, T[ti], k_iter);
                output(1, 5, 1, 1, 5, 1, T[ti], q_j3);
            }
            // Нахождение q(j-2) через 2 слойную
            else if (ti == 1)
            {
                schemeSteps.deltaT = T[ti] - T[ti - 1];
                findQ2(schemeSteps, T[ti]);
                printf_s("\nt:%d = %.5lf, k_iter = %d\n", ti, T[ti], k_iter);
                output(1, 5, 1, 1, 5, 1, T[ti], q_j2);
            }
            // Нахождение q(j-1) через 3 слойную
            else
            {
                schemeSteps.deltaT = T[ti] - T[ti - 2];
                schemeSteps.deltaT0 = T[ti] - T[ti - 1];
                schemeSteps.deltaT1 = T[ti - 1] - T[ti - 2];
                findQ1(schemeSteps, T[ti]);
                printf_s("\nt:%d = %.5lf, k_iter = %d\n", ti, T[ti], k_iter);
                output(1, 5, 1, 1, 5, 1, T[ti], q_j1);
            }
        }
    }
    else
    {
        initApproximation();
        printf_s("\nt:%d = %.5lf, k_iter = %d\n", 0, T[0], k_iter);
        output(1, 5, 1, 1, 5, 1, T[0], q_j3);
        printf_s("\nt:%d = %.5lf, k_iter = %d\n", 1, T[1], k_iter);
        output(1, 5, 1, 1, 5, 1, T[1], q_j2);
        printf_s("\nt:%d = %.5lf, k_iter = %d\n", 2, T[2], k_iter);
        output(1, 5, 1, 1, 5, 1, T[2], q_j1);
    }

    for (int ti = 3; ti < Nt; ti++)
    {
        schemeSteps.deltaT = T[ti] - T[ti - 3];
        schemeSteps.deltaT0 = T[ti] - T[ti - 1];
        schemeSteps.deltaT1 = T[ti - 1] - T[ti - 2];
        schemeSteps.deltaT2 = T[ti - 2] - T[ti - 3];
        findQ(schemeSteps, T[ti]);
        printf_s("\nt:%d = %.5lf, k_iter = %d\n", ti, T[ti], k_iter);
        output(0, 10, 2.5, 0, 10, 2.5, T[ti], q_j);

        swap(q_j2, q_j3);
        swap(q_j1, q_j2);
        swap(q_j, q_j1);

        clearVector(q_j);
    }

    double norm = normL2(T[Nt - 1], q_j1);
    printf_s("\n t = %.5lf, norm = %.15lf\n", T[Nt - 1], norm);
}

int main()
{
    setlocale(LC_ALL, "");
    inputGrid();
    inputPartitions();
    inputTimeLayers();
    inputKraev();
    generatePortrait();
    calcSolution();
}
