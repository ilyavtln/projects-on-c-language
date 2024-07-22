#include <iostream>
#include <fstream>
#include <cmath>

#define strin "%lf"
#define strout "%.15lf "

typedef double type;
typedef double summa;

using namespace std;

int n; //Размер 
int m; //Элементов в al и au
int* ia; //Индексный массив
type* au; //Верхний треугольник
type* al; //Нижний треугольник
type* di; ///Диагональ
type* f; //Исходный вектор
type* x; //Вектор x
type* y; //Вектор y 
type** M; //Плотная матрицы

void input()
{
	ifstream size("size.txt");
	size >> n;
	size.close();

	ifstream matrix("matrix.txt");
	//Индексный массив
	ia = new int[n + 1];
	for (int i = 0; i < n + 1; i++) matrix >> ia[i];
	if (ia[0] == 1) //Для чтения файлов программ на Фортране
		for (int i = 0; i <= n; i++) ia[i]--;
    m = ia[n]; 
	//Нижний треугольник
	al = new type[m];
	for (int i = 0; i < m; i++) matrix >> al[i];
	//Диагональ
	di = new type[n];
	for (int i = 0; i < n; i++) matrix >> di[i];
	//Верхний треугольник
	au = new type[m];
	for (int i = 0; i < m; i++) matrix >> au[i];
	matrix.close();

	ifstream vector("vector.txt");
	f = new type[n];
	for (int i = 0; i < n; i++) vector >> f[i];
	vector.close();

	x = f;
	y = f;
}

void minusK(int ex)
{
    di[0] -= pow(10, -ex);
    f[0] -= pow(10, -ex);
}

void CalcLU()
{
    for (int i = 0; i < n; i++)
    {
        int i0 = ia[i];
        int i1 = ia[i + 1];
        type sum = 0;
        for (int k = i0, j = i - (i1 - i0); k < i1; k++, j++)
        {
            int j0 = ia[j];
            int j1 = ia[j + 1];

            int ik = i0;
            int kj = j0;

            int kij = k - i0 - (j1 - j0);
            if (kij < 0)
                kj -= kij;
            else
                ik += kij;

            type suml = 0;
            type sumu = 0;
            for (ik; ik < k; ik++, kj++)
            {
                suml += al[ik] * au[kj];
                sumu += au[ik] * al[kj];
            }
            al[k] -= suml;
            au[k] -= sumu;
            au[k] /= di[j];
            sum += al[k] * au[k];
        }
        di[i] -= sum;
    }
}

void CalcY() //Ly=f
{
    for (int i = 0; i < n; i++)
    {
        summa sum = 0;
        int i0 = ia[i];
        int i1 = ia[i + 1];
        for (int j = i - (i1 - i0), k = i0; k < i1; j++, k++)
            sum += al[k] * y[j];
        y[i] = (f[i] - sum) / di[i];
    }
}

void CalcX() //Ux = y
{
    for (int i = n - 1; i >= 0; i--)
    {
        summa sum = 0;
        int i0 = ia[i];
        int i1 = ia[i + 1];
        type xi = x[i] = y[i];
        for (int j = i - (i1 - i0), k = i0; k < i1; j++, k++)
        {
            sum = au[k] * xi;
            y[j] -= sum;
        }
    }
}

void LU()
{
    CalcLU();
    CalcY();
    CalcX();
}

void makeGilbert(int k)
{
    type** G = new type* [k];
    for (int i = 0; i < k; i++)
    {
        G[i] = new type [k];
        for (int j = 0; j < k; j++)
        {
            G[i][j] = 1.0 / (i + j + 1);
        }
    }
    //Вывод размера матрицы
    FILE* size = NULL;
    fopen_s(&size, "size.txt", "w");
    fprintf_s(size, "%d", k);
    fclose(size);

    //Вывод матрицы Гильберта
    FILE* matrix = NULL;
    fopen_s(&matrix, "matrix.txt", "w");
    //Заполнение индексного массива
    fprintf_s(matrix, "%d %d ", 0, 0);
    int temp = 0;
    for (int i = 1; i < k; i++)
    {
        temp += i;
        fprintf_s(matrix, "%d ", temp);
    }
    fprintf_s(matrix, "\n", NULL);
    //Заполнение массива al
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < i; j++)
        {
            fprintf_s(matrix, strout, G[i][j]);
        }
    }
    fprintf_s(matrix, "\n", NULL);
    //Заполнение диагонали
    for (int i = 0; i < k; i++)
        fprintf_s(matrix, strout, G[i][i]);
    fprintf_s(matrix, "\n", NULL);
    //Заполнение массива au
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < i; j++)
        {
            fprintf_s(matrix, strout, G[i][j]);
        }
    }
    fclose(matrix);
    //Вывод вектора в файл
    FILE* vector = NULL;
    fopen_s(&vector, "vector.txt", "w");
    type* r = new type[k];
    for (int i = 0; i < k; i++)
    {
        summa sum = 0;
        for (int j = 0; j < k; j++)
            sum += G[i][j] * (j + 1);
        r[i] = sum;
        fprintf_s(vector, strout, r[i]);
    }
    fclose(vector);
    //Очистка памяти
    delete[] r;
    for (int i = 0; i < k; i++)
        delete G[i];
    delete[] G;
}

void toPlot() //Создание матрицы в плотном формате
{
    ofstream plotA("plotA.txt");
    M = new type * [n];
    for (int i = 0; i < n; i++)
        M[i] = new type[n];

    for (int i = 0; i < n; i++)
        M[i][i] = di[i];

    for (int i = 1; i < n; i++)
    {
        int i0 = ia[i];
        int i1 = ia[i + 1];
        for (int j = 0; j < i - (i1 - i0); j++)
            M[i][j] = 0.0;
        for (int j = i - (i1 - i0); j < i; j++)
            M[i][j] = al[ia[i] + j - i + (i1 - i0)];
        for (int j = 0; j < i - (i1 - i0); j++)
            M[j][i] = 0.0;
        for (int j = i - (i1 - i0); j < i; j++)
            M[j][i] = au[ia[i] + j - i + (i1 - i0)];
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            plotA << M[i][j] << " ";
        }
        plotA << "\n";
    }
    plotA.close();
}

void Gauss() // с выбором ведущего элемента
{
    for (int i = 0; i < n - 1; i++) //прямой ход: приводим к верхнетреугольному виду
    {
        type max = abs(M[i][i]);
        int ind = i;
        for (int j = i + 1; j < n; j++) // поиск наибольшего по модулю в столбце
        {
            if (abs(M[j][i]) > max)
            {
                max = abs(M[j][i]);
                ind = j;
            }
        }
        if (max > 0)
        {
            if (ind != i) // замена строк, если ведущий элемент не на диагонали
            {
                swap(M[i], M[ind]);
                swap(f[i], f[ind]);
            }
        }
        if (M[i][i] != 0) // проверка что на диагонали не 0
        {
            for (int j = i + 1; j < n; j++)
            {
                type b = M[j][i] / M[i][i];
                M[j][i] = 0;
                f[j] -= b * f[i];
                for (int m = i + 1; m < n; m++)
                    M[j][m] -= b * M[i][m];
            }
        }
    }

    x = new type[n];

    for (int i = n - 1; i >= 0; i--) // поиск вектора x
    {
        type sum = 0;
        for (int j = i + 1; j < n; j++)
            sum += M[i][j] * x[j];
        x[i] = (f[i] - sum) / M[i][i];
    }
}

void outputX()
{
    FILE* res;
    fopen_s(&res, "output.txt", "w");
    for (int i = 0; i < n; i++)
        fprintf_s(res, strout, x[i]);
    fclose(res);
}

void output()
{
    FILE* lual;
    FILE* luau;
    FILE* ludi;
    fopen_s(&lual, "matrixL.txt", "w");
    fopen_s(&luau, "matrixU.txt", "w");
    fopen_s(&ludi, "diagL.txt", "w");
    for (int i = 0; i < n; i++)
        fprintf_s(ludi, strout, di[i]);
    for (int i = 0; i < m; i++)
    {
        fprintf_s(lual, strout, al[i]);
        fprintf_s(luau, strout, au[i]);
    }
    fclose(lual);
    fclose(luau);
    fclose(ludi);
}


void printMenu()
{
	cout << "Выберите действие:" << endl;
	cout << "1. Решить, используя LU разложение" << endl;
    cout << "2. Исследование на обусловленность" << endl;
	cout << "3. Матрица Гильберта, LU разложение" << endl;
	cout << "4. Решить, используя метод Гаусса" << endl;

	int answer;
	cin >> answer;
	switch (answer)
	{
		case 1:
			input();
            minusK(3);
            LU();
            outputX();
            output();
			break;
        case 2:
            input();
            minusK(3); //число к
            LU();
            outputX();
            break;
		case 3:
            makeGilbert(12); //Создание матрицы Гильберта размера k
            input();
            LU();
            outputX();
            output();
			break;
		case 4:
            input();
            minusK(0); //число к
            toPlot();
            Gauss();
            outputX();
			break;
	}
}

int main()
{
    setlocale(LC_ALL, "Russian");
    setlocale(LC_NUMERIC, "C");
	printMenu();
    return 0;
}
