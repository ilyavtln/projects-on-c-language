#include <windows.h>
#include <stdio.h>
#include <string.h>

#define BUF_SIZE 256

typedef int(*ImportFunction)(char*);
LRESULT CALLBACK WindowFunc(HWND hWind, UINT msg, WPARAM wParam, LPARAM lParam);

LPCSTR szClassName = "MyClass";
LPCSTR szTitle = "РГЗ по управлению ресурсами";
LPCSTR author = "Ватлин Илья, ПМ-13";

char Task[] = "Определить высоту экрана в пикселях:";
char Info[BUF_SIZE];

DWORD WINAPI ThreadFunc(void *)
{
	ImportFunction DLLInfo;
	HINSTANCE hinstLib = LoadLibrary(TEXT("info.dll"));
	if (hinstLib == NULL)
	{
		return EXIT_FAILURE;
	}

	DLLInfo = (ImportFunction)GetProcAddress(hinstLib, "Information");
	int flag = DLLInfo(Info);
	FreeLibrary(hinstLib);

	if (flag == NULL)
	{
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

LRESULT CALLBACK WindowFunc(HWND hWind, UINT msg, WPARAM wParam, LPARAM lParam)
{
	HANDLE hThread;
	DWORD IDThread;
	PAINTSTRUCT ps;
	HDC hDC;

	switch (msg)
	{
	case WM_CREATE:
		hThread = CreateThread(NULL, 0, ThreadFunc, NULL, 0, &IDThread);
		if (hThread != NULL)
		{
			MessageBox(hWind, author, szTitle, NULL);
			WaitForSingleObject(hThread, INFINITE);
			CloseHandle(hThread);
		}
		break;
	case WM_PAINT:
		hDC = BeginPaint(hWind, &ps);
		TextOut(hDC, 10, 10, Task, lstrlen(Task));
		TextOut(hDC, 10, 40, Info, lstrlen(Info));
		EndPaint(hWind, &ps);
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hWind, msg, wParam, lParam);
	}
	return EXIT_SUCCESS;
}

int WINAPI WinMain(HINSTANCE hThisInst, HINSTANCE hPrevInst, LPSTR str, int nWinMode)
{
    MSG msg;
    WNDCLASS wcl;
    HWND hWnd;

	wcl.hInstance = hThisInst;
	wcl.lpszClassName = szClassName;
	wcl.lpfnWndProc = WindowFunc;
	wcl.style = CS_HREDRAW | CS_VREDRAW;
	wcl.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	wcl.hCursor = LoadCursor(NULL, IDC_ARROW);
	wcl.lpszMenuName = NULL;
	wcl.cbClsExtra = 0;
	wcl.cbWndExtra = 0;
	wcl.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
	RegisterClass(&wcl); // Регистрация класса окна
	hWnd = CreateWindow(szClassName,
						szTitle, 
						WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN | WS_CLIPSIBLINGS, 
						0, 0, 500, 120,
						HWND_DESKTOP, 
						NULL, 
						hThisInst, 
						NULL);
	ShowWindow(hWnd, nWinMode);
	UpdateWindow(hWnd);
	while (GetMessage(&msg, NULL, 0, 0)) // Получение сообщения
	{
		TranslateMessage(&msg); // Преобразование виртуальных кодов клавиш в ASCII-значения
		DispatchMessage(&msg); // Посылка сообщения в нужную оконную процедуру
	}
	return msg.wParam;
}