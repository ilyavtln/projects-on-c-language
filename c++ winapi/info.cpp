#include <windows.h>
#include <stdio.h>

extern "C" _declspec(dllexport) int Information(char* InfoString)
{
    // Получаем высоту экрана в пикселях
    int height = GetSystemMetrics(SM_CYSCREEN);

    // Записываем в буфер
    sprintf(InfoString, "Высота экрана в пикселях: %d px", height);

    return EXIT_SUCCESS;
}