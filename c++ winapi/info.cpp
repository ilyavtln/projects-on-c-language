#include <windows.h>
#include <stdio.h>

extern "C" _declspec(dllexport) int Information(char* InfoString)
{
    // �������� ������ ������ � ��������
    int height = GetSystemMetrics(SM_CYSCREEN);

    // ���������� � �����
    sprintf(InfoString, "������ ������ � ��������: %d px", height);

    return EXIT_SUCCESS;
}