#pragma warning(disable: 4996)
#pragma comment (lib,"Ws2_32.lib")
#include <iostream>
#include <string>
#include <WinSock2.h>

using namespace std;

int main()
{
    //�������� ����������
    WSADATA wsaData;
    //������ ���������� WinSock
    WORD sockVer = MAKEWORD(2, 2);
    //������������� ��������
    int wsastarup = WSAStartup(sockVer, &wsaData);
    if (wsastarup != 0)
    {
        cout << "WSAStartup error: " << wsastarup << endl;
        return wsastarup;
    }

    // �������� ��������������� ������
    SOCKET ServerSocket = INVALID_SOCKET;
    ServerSocket = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);
    if (ServerSocket == INVALID_SOCKET)
    {
        cout << "socket creation error: " << WSAGetLastError() << endl;
        WSACleanup();
        return 1;
    }

    //����� ������
    SOCKADDR_IN sin;
    //��������� ����������
    sin.sin_family = PF_INET;
    //port
    sin.sin_port = htons(2003);
    //ip �����
    sin.sin_addr.s_addr = INADDR_ANY;

    // �������� ������ ������
    int retVal = bind(ServerSocket, (LPSOCKADDR)&sin, sizeof(sin));
    if (retVal == SOCKET_ERROR)
    {
        cout << "socket binding error: " << WSAGetLastError() << endl;
        closesocket(ServerSocket);
        WSACleanup();
        return 1;
    }

    cout << "Server started at " << inet_ntoa(sin.sin_addr) << ", port " << htons(sin.sin_port) << endl;

    do
    {
        //������������� ����� �����
        retVal = listen(ServerSocket, SOMAXCONN);
        if (retVal == SOCKET_ERROR)
        {
            cout << "socket listening error: " << WSAGetLastError() << endl;
            closesocket(ServerSocket);
            WSACleanup();
            return 1;
        }
        //���� �������
        SOCKET clientSocket;
        SOCKADDR_IN from;
        int fromlen = sizeof(from);
        int maxLen = 1024;
        char* buf = new char[maxLen] {};
        string output = "";
        // ����� ����������
        clientSocket = accept(ServerSocket, (SOCKADDR*)&from, &fromlen);
        if (clientSocket == INVALID_SOCKET)
        {
            cout << "connection accept error: " << WSAGetLastError() << endl;
            closesocket(ServerSocket);
            WSACleanup();
            return 1;
        }
        cout << "New connection accepted from " << inet_ntoa(from.sin_addr) << ", port " << htons(from.sin_port) << endl;
        // ��������� ������
        retVal = recv(clientSocket, buf, maxLen, 0);
        if (retVal == SOCKET_ERROR)
        {
            cout << "data receiving error: " << WSAGetLastError() << endl;
            closesocket(clientSocket);
            WSACleanup();
            return 1;
        }
        else if (retVal > 0)
        {
            //��������� �������� ip � port �������
            string ip = inet_ntoa(sin.sin_addr);
            string port = to_string(htons(sin.sin_port));
            for (int i = 0; i < retVal; i++)
            {
                output += buf[i];
                // ��������� ip ������� � ���� � ����� ������� �����������
                if (buf[i] == '.' || buf[i] == '!' || buf[i] == '?')
                {
                    output += "(" + ip + ":" + port + ")";
                }
            }
            // ���������� ���������
            int sendResult = send(clientSocket, output.c_str(), output.length(), 0);
            if (sendResult == SOCKET_ERROR) cout << "send error: " << WSAGetLastError() << endl;
            closesocket(clientSocket);
        }
        else cout << "connection closed" << endl;
    } while (retVal > 0);
    WSACleanup();
    return 0;
}
