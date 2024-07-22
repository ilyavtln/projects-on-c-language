#pragma comment (lib,"Ws2_32.lib")
#pragma warning(disable: 4996)
#include <winsock2.h>
#include <string>
#include <iostream> 

using namespace std;

int main()
{
	//�������� ����������
	WSADATA wsaData;
	//������ ���������� WinSock
	WORD sockVer = MAKEWORD(2, 2);
	int retVal = 0;
	WSAStartup(sockVer, (LPWSADATA)&wsaData);
	LPHOSTENT hostEnt;
	hostEnt = gethostbyname("localhost");
	if (!hostEnt)
	{
		printf("Unable to collect gethostbyname\n");
		WSACleanup();
		system("pause");
		return 1;
	}
	//������� �����
	SOCKET clientSock = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
	if (clientSock == SOCKET_ERROR)
	{
		printf("Unable to create socket\n");
		WSACleanup();
		system("pause");
		return 1;
	}
	string ip;
	int port;
	cout << "ip: ";
	cin >> ip;
	cout << "port: ";
	cin >> port;
	cin.ignore();
	SOCKADDR_IN serverInfo; //��������� ����� ��� ������ � ����������� 
	serverInfo.sin_family = AF_INET;
	serverInfo.sin_addr.S_un.S_addr = inet_addr(ip.c_str());
	serverInfo.sin_port = htons(port);
	//�������� ������������� � ������� �� ip � port
	retVal = connect(clientSock, (LPSOCKADDR)&serverInfo, sizeof(serverInfo));
	if (retVal == SOCKET_ERROR)
	{
		printf("Unable to connect\n");
		int ceck = WSAGetLastError();
		WSACleanup();
		system("pause");
		return 1;
	};
	cout << "Connection made sucessfully" << endl;
	char pBuf[256];
	cout << "message: ";
	gets_s(pBuf);
	//�������� ������ �� ������
	retVal = send(clientSock, pBuf, strlen(pBuf), 0);
	if (retVal == SOCKET_ERROR)
	{
		printf("Unable to send\n");
		WSACleanup();
		system("pause");
		return 1;
	}

	char szResponse[256]{};
	//�������� �������� ����� �� �������
	retVal = recv(clientSock, szResponse, 256, 0); //��������� ������ �� ������
	if (retVal == SOCKET_ERROR)
	{
		printf("Unable to recv\n");
		WSACleanup();
		system("pause");
		return 1;
	}
	char* Resp = szResponse;
	printf("%s\n", Resp);
	closesocket(clientSock);
	WSACleanup();
	system("pause");
	return 0;
}
