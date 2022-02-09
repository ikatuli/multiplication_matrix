#include <fstream>
#include <iostream>
#include <string>

std::pair<int,int> mn (std::ifstream &in)
{
	std::string m;
	std::string n;
	getline(in, m, ' ');
	getline(in, n, '\n');
	return std::pair<int,int>(std::stoi(m),std::stoi(n));
}

void input(std::ifstream &in,int M[],int m,int n)
{
	std::string item;
	for (int i=0;i<m;i++)
		{
			for (int j=0;j<n;j++)
			{
				if (j+1!=n)	
					getline(in, item, ' '); 
				else
				   getline(	in, item, '\n');
				M[(i*n+j)]=std::stoi(item);
			}
		}
}

void print(int M[], int m,int n)
{
	for (int i=0;i<m;i++)
		{
			for (int j=0;j<n;j++)
			{
				std::cout << M[(i*n+j)] << " ";
        }
        std::cout << std::endl;
    }

}

int main(int argc, char* argv[])
{
	std::string item;
	std::ifstream in("input.txt");

	auto [m,n]=mn(in); // Получаем размерность матрицы A из файла

	int A[m*n];

	if (in.is_open()){
		input(in,A,m,n);
    }
	print(A,m,n);

	auto [n1,k]=mn(in); // Получение матрицы B
	int B[n*k];

	if (in.is_open()){
		input(in,B,n,k);
    }

	print(B,n,k);

    in.close();

	return 0;
}
