#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <mpi.h>

typedef struct {
    MPI_Comm grid_comm;
	MPI_Comm row_comm;  
	MPI_Comm col_comm;
    int n_proc;         // Количество процессов
    int grid_dim;       // размер сетки = sqrt(n_proc) */
    int my_row;         
	int my_col;         
	int my_rank;        // Номер процесса
} GridInfo;


void grid_init(GridInfo *grid)
{
    int old_rank;
    int dimensions[2];
    int wrap_around[2];
    int coordinates[2];
    int free_coords[2];

    /*Получаем информацию перед наложением  cart_grid */
    MPI_Comm_size(MPI_COMM_WORLD, &(grid->n_proc)); //Сохраняем количество процессов
    MPI_Comm_rank(MPI_COMM_WORLD, &old_rank); // Номер процесса

    grid->grid_dim = (int)sqrt(grid->n_proc); // корень от количества процессов
											  //
    if (grid->grid_dim * grid->grid_dim != grid->n_proc) { //  Вызывает ошибку, если обратно не получается количество процессов.
        std::cout << "Количество процессов должно быть таково, что бы его корень был целым числом \n";
        exit(-1);
    }
    /* 	set the dimensions */
    dimensions[0] = dimensions[1] = grid->grid_dim;
    wrap_around[0] = wrap_around[1] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, 1, &(grid->grid_comm)); //Изменяем топологию.
    MPI_Comm_rank(grid->grid_comm, &(grid->my_rank));
    MPI_Cart_coords(grid->grid_comm, grid->my_rank, 2, coordinates);
    grid->my_row = coordinates[0];
    grid->my_col = coordinates[1];

    /* Коммуникаторы для строк */
    free_coords[0] = 0;
    free_coords[1] = 1;
    MPI_Cart_sub(grid->grid_comm, free_coords, &(grid->row_comm));

    /* Коммуникаторы для столбцов */
    free_coords[0] = 1;
    free_coords[1] = 0;
    MPI_Cart_sub(grid->grid_comm, free_coords, &(grid->col_comm));
}

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

void matrix_dot(int *A, int *B, int *C, int size)
{ //Умножения блоков
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            for (int k = 0; k < size; ++k) {
                C[i * size + j] += A[i * size + k] * B[k * size + j];
            }
        }
    }
}

void FoxAlgorithm(int *A, int *B, int *C, int size, GridInfo *grid)
{
    int *buff_A = (int*)calloc(size * size, sizeof(int));
    MPI_Status status;
    int root;
    int src = (grid->my_row + 1) % grid->grid_dim;
    int dst = (grid->my_row -1 + grid->grid_dim) % grid->grid_dim;

    for (int stage = 0; stage < grid->grid_dim; ++stage) {
        root = (grid->my_row + stage) % grid->grid_dim;
        if (root == grid->my_col) {
            MPI_Bcast(A, size * size, MPI_INT, root, grid->row_comm);
            matrix_dot(A, B, C, size);
        } else {
            MPI_Bcast(buff_A, size * size, MPI_INT, root, grid->row_comm);
            matrix_dot(buff_A, B, C, size);
        }
        MPI_Sendrecv_replace(B, size * size, MPI_INT, dst, 0, src, 0, grid->col_comm, &status);
    }
}

void square(int M[],int *pM,int m,int n, int size) //Копируем существующую матрицу в новую, квадратную.
{
	for (int i=0;i<size;i++)
	{
		for (int j=0;j<size;j++)
		{
			if ((i<m) and (j<n))
			{
				pM[(i*size+j)]=M[(i*n+j)];
			}
			else
			{
				pM[(i*size+j)]=0;
			}
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc < 3) {
		std::cout << "Укажите путь к файлам с матрицами\n";
        exit(-1);
		}

	int *pA, *pB, *pC;
	int *local_pA, *local_pB, *local_pC;
	int matrix_size;

	MPI_Init(&argc, &argv); // Инициализация MPI
	GridInfo grid;
    grid_init(&grid);

	//Для главного процесса
	if (grid.my_rank == 0) {
				
		std::ifstream in_A(argv[1]);
		
		auto [m,n]=mn(in_A); // Получаем размерность матрицы A из файла
						   
		int A[m*n];
		
		if (in_A.is_open()){
			input(in_A,A,m,n);
		}
		in_A.close();
		
		std::cout<<"Матрица A:\n";
		print(A,m,n);

		std::ifstream in_B(argv[2]);
		auto [n1,k]=mn(in_B); // Получение матрицы B

		int B[n*k];
		
		if (in_B.is_open()){
			input(in_B,B,n,k);
		}
		in_B.close();

		std::cout<<"Матрица B:\n";
		print(B,n,k);

		// Находи наибольшую сторону.

		//Не хочу ликовать algorithm. Напишу замену max
		if (n>m)
		{
			if (n>k) {matrix_size=n;}
			else {matrix_size=k;}
		}
		else 
		{
			if (m>k) {matrix_size=m;}
			else {matrix_size=k;}
		}

		while (matrix_size % grid.grid_dim != 0) //Квадратная матрица должна делится нацело на размер сетки
		{
			matrix_size+=1;
		}

		for (int thread=1; thread < grid.n_proc;thread++) //Отправляем размер матрицы в вспомогательные процессы
		{
			MPI_Send(&matrix_size,1,MPI_INT, thread, 0, MPI_COMM_WORLD);
		}

		//Превращаем матрицы в квадраты
		
		pA  = (int *)malloc(matrix_size * matrix_size * sizeof(int));
		square(A,pA,m,n,matrix_size);
		pB  = (int *)malloc(matrix_size * matrix_size * sizeof(int));
		square(B,pB,n,k,matrix_size);
		pC  = (int *)malloc(matrix_size * matrix_size * sizeof(int));

	}
	else //Получаем размер матрицы
	{
		MPI_Status status;
		MPI_Recv(&matrix_size,1, MPI_INT,MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	}


	int local_matrix_size = matrix_size / grid.grid_dim;
	local_pA  = (int *)malloc(local_matrix_size * local_matrix_size * sizeof(int));
	local_pB  = (int *)malloc(local_matrix_size * local_matrix_size * sizeof(int));
	local_pC  = (int *)malloc(local_matrix_size * local_matrix_size * sizeof(int));
	//Выделяем на них память

	MPI_Datatype blocktype, type;
    int array_size[2] = {matrix_size, matrix_size};
  	int subarray_sizes[2] = {local_matrix_size, local_matrix_size};
  	int array_start[2] = {0, 0};
  	MPI_Type_create_subarray(2, array_size, subarray_sizes, array_start,MPI_ORDER_C, MPI_INT, &blocktype);
    MPI_Type_create_resized(blocktype, 0, local_matrix_size * sizeof(int), &type);
  	MPI_Type_commit(&type);	
	

	int displs[grid.n_proc];
  	int sendcounts[grid.n_proc];
    if (grid.my_rank == 0) {
        for (int i = 0; i < grid.n_proc; ++i) {
            sendcounts[i] = 1;
        }
        int disp = 0;
        for (int i = 0; i < grid.grid_dim; ++i) {
          	for (int j = 0; j < grid.grid_dim; ++j) {
            		displs[i * grid.grid_dim + j] = disp;
            		disp += 1;
            }
            disp += (local_matrix_size - 1) * grid.grid_dim;
        }
    }

	MPI_Scatterv(pA, sendcounts, displs, type, local_pA, local_matrix_size * local_matrix_size, MPI_INT , 0, MPI_COMM_WORLD);
  	MPI_Scatterv(pB, sendcounts, displs, type, local_pB, local_matrix_size * local_matrix_size, MPI_INT , 0, MPI_COMM_WORLD);

	//Время
	double start_time, end_time;
    MPI_Barrier(grid.grid_comm);
    if (grid.my_rank == 0) {
        start_time = MPI_Wtime();
    }

	//Самое важное
	FoxAlgorithm(local_pA, local_pB, local_pC, local_matrix_size, &grid);

	MPI_Barrier(grid.grid_comm);
    if (grid.my_rank == 0) {
        end_time = MPI_Wtime() - start_time;
    }

	// Собираем куски матрицы по всем процессам
	
	MPI_Gatherv(local_pC, local_matrix_size*local_matrix_size, MPI_INT, pC, sendcounts, displs, type, 0, MPI_COMM_WORLD);

	if (grid.my_rank == 0) {
		std::cout << "Время: " << end_time << " Количество процессов:" <<  grid.n_proc << " Размер матрицы:" << matrix_size<<"\n Результат:\n";
		print(pC,matrix_size,matrix_size);
	}

	free (local_pA);// Выкидываем матрицы
	free (local_pB);
	free (local_pC);

	MPI_Finalize();
	return 0;
}
