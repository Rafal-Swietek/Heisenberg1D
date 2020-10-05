#include "parallelizing.h"

int ProcNum; // number of processes
int ProcRank; // rank of current process

void DataDistribution(mat& matrix, vec& pProcRows, vec& initial_vec, int& Size, int& RowNum) {
	int* pSendNum; // the number of elements sent to the process
	int* pSendInd; // the index of the first data element sent to the process
	int RestRows = Size; // Number of rows, that haven’t been distributed yet 
	//MPI_Bcast(initial_vec, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
void ProcessInitialization(vec& pProcRows, vec& pProcResult, int& Size, int& RowNum) {
	
	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int RestRows = Size;
	for (int i = 0; i < ProcRank; i++)
		RestRows = RestRows - RestRows / (ProcNum - i);
	RowNum = RestRows / (ProcNum - ProcRank);
	pProcRows = mat(RowNum, Size, fill::zeros);
	pProcResult = vec(RowNum);
}

void parallel_matrix_vector_product(mat& matrix, int& N, int& M, vec& initial_vec, int& N_in, vec& output_vec, int& N_out, int argc, char* argv[]) {
	vec pProcRows;
	mat pProcResult;
	int RowNum;
	if (N_in != N || N_out != M) throw "Dimensions do not agree!";
	else {
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
		MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
		DataDistribution(matrix, pProcRows, initial_vec, N, RowNum);

		MPI_Finalize();
	}
}