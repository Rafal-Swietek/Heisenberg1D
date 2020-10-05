#pragma once
#include "Heisenberg1D.h"

void DataDistribution(mat& matrix, vec& pProcRows, vec& initial_vec, int& Size, int& RowNum);
void parallel_matrix_vector_product(mat& matrix, int& N, int& M, vec& initial_vec, int& N_in, vec& output_vec, int& N_out, int argc, char* argv[]);
void ProcessInitialization(vec& pProcRows, vec& pProceResult, int& Size, int& RowNum);