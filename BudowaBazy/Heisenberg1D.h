#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <complex>
#include <clocale>
#include <armadillo>
#include <mpi.h>
#include <Windows.h>
#include <thread> //multithreading
#include <mutex>
#include <omp.h>
#include <execution>
//#include <atomic> // atomic variable type allows multithreading without mutex locking

using namespace std;
using namespace arma;

class Heisenberg1D{
private:
	int L; //chain length
	double Delta; //anisotropic term
	double J; //exchange integral
	int S_up; //number of up-spins

	int N; //number of states - if using blocks its binomial(n,i) - i-th block, otherwise N=2^L
	vector<int> base_vector;
	Col<int> mapping; //generates the mapping of the base_vector number to the number in the block
	Col<int> mapping_inv;
	
	cx_vec phi_T;

	mat H;
	mat H_magnetic;
	mat eigenvectors;
	vec eigenvalues;

	mat Krylov_space;
	vec Lanczos_eigenVal;
	mat Lanczos_eigenVec;;
	mat H_Lanczos;
	//cx_mat unitary_evolution;

	vector<vector<string>> results; // to store data until end of multithreading

public:
	Heisenberg1D(int L, double J, double Delta); //Constructor for total Hamiltonian
	Heisenberg1D(int L, int S_up, double J, double Delta); //Constructor for subblock Hamiltonian
	~Heisenberg1D();

	vector<int> first_vec;
	mat get_hamil();
	cx_vec coeff;
	vector<double> Lanczos_gs;

	void Hamiltonian();
	void Magnetic_perturbance(double q, double h);
	void Diagonalization();

	void Matrix_vector_multiply(mat Matrix, vec input_vector, vec output_vector);
	void Lanczos_Diagonalization(int m);
	void Lanczos_gs_energy();
	void Lanczos_convergence(int m, stringstream& ss);
	void print_eigenvalues(ofstream& fileE);

	//Time evolution, magnetization, correlation etc.
		cx_vec Time_evolution(cx_vec wavefuncton, double dt);
		cx_vec Time_evolution_stationary(cx_vec wavefuncton, double time);
		void generate_coefficients_vector(cx_vec domain_wall_vec);
		void Magnetization_evolution(string del);
		void Magnetization_1thread(const int beg, const int end, vector<double> time_vec, ofstream& file_evolution);
		//void Magnetization_evolution_total(string del);
		void print_evolution(string del);
		void Overlap(double time_end, double neel_temperature, string del);
		void Spin_Correlations(string del);

		void Quantum_fidelity(cx_vec wavefunction); //here change to have wavefunction as input vector
	//--------------------------------------

	// Linear responce
		void spinStructureFactor_lin_responce(double temperature);
		void spinStructureFactor_T_0();
		void spinStructureFactor_FFT(double temperature);
	//-----------------------

		double get_average_gap();
		void SSF_magnetic_field(double h);
		void SSF_Lanczos();
		void Build_Lanczos_Hamil(cx_vec initial_vec, int M);
		void Lanczos_Heat_capacity_Randomized(int M);


	// Tools
		void generate_mapping_subblock();
		void generate_mapping_total();

		void print_Hamiltonian_subblocks(wofstream &file);
		void print_Hamiltonian_total(wofstream& file);

		vector<int> int_to_basevector(int idx); //converges int to eigenvector with S_up up-spins
		int basevector_to_int(vector<int> vec); //converges binar vector to index in given subblock
	//-----------------------------

	// Trivial quantities
		Col<int> check_eigenvalues();
		double get_energy_gap();
		double get_domain_wall_state_temp(double T_step, double T_end, string del);
		void Heat_Capacity(double step, double T_0, double T_end, ofstream& savefile);
	//----------------------------
};

vector<int> int_to_binary(int idx, int L); //converges int to binary code of length N
int binary_to_int(vector<int> vec); //converges vector with binary code to decimal system

long long int factorial(int n);
long long int Binomial(int n, int k); //Binomial coefficient n over k

void Compute_subblocks(int L, double J, double Delta, string del); //main program to execude
void Compute_total(int L, double J, double Delta, string del); //main program to execude



