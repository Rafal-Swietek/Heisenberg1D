#include "Heisenberg1D.h"

using namespace std;

void test(int L, double J, double Delta, ofstream& fileEigenvalues, string del) {

	Heisenberg1D Object(L, L / 2, J, Delta);
	Object.Hamiltonian();
	//Print Hamiltonian
		/*//create a locale that has the ctype category  copied from the "en_US.UTF-8"
			locale loc = locale(locale(), "en_US.UTF8", locale::ctype);
			wofstream file;
			file.open("Hamiltonian_subblocks_L=" + to_string(L) + "_Delta=" + del + ".txt");
			file.imbue(loc); // adding the locale to the stream
		//-------------------------------------------------------
		Object.print_Hamiltonian_subblocks(file);
		file.close();*/
	//----------------
	cout << "\tHamiltonian L=" << L << " created!\n";

	//Diagonalization
	Object.Diagonalization();
	cout << "\tHamiltonian L=" << L << " diagonalized!!" << endl;
	//----------------

	/*Object.print_eigenvalues(fileEigenvalues);
		cout << "Eigenvalues correct?:\n";
		Col<int> check = Object.check_eigenvalues();
		cout << check;
		cout << "\nOh yeah";*/
	
	Object.~Heisenberg1D();
}

void matrix_vector_mult(mat A, vec B, vec C, int dimension) {
	C.zeros();
#pragma omp parallel for shared(A,B,C,dimension) collapse(2)
	for (int j = 0; j < dimension; j++) {
		double temp = 0;
#pragma omp parallel for shared(A,B,dimension, j) reduction(+: temp)
		for (int k = 0; k < dimension; k++) {
			temp += A(j, k) * B(k);
		}
		C(j) = temp;
	}
}

void main(int argc, char* argv[]) {

	time_t t_start, t_end;
	time(&t_start);

	int L = 10; //chain length
	double J = 1; //exchange integral
	double Delta = 1.0; // anisotropic term
	
	/*int N = pow(2, L);

	mat A = randu(N, N);
	vec B = randu(N);
	vec C = randu(N);
	for (int j = 0; j < N; j++) {
		C = A * B;
		//matrix_vector_mult(A, B, C, N);
	}*/


	TCHAR Npath[MAX_PATH];
	GetCurrentDirectory(MAX_PATH, Npath);
	CreateDirectory("./Results/", NULL);
	SetCurrentDirectory("./Results/");
	/*if (J < 0) {
		CreateDirectory("./Results/FM/", NULL);
		SetCurrentDirectory("./Results/FM/");
	}
	else{
		CreateDirectory("./Results/AFM/", NULL);
		SetCurrentDirectory("./Results/AFM/");
	}*/
		//text values of parameters
			stringstream s1;
			s1 << fixed << setprecision(0) << Delta;
			string del = s1.str();
		//--------------
		//test(L, J, Delta, fileEigenvalues, del);
		//Compute_subblocks(L, J, Delta, del);

		Compute_total(L, J, Delta, del);
		//cout << "Delta=" + del + "complete" << endl;

		/*Delta = 1.0;
		//text values of parameters
		s1.str(string{});
		s1 << fixed << setprecision(0) << Delta;
		del = s1.str();
		//--------------
		//test(L, J, Delta, fileEigenvalues, del);
		Compute_subblocks(L, J, Delta, del);
		cout << "Delta=" + del + "complete" << endl;

		Delta = 2.0;
		//text values of parameters
			s1.str(string{});
			s1 << fixed << setprecision(0) << Delta;
			del = s1.str();
		//--------------
		//test(L, J, Delta, fileEigenvalues, del);
		Compute_subblocks(L, J, Delta, del);
		cout << "Delta=" + del + "complete" << endl;*/

		//Separating execude code (Compute) into threads 
		//- each threat generates, diagonalizes and process a Hamiltonian with L atoms
			/*vector<thread> thread; //maximal size of 16
			thread.reserve(8); // L = 2 - 16 every 2
			for (L = 2; L <= 16; L+=2) {
				if(L % 2 == 0) thread.emplace_back(Compute, L, J, Delta, ref(fileGapeven), ref(fileEigenvalues), del); //add thread to vector of threads
				else thread.emplace_back(Compute, L, J, Delta, ref(fileGapodd), ref(fileEigenvalues), del);
				cout << "Thread with L=" << L << " has joined the party" << endl;
			}
			for (auto& t : thread) t.join(); //Wait for all threads to finish

			for (auto& t : thread) t.~thread(); //destroy all threads - no memory leackage

			//Last execude exceeeds memory -> separetely calculated
			// calculating as thread -> 2GB RAM to generate H, but separately 1GB
			//L = 13;
			//Compute(L, J, Delta, ref(fileGapodd), ref(fileEigenvalues), del);
			thread.clear();*/
		//------------------------------------

	// SIMPLE Hamiltonian build and print with arrow coded basevectors
	/*	//create a locale that has the ctype category  copied from the "en_US.UTF-8"
			locale loc = locale(locale(), "en_US.UTF8", locale::ctype);
			wofstream file, file2;
			file.open("Hamiltonian_subblocks_L=" + to_string(L) + "_Delta="+del+".txt");
			file.imbue(loc); // adding the locale to the stream
		//-------------------------------------------------------
		// Hamiltonian in subblocks for given coefficents
			vector<double> eigenvalues;
			for (int S_up = 0; S_up <= L; S_up++) {
				Heisenberg1D Object(L, S_up, J, Delta); // S_up depicts the number of up-spins
				Object.Hamiltonian_subblocks();
				Object.print_Hamiltonian_subblocks(file); //printing with representation of base vectors (as uparrow and downarrow)
				 //file << object.get_hamil(); //normal printing
			}
			file.close();
			fileC.close();
			fileGap.close();
		//----------------------------------------------------

		//Total Hamiltonian with given coefficients
			file2.open("Hamiltonian_total_L=" + to_string(L) + "_Delta=" + del + ".txt");
			file2.imbue(loc); // add the locale to the stream
			Heisenberg1D Object(L, J, Delta);
			Object.Hamiltonian_total();
			Object.print_Hamiltonian_total(file2); //printing with representation of base vectors (as uparrow and downarrow)
			//file2 << object.get_hamil(); //normal printing
			file2.close();
		//----------------------------------------------------*/
		
	SetCurrentDirectory(Npath);

	//Calculate execution time
		time(&t_end);
		int hours, minutes, seconds;
		hours = static_cast<int>((t_end - t_start) / 3600);
		minutes = static_cast<int>((t_end - t_start - 3600*hours) / 60);
		seconds = t_end - t_start - 3600 * hours - 60 * minutes;
		cout << "\nProgram executed in:\t" << fixed << hours << "h " << minutes << "min " << seconds << "sec" << setprecision(2) << endl;
}





/*
//--------------------------------------------------------------
//----------------------------COMMENTS--------------------------
//--------------------------------------------------------------
 - if L -> 30, change int L to long int L, because out of int bounds
 - expmat to matrix exponential: U=exp(-i/h*Ht), if symmetric: expmat_sym




*/