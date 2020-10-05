#include "Heisenberg1D.h"


mutex my_mutex;
complex<double> i = 1i;
double pi = 3.141592653;

//Destructor
Heisenberg1D::~Heisenberg1D() {
	H.~Mat();
	eigenvectors.~Mat();
	H_magnetic.~Mat();
	H_Lanczos.~Mat();
	Krylov_space.~Mat();
	Lanczos_eigenVal.~vec();
	Lanczos_gs.~vector();

	eigenvalues.~vec();
	mapping.~Col();
	mapping_inv.~Col();

	coeff.~cx_vec();
	phi_T.~cx_vec();

	base_vector.~vector();
	first_vec.~vector();
}
//-----------------------

//--------------------------------------------------------------------------------------------------
//----------------------------------------SUBBLOCKS HAMILTONIAN-------------------------------------
//--------------------------------------------------------------------------------------------------
// Constructor for Hamiltonian separated to blocks
Heisenberg1D::Heisenberg1D(int L, int S_up, double J, double Delta) {
	this->L = L; //number of sites
	this->S_up = S_up; //number of up-spins	
	this->N = Binomial(L, S_up); // number of states in given block
	this->Delta = Delta; this->J = J;

	this->base_vector = vector<int>(L); //spin on each site, for instance |0100> for L=4
	this->H = mat(N,N, fill::zeros); //hamiltonian

	this->first_vec = vector<int>(L, 0); //vector with S_up ones in the beggininng, vector filled with 0
	for (int j = 0; j < S_up; j++) { first_vec[L-j-1] = 1; } //setting last S_up bits to 1 (binary is inverse order of vector)
	this->mapping = Col<int>(std::pow(2, L)); mapping.fill(-1);
	this->mapping_inv = Col<int>(N); mapping_inv.fill(-1);
	generate_mapping_subblock();
}
//------------------------

//----------------------------------------------------------------------------------------------
//----------------------------------------TOTAL HAMILTONIAN-------------------------------------
//----------------------------------------------------------------------------------------------
// Constructor for total Hamiltonian (no subblocks)
Heisenberg1D::Heisenberg1D(int L, double J, double Delta) {
	this->L = L; //number of sites
	this->N = std::pow(2, L);
	this->Delta = Delta; this->J = J;

	this->base_vector = vector<int>(L); //spin on each site, for instance |0100> for L=4
	this->H = mat(N, N, fill::zeros); //hamiltonian

	int S_up = L / 2;
	this->first_vec = vector<int>(L, 0); //vector with S_up ones in the beggininng, vector filled with 0
	for (int j = 0; j < S_up; j++) { first_vec[L-j-1] = 1; } //setting last S_up bits to 1 (binary is inverse order of vector)

	this->mapping = Col<int>(std::pow(2, L)); mapping.fill(-1);
	this->mapping_inv = Col<int>(N); mapping_inv.fill(-1);
	generate_mapping_total();
}
//-------------------------

//----------------------------------------------------------------------------------------------
//----------------------------------------BUILD HAMILTONIAN-------------------------------------
//----------------------------------------------------------------------------------------------

// Generating Hamiltonian
void Heisenberg1D::Hamiltonian() {
	int s_i, s_j; //i=j, j=j+1
	int idx = 0; //indices equivalent to spin-filp due to kinetic term
	int PBC; //allows periodic boundary conditions
	for (int k = 0; k < N; k++) {
		base_vector = int_to_binary(mapping_inv(k), L);
		for (int j = 0; j < L - 1; j++) { //to delete PBC change to: j < L-1
			if (j == L - 1) PBC = 0;
			else PBC = j + 1;
			s_i = base_vector[j]; // i
			s_j = base_vector[PBC]; // i+1
			H(k,k) += J*Delta * (s_i - 0.5) * (s_j - 0.5); // "potential" term - diagonal part, spin s_z_i = eigenvector[i] - 1/2
			if (s_i == 0 && s_j == 1) { // S_i^+ S_i+1^-
				base_vector[j] = 1; //spin filp
				base_vector[PBC] = 0;
				//idx = basevector_to_int(base_vector); //getting index of the resulting eigenvector
				idx = mapping(binary_to_int(base_vector));
				base_vector[j] = 0; //in order to preserve eigenvector for next iteration
				base_vector[PBC] = 1;
				//idx = mapping( mapping_inv(k) - std::pow(2, j) ); //  we now what changed, therefore the decimal representation has to changed by +-2^j
				H(idx, k) += J / 2.;
			}
			if (s_i == 1 && s_j == 0) { // S_i^- S_i+1^+
				base_vector[j] = 0; //spin filp (inverse to the previous)
				base_vector[PBC] = 1;
				//idx = basevector_to_int(base_vector);
				idx = mapping(binary_to_int(base_vector));
				base_vector[j] = 1; //in order to preserve eigenvector for next iteration
				base_vector[PBC] = 0;
				//idx = mapping( mapping_inv(k) + std::pow(2, j) ); //  we now what changed, therefore the decimal representation has to changed by +-2^j
				H(idx, k) += J / 2.;
			}
		}
	}
}
//----------------------------------------------------\

// Generates the part of the hamiltonian with the magnetic field
void Heisenberg1D::Magnetic_perturbance(double q, double h) {
	this->H_magnetic = mat(N, N, fill::zeros);
	for (int j = 0; j < N; j++) {
		base_vector = int_to_binary(mapping_inv(j), L);
		for (int k = 0; k < L; k++) {
			H_magnetic(j, j) += h*cos(q * (double)k) * (base_vector[k] - 0.5);
		}
	}
}
//------------------------------------------------

// Printing Hamiltonian subblocks with S_up spin up
void Heisenberg1D::print_Hamiltonian_subblocks(wofstream& file) {
	double s_z = -L / 2. + S_up;
	file << "S^z_tot = " << s_z << "\n";
	file << "\t";
	for (int k = 0; k < N; k++) {
		base_vector = int_to_binary(mapping_inv(k), L);
		for (int j = 0; j < L; j++) {
			if (base_vector[j] == 0) { file << L'\u2193'; }
			if (base_vector[j] == 1) { file << L'\u2191'; }
		}
		file << "\t";
	}
	file << "\n";
	for (int k = 0; k < N; k++) {
		base_vector = int_to_binary(mapping_inv(k), L);
		for (int j = 0; j < L; j++) {
			if (base_vector[j] == 0) { file << L'\u2193'; }
			if (base_vector[j] == 1) { file << L'\u2191'; }
		}
		file << "\t";
		for (int j = 0; j < N; j++) {
			file << H(k,j) << "\t";
		}
		file << "\n";
	}
	file << " ________________________________________________\n\n";
}
//----------------------------------------------------

// Printing total Hamiltonian
void Heisenberg1D::print_Hamiltonian_total(wofstream& file) {
	file << "\t";
	for (int k = 0; k < N; k++) {
		base_vector = int_to_binary(k, L);
		for (int j = 0; j < L; j++) {
			if (base_vector[j] == 0) { file << L'\u2193'; }
			if (base_vector[j] == 1) { file << L'\u2191'; }
		}
		file << "\t";
	}
	file << "\n";
	for (int k = 0; k < N; k++) {
		base_vector = int_to_binary(k, L);
		for (int j = 0; j < L; j++) {
			if (base_vector[j] == 0) { file << L'\u2193'; }
			if (base_vector[j] == 1) { file << L'\u2191'; }
		}
		file << "\t";
		for (int j = 0; j < N; j++) {
			file << H(k, j) << "\t";
		}
		file << "\n";
	}
}
//----------------------------------------------------

//converges int to eigenvector with S_up up-spins - method: every permutation of the vector |00...0 1...11>
vector<int> Heisenberg1D::int_to_basevector(int idx) {
	vector<int> vec(first_vec);
	for (int j = 0; j < idx; j++) { //takes the idx-th permutation starting from initial state |00...0 11....1>
		next_permutation(vec.begin(), vec.end());
	}
	return vec;
}
//----------------------------------------------------

//converges binar vector to index in given subblock - generating the previous permutation until th efirst vector is obtained
int Heisenberg1D::basevector_to_int(vector<int> vec) {
	int val = 0;
	while (vec == first_vec) { //takes the previous permutation until it reaches the initial state |00...0 11....1>
		prev_permutation(vec.begin(), vec.end());
		val++;
	}
	return val;
}
//----------------------------------------------------

//generates the vector, which maps the base_vector index to the index in given subblock
void Heisenberg1D::generate_mapping_subblock() {
	int idx = 0; int k;
	do{ //first_vec used temporary, after loop returns to initial state
		k = binary_to_int(first_vec);
		mapping(k) = idx;
		try {
			mapping_inv(idx) = k;
		}
		catch(const out_of_range& e){
			std::cout << "Inverse mapping out of range";
		}
		idx++;
	} while (next_permutation(first_vec.begin(), first_vec.end()));
	//std::cout << mapping.t() << endl;
	//old version
	/*int sz; //sum of spins
	int idx = 0; //index in subblock
	for (int i = 0; i < mapping.size(); i++) {
		base_vector = int_to_binary(i,L); //here base_vector is just used temporary
		sz = 0;
		for (int j = 0; j < L; j++) {
			sz += base_vector(j);
		}
		if (sz == L / 2) {
			mapping(i) = idx; 
			idx++;
		}
	}*/
}
//----------------------------------------------------

//generates the vector, which maps for the total hamiltonian: map(i) = i
void Heisenberg1D::generate_mapping_total() {
	for (int j = 0; j < N; j++) { mapping(j) = j; mapping_inv(j) = j; }
}
//----------------------------------------------------

// Conversion of int to binary vector - using modulo operator
vector<int> int_to_binary(int idx, int L){
	vector<int> vec(L);
	int temp = idx ;
	for (int k = 0; k < L; k++) {
		vec[k] = temp % 2; 
		temp = static_cast<int>(temp / 2.);
	}
	return vec;
}
// Conversion of binary vector to int
int binary_to_int(vector<int> vec){
	int val = 0;
	for (int k = 0; k < vec.size(); k++) {
		val += vec[k] * std::pow(2, k);
	}
	return val;
}
//----------------------------------------------------

//----------------------------------------------------------------------------------------------
//--------------------------------------Diagonalization & exercises-----------------------------
//----------------------------------------------------------------------------------------------
void Heisenberg1D::Matrix_vector_multiply(mat Matrix, vec input_vector, vec output_vector) {

#pragma omp parallel for shared(Matrix, input_vector, output_vector)
	for (int k = 0; k < input_vector.size(); k++) {
		double temp = 0;
#pragma omp parallel for shared(Matrix, input_vector, output_vector) reduction(+: temp)
		for (int j = 0; j < input_vector.size(); j++) {
			temp += Matrix(k, j) * input_vector(k);
		}
		output_vector(k) = temp;
	}

}

//Diagonalizes the hamiltonian
void Heisenberg1D::Diagonalization() {
	this->eigenvalues = vec(N); //eigenvalues
	this->eigenvectors = mat(N, N, fill::zeros); //eigenvectors
	try {
		eig_sym(eigenvalues, eigenvectors, H); 
	}
	catch (const bad_alloc & e) {
		std::cout << "Memory exceeded" << e.what() << "\n";
		std::cout << H.size() * sizeof(H(0, 0)) << "\n";
		exit(123);
	}
	// Some shit to print most probable base_vector as ground state
		int tmp = 0;
		vec temp = eigenvectors.col(0); //ground state
		for (int j = 1; j < N; j++) //get biggest base state from ground state
			if (temp(j) > temp(tmp))
				tmp = j;
		vector<int> some_vec = int_to_binary(mapping_inv(tmp), L);
		std::cout << "\nBiggest state in ground state is:" << "\t\t";
		for (int j = 0; j < L; j++) std::cout << some_vec[j];
		std::cout << endl;
			/*locale loc = locale(locale(), "en_US.UTF8", locale::ctype);
			wcout.imbue(loc);
		wcout << "\nBiggest state in ground state is:" << "\t\t";
		SetConsoleOutputCP(CP_UTF8);
		for (int j = 0; j < L; j++) {
			if(some_vec[j] == 0) wcout << L"\x2193";
			else  wcout << L"\x2191";
		}
		wcout << "\n" << endl;*/
		temp.~vec(); some_vec.~vector();
		//cout << "Hamil:\n" << H << endl;
		//cout << "Energies:\n" << eigenvalues << endl;
	//-------------------------------------
		/*ofstream file("Energies_Delta=" + to_string((int)Delta) + ".txt");
		for (int k = 0; k < eigenvalues.size(); k++) {
			file << pow(2, L) + 1 << "\t\t" << eigenvalues(k) << endl;
		}
		file.close();*/

}
//----------------------------------------------------

//Prints eigenvalues to file
void Heisenberg1D::print_eigenvalues(ofstream& fileE){
	for (int k = 0; k < eigenvalues.size(); k++) {
		fileE << L << "\t" << eigenvalues(k) << endl;
	}
}
//--------------------------------

//Checks if eigenvectors and eigenvalues are correct
Col<int> Heisenberg1D::check_eigenvalues() {
	Col<int> result(N);
	mat energy_mat(N, N, fill::zeros);
	energy_mat = eigenvectors.i() * H * eigenvectors;
	for (int k = 0; k < N; k++) {
		if ( (eigenvalues(k) - energy_mat(k, k)) < 1e-12) result(k) = 1;
		else result(k) = 0;
	}
	energy_mat.~Mat();
	return result;
 }
//----------------------------------------------------

//Calculates the energy gap in the system
double Heisenberg1D::get_energy_gap() {
	sort(eigenvalues,"ascend"); //sorting ascending
	int k = 1;
	while ( fabs(eigenvalues(0) - eigenvalues(k)) < 1e-12) {
		k++; //lowest level could be degenerated
		if (k >= eigenvalues.size()) return 0;
	}
	return (eigenvalues(k) - eigenvalues(0)) / (L + 0.0);
}
//----------------------------------------------------


//Calculates the Heat Capacity of the system
void Heisenberg1D::Heat_Capacity(double T_step, double T_0, double T_end, ofstream& savefile) {
		double T = T_0; //temperature
		double energy_av; //average of energy E
		double energy2_av; //average of E^2
		while (T <= T_end) {
			//ED Heat capacity
				double Partition_Function = 0; 
				energy_av = 0; energy2_av = 0;
				for (int j = 0; j < N; j++) {
					Partition_Function += std::exp(-eigenvalues(j) / T); //partition function(T)
					energy_av += eigenvalues(j) * std::exp(-eigenvalues(j) / T); //average energy(T)
					energy2_av += eigenvalues(j) * eigenvalues(j) * std::exp(-eigenvalues(j) / T); //average energy^2(T)
				}
				energy_av = energy_av / Partition_Function;
				energy2_av = energy2_av / Partition_Function;
				double heat_capacity = (energy2_av - energy_av * energy_av) / T / T / (L + 0.0);
			savefile << T << "\t\t" << heat_capacity << endl; //save heat capacity to file
			T += T_step;
		}
}
//----------------------------------------------------


// CALCULATES THE AVERAGE ENERGY OF THE SYSTEM and evaluates the temperature of the domain-wall state
double Heisenberg1D::get_domain_wall_state_temp(double T_step, double T_end, string del) {
	ofstream file_energy("Average energy_Delta=" + del + ".txt");
	
	vec domain_wall_state(N, fill::zeros);
	int idx = mapping(binary_to_int(first_vec));
	domain_wall_state(idx) = 1;
	double energy_neel = dot(domain_wall_state, H * domain_wall_state) / L;

	std::cout << energy_neel << endl;
	domain_wall_state.~vec();

	double T = -T_end;
	int k = 0;
	double domain_wall_temperature = DBL_MAX; // T = infinity if average energy does not reach domain_wall_state energy
	while(T < 0){
		double Partition_Function = 0, av_energy = 0;
		for (int j = 0; j < N; j++) {
			Partition_Function += std::exp(-eigenvalues(j) / T); //partition function(T)
			av_energy += eigenvalues(j) * std::exp(-eigenvalues(j) / T); //average energy(T)
		}
		av_energy = av_energy / Partition_Function / (L + 0.0);
		if (fabs(energy_neel - av_energy) < T_step) domain_wall_temperature = T;
		file_energy << T << "\t" << av_energy << "\t" << energy_neel << endl;
		T += T_step;
	}
	//if( (energy_neel - av_energy(av_energy.size() - 1) ) < 1e-2 ) domain_wall_temperature = DBL_MAX;
	file_energy.close();
	return domain_wall_temperature;
}
//---------------------------------------------

//Generates coefficients vetor, which keeps the dot product of domain-wall state and each eigenvector
void Heisenberg1D::generate_coefficients_vector(cx_vec wavefunction) {
	coeff = cx_vec(N,fill::zeros);
	vec zero(N, fill::zeros);
	for (int j = 0; j < N; j++) {
		vec base = eigenvectors.col(j);
		cx_vec base_vec(base, zero);
		coeff(j) = cdot(base_vec, wavefunction);

		base.clear(); base_vec.clear();
	}

	zero.~vec();
}
//---------------------------------------------

//Time evolution of input wavefunction (vector) for time-independent hamitlonian
cx_vec Heisenberg1D::Time_evolution_stationary(cx_vec wavefunction, double time) {
	//cx_vec domain_wall_vec = wavefunction; //domain_wall_state |00...011..1>
	wavefunction.zeros();
	for (int k = 0; k < N; k++) {
		wavefunction += exp(-1i * time * eigenvalues(k)) * coeff(k) * eigenvectors.col(k);// cdot(base_vec, domain_wall_vec);
	}
	//domain_wall_vec.~cx_vec();

	return wavefunction;
}
//---------------------------------------------

//Time evolution of input wavefunction (vector) for time-dependent hamitlonian
cx_vec Heisenberg1D::Time_evolution(cx_vec wavefunction, double dt) {
	cx_vec fun1(wavefunction), fun2(wavefunction), fun3(wavefunction); // fun4(wavefunction);
	fun1 = -1i * dt * H * wavefunction;
	fun2 = -1i * dt * H * fun1 / 2;
	fun3 = -1i * dt * H * fun2 / 3;
	wavefunction = wavefunction + fun1 + fun2 + fun3;
	fun1.~cx_vec(); fun2.~cx_vec(); fun3.~cx_vec();
	return wavefunction;
}
//---------------------------------------------

// Prints to file evolution of domain_wall_state wavefunction
void Heisenberg1D::print_evolution(string del) {

	cx_vec wavefunction(N, fill::zeros);
	int idx = mapping(binary_to_int(first_vec));
	wavefunction(idx) = 1; //first base_vector is 0011 = first_vec

	//create a locale that has the ctype category  copied from the "en_US.UTF-8"
		locale loc = locale(locale(), "en_US.UTF8", locale::ctype);
		wofstream file("Evolution visualization_Delta=" + del + ".txt");
		file.imbue(loc); // adding the locale to the stream
	//-------------------------------------------------------

	double t = 0, dt = 0.01;
	vector<int> base_vec(L);
	int p = 0;
	while (t <= 2 * L) {
		if ( static_cast<int>( t/dt) == 50*p) {
			file << t << "\t";
			for (int j = 0; j < N; j++) {
				if (wavefunction(j) != 0.0) {
					base_vec = int_to_binary(mapping_inv(j), L);
					file << fixed << setprecision(2) << wavefunction(j) << "*|";
					for (int k = 0; k < L; k++) {
						if (base_vec[k] == 0) { file << L'\u2193'; }
						if (base_vec[k] == 1) { file << L'\u2191'; }
					}
					file << ">";
					if (j < N - 1) file << "  +  ";
				}
			}
			file << "\n"; 
			p++;
		}
		wavefunction = Time_evolution(wavefunction, dt);
		t += dt;
	}
	file.close();
	loc.~locale();
}
//---------------------------------------------

//Calculates the magnetization on each site over time evolution -> output to file
void Heisenberg1D::Magnetization_evolution(string del) {

	double dt = 0.1;
	double time = 0;

	ofstream file_evolution("Unitary evolution_Delta=" + del + ".txt");
	ofstream file_evolution2("Neel state evolution_Delta=" + del + ".txt");

	//normal
	/*cx_vec wavefunction(N, fill::zeros);
	int idx = mapping(binary_to_int(first_vec));
	wavefunction(idx) = 1; //first base_vector is 0011 = first_vec
	cx_vec domain_wall_vec = wavefunction;
	generate_coefficients_vector(domain_wall_vec);
	for (time = 0; time <= L; time += dt) { // lentgh is 2*L/dt = 200 * L
		//wavefunction = Time_evolution_stationary(domain_wall_vec, time); //initial vector stays the same
		wavefunction.zeros();
		for (int k = 0; k < N; k++) {
			wavefunction += exp(-1i * time * eigenvalues(k)) * coeff(k) * eigenvectors.col(k);// cdot(base_vec, domain_wall_vec);
		}
		if (static_cast<int>(time / dt) % 5 == 0) file_evolution2 << " 'hed'\n " << endl; //header only for separation in file, easier to plot with gnuplot

		for (int k = 0; k < L; k++) {
			double Sz = 0;
			for (int j = 0; j < N; j++) {
				vector<int> vect = int_to_binary(mapping_inv(j), L); //jth base_vector
				Sz = Sz + (vect[k] - 0.5) * abs(wavefunction(j)) * abs(wavefunction(j)); //otput should be real, but C++ forbids such typecast
				vect.~vector();
			}
			file_evolution << k + 1 << "\t" << time << "\t" << Sz << endl;
			if (static_cast<int>(time / dt) % 5 == 0) { //write to file every 0.1 time evolution
				file_evolution2 << k + 1 << "\t" << time << "\t" << Sz << endl;
			}
		}
		//wavefunction = Time_evolution(wavefunction, dt);
	}
	domain_wall_vec.~cx_vec();
	wavefunction.~cx_vec();*/

	//Parallel with openMP
	/*vector<double> time_vec;
	while (time <= L) {
		time_vec.push_back(time);
		time += dt;
	}
	int nloop = time_vec.size();
	results = vector<vector<string>>(nloop);
	cx_vec wavefunction(N, fill::zeros);
	int idx = mapping(binary_to_int(first_vec));
	wavefunction(idx) = 1; //first base_vector is 0011 = first_vec
	cx_vec domain_wall_vec = wavefunction;
	generate_coefficients_vector(domain_wall_vec);
	for (int m = 0; m < time_vec.size(); m++) {
		//wavefunction = Time_evolution_stationary(domain_wall_vec, time); //initial vector stays the same
			wavefunction.zeros();
#pragma loop(hint_parallel(16))
			for (int k = 0; k < N; k++) {
				wavefunction += exp(-1i * time_vec[m] * eigenvalues(k)) * coeff(k) * eigenvectors.col(k);// cdot(base_vec, domain_wall_vec);
			}

			results[m] = vector<string>(L);

			for (int k = 0; k < L; k++) {
				double Sz = 0;
				for (int j = 0; j < N; j++) {
					vector<int> vect = int_to_binary(mapping_inv(j), L); //jth base_vector
					Sz = Sz + (vect[k] - 0.5) * abs(wavefunction(j)) * abs(wavefunction(j)); //otput should be real, but C++ forbids such typecast
					vect.~vector();
				}
				results[m][k] = to_string(k + 1) + "\t" + to_string(time_vec[m]) + "\t" + to_string(Sz);
				//file_evolution << k + 1 << "\t" << time_vec[m] << "\t" << Sz << endl;
			}
			//wavefunction = Time_evolution(wavefunction, dt);
			domain_wall_vec.~cx_vec();
			wavefunction.~cx_vec();
	}

	for (int j = 0; j < nloop; j++) {
		for (int k = 0; k < L; k++) {
			file_evolution << results[j][k] << endl;
		}
	}
	results.~vector();
	time_vec.~vector();*/

	//Parallel using <thread> and <mutex>
	const size_t nthreads = thread::hardware_concurrency();
	std::cout << "Using " << nthreads << " threads" << endl;
	vector<thread> threads(nthreads);
	vector<double> time_vec;
	while (time <= 2*L) {
		time_vec.push_back(time);
		time += dt;
	}
	int nloop = time_vec.size();
	results = vector<vector<string>>(nloop);
	int thread_begin = 0, thread_end = -1;
	for (int t = 0; t < nthreads; t++) {
		thread_begin = t * nloop / nthreads;
		thread_end = ( (t + 1) == nthreads ? nloop : (t + 1) * nloop / nthreads );
		threads[t] = thread(&Heisenberg1D::Magnetization_1thread, this, thread_begin, thread_end, time_vec, ref(file_evolution));
	}
	for (auto& t : threads) t.join();
	for (auto& t : threads) t.~thread();

	for (int j = 0; j < nloop; j++) {
		for (int k = 0; k < L; k++) {
			file_evolution << results[j][k] << endl;
		}
	}
	time_vec.~vector();
	results.~vector();

	file_evolution.close();
	file_evolution2.close();
}

void Heisenberg1D::Magnetization_1thread(const int beg, const int end, vector<double> time_vec, ofstream& file_evolution) {
	
	cx_vec wavefunction(N, fill::zeros);
	my_mutex.lock();
	int idx = mapping(binary_to_int(first_vec));
	my_mutex.unlock();
	wavefunction(idx) = 1; //first base_vector is 0011 = first_vec
	cx_vec domain_wall_vec = wavefunction;
	my_mutex.lock();
	generate_coefficients_vector(domain_wall_vec);
	my_mutex.unlock();

	for (int m = beg; m < end; m++) {
		my_mutex.lock();
			results[m] = vector<string>(L); //allocate vector of length L on m-the vector position -> matrix L x timez_vec.size()
			//mat eigenVec = eigenvectors; //creating matrix for each thread - RAM will explode
			//cx_vec coefficients = coeff;
		my_mutex.unlock();
		wavefunction.zeros();
		//Evolution
		for (int k = 0; k < N; k++) {//here problems with shared coeff and eigenvector, defininf mutex in loop slowers down the program execution
			wavefunction += exp(-1i * time_vec[m] * eigenvalues(k)) * coeff(k) * eigenvectors.col(k);// cdot(base_vec, domain_wall_vec);
		}
		//eigenVec.~Mat(); coefficients.~cx_vec();
		//if (static_cast<int>(time / dt) % 20 == 0) file_evolution2 << " 'hed'\n " << endl; //header only for separation in file, easier to plot with gnuplot
		for (int k = 0; k < L; k++) {
			double Sz = 0;
			for (int j = 0; j < N; j++) {
				vector<int> vect = int_to_binary(mapping_inv(j), L); //jth base_vector
				Sz = Sz + (vect[k] - 0.5) * abs(wavefunction(j)) * abs(wavefunction(j)); //output should be real, but C++ forbids such typecast
				vect.~vector();
			}
			my_mutex.lock();
				results[m][k] = to_string(k + 1) + "\t" + to_string(time_vec[m]) + "\t" + to_string(Sz);
			//file_evolution << k + 1 << "\t" << time_vec[m] << "\t" << Sz << endl;
			my_mutex.unlock();
			/*if (static_cast<int>(time / dt) % 20 == 0) { //write to file every 0.1 time evolution
				file_evolution2 << k + 1 << "\t" << time_vec[m] << "\t" << Sz << endl;
			}*/
		}
	}
	domain_wall_vec.~cx_vec();
	wavefunction.~cx_vec();
}
//-----------------------------------------------------------------------------



//Overlap domain-wall state (t) and state definde by domain-wall temperature: phi_T=\sum_n exp(-\beta e_n) |n>
void Heisenberg1D::Overlap(double time_end, double domain_wall_temperature, string del) {
	
	ofstream file("Overlap time evolution_Delta=" + del + ".txt");
	double t = 0.00, dt = 0.1;
	complex<double> overlap;

	cx_vec wavefunction(N, fill::zeros);
	int idx = mapping(binary_to_int(first_vec)); //position in wavefunction for domain-wall basevector
	wavefunction(idx) = 1; //domain-wall state
	cx_vec domain_wall_vec = wavefunction;

	// Declaration of phi_T, no need because analitycal formula
	/*double partition_function = 0;
	this->phi_T = cx_vec(N, fill::zeros);
	vec temp(N, fill::zeros), zero(N,fill::zeros);
	for (int j = 0; j < N; j++) {
		partition_function += exp(-eigenvalues(j) / domain_wall_temperature);
		temp += exp(-eigenvalues(j) / domain_wall_temperature) * eigenvectors.col(j);
	}
	phi_T = cx_vec(temp, zero);
	temp.~vec(); zero.~vec();
	phi_T = phi_T / partition_function;
	cx_vec varphi_T(phi_T);*/
	double partition_function = 0;
	for (int j = 0; j < N; j++) partition_function += exp(-eigenvalues(j) / domain_wall_temperature);

	generate_coefficients_vector(wavefunction); //  --> cx_vec coeff
	while (t <= time_end) {
		overlap = 0;
		//cdot(wavefunction, base_vec_complex) = coeff(j)* = conj(coeff(j))
		for (int j = 0; j < N; j++) 
			overlap += exp((-1 / domain_wall_temperature + 1i * t) * eigenvalues(j)) * conj(coeff(j)) / partition_function;//no need for time evoultion of wavefunction

		file << t << "\t\t" << real(overlap) << "\t\t" << imag(overlap) << endl;
		t += dt;
	}
	wavefunction.~cx_vec();
	file.close();
	std::cout << "Overlap finished" << endl;
}
//----------------------------------------------------

// Calculates the spin correlation function at T=0
void Heisenberg1D::Spin_Correlations(string del) {

	cx_vec wavefunction(N, fill::zeros);
	int idx = mapping(binary_to_int(first_vec));
	wavefunction(idx) = 1; //domain-wall state
	cx_vec domain_wall_state(wavefunction);
	generate_coefficients_vector(domain_wall_state);

	double dt = 2;
	double t = 0;
	ofstream file("Spin correlator_Delta="+del+".txt");
	int a = 0; // (L/2)-th site is on position a, because vector starts from 0
	vector<vector<double>> Correlator(L); //correlation function <S_i * S_i+1> - <S_i><S_i+1> (z component)

	int p = 0;
	while (t <= 4*L) {
		wavefunction = Time_evolution_stationary(domain_wall_state, t);

		file << "'hed'\n" << endl;
		for (int r = 0; r < L-1; r++) { //corelation through whole chain
			Correlator[r].resize(p + 1);
			double S12 = 0, S1 = 0, S2 = 0;

			for (int j = 0; j < N; j++) {
				vector<int> vect = int_to_binary(mapping_inv(j), L);
				S12 += (vect[r] - 0.5) * (vect[r+1] - 0.5) * abs(wavefunction(j)) * abs(wavefunction(j));
				S1 += (vect[r] - 0.5) * abs(wavefunction(j)) * abs(wavefunction(j));
				S2 += (vect[r+1] - 0.5) * abs(wavefunction(j)) * abs(wavefunction(j));
			}
			Correlator[r][p] = S12 - S1 * S2;
			file << r + 1 << "\t" << Correlator[r][p] << "\t" << t << endl;
		}
		file << "\n" << endl; 
		t += dt;
		p++;
	}
	Correlator.~vector();
	file.close();
}
//-----------------------------------------

// Calculates the quantum fidelty: matrix element of unitary evolution operator
void Heisenberg1D::Quantum_fidelity(cx_vec wavefunction) {
	double time = 0, dt = 0.1;
	cx_double Loschmidt_echo;
	/*vector<int> neel_vec(L);
	for (int k = 0; k < L; k++) {
		if (k % 2 == 0) neel_vec[k] = 1;
		else neel_vec[k] = 0;
	}
	cx_vec wavefunction(N, fill::zeros);
	int idx; ofstream echo_file;
	if (whichInitVector == 0) {
		idx = mapping(binary_to_int(neel_vec));
		echo_file.open("Loschmidt_echo_neel-state_delta=" + del + ".txt");
	}
	else {
		idx = mapping(binary_to_int(first_vec));
		echo_file.open("Loschmidt_echo_domain-wall-state_delta=" + del + ".txt");
	}

	wavefunction(idx) = 1; //first base_vector is 0011 = first_vec*/
	ofstream echo_file("Loschmidt_echo.txt");
	cx_vec wavefunction_0 = wavefunction;
	generate_coefficients_vector(wavefunction); //  --> coeff(n) = <n|psi(0)>
	
	while (time <= 5 * L) {
		//wavefunction = Time_evolution_stationary(wavefunction_0, time);
		//Loschmidt_echo = abs(cdot(wavefunction_0, wavefunction)); // |<psi(0), psi(t)>|
		Loschmidt_echo = 0;
		for (int k = 0; k < N; k++) Loschmidt_echo += exp(-1i * time * eigenvalues(k)) * abs(coeff(k)) * abs(coeff(k));
		echo_file << time << "\t\t" << abs(Loschmidt_echo) * abs(Loschmidt_echo) << endl; // |<psi(0), psi(t)>|^2
		time += dt;
	} 
	//neel_vec.~vector();
	echo_file.close();
	std::cout << "Loschmidt echo finished" << endl;
}
//-----------------------------------------

//----------------------------------------------------------------------------------------------
//----------------------------------------Linear responce------------------------------
//----------------------------------------------------------------------------------------------

void Heisenberg1D::spinStructureFactor_lin_responce(double temperature) {
	double partition_function = 0;
#pragma omp parallel for shared(temperature) reduction(+: partition_function)
	for (int j = 0; j < N; j++) partition_function += exp(-eigenvalues(j) / temperature);
	double q, omega = 0; // eigenvalues(0) - eigenvalues(N - 1);
	double domega = 0.01;
	
	vec zero(N, fill::zeros);
	vector<double> omega_vec;
	while (omega <= eigenvalues(N - 1) - eigenvalues(0)) {
		omega_vec.push_back(omega);
		omega += domega;
	}
	int nloop = omega_vec.size();

	stringstream s1, s2;
	s1 << fixed << setprecision(3) << temperature;
	s2 << fixed << setprecision(0) << Delta;
	string del = s2.str(); string T = s1.str(); 
	ofstream SpinFactorFile;
	if( temperature == DBL_MAX)  SpinFactorFile.open("SSF_T=inf_Delta=" + del + ".txt");
	else if (temperature == -DBL_MAX)  SpinFactorFile.open("SSF_T=-inf_Delta=" + del + ".txt");
	else SpinFactorFile.open("SSF_T="+T+"_Delta="+del+".txt");
	//ofstream SpinFactorFile("SSF_Delta=" + del + ".txt");
	//ofstream file2("SSF-snap_T=" + T + "_Delta=" + del + ".txt");

	vector<vector<string>> resultSF(L + 1);
	vector<vector<double>> SF(L + 1);
	for (int ID = 0; ID <= L; ID++) { //sum over q
		resultSF[ID] = vector<string>(nloop);
		SF[ID] = vector<double>(nloop, 0.0);
		q = 2 * pi / L * ID;
		cx_mat Sq(mat(N, N, fill::zeros), mat(N, N, fill::zeros));
		mat Sz(N, N, fill::zeros); //z-spin matrix for k^th site : S_k^z
		for (int k = 0; k < L; k++) {//Calculate  S_q^z
#pragma omp parallel for shared(Sz) schedule(static) // approx 3x faster & correct results
			for (int p = 0; p < N; p++) {
				vector<int> vect = int_to_binary(mapping_inv(p), L); //j^th base_vector
				Sz(p, p) = (double)vect[k] - 0.5;
			}
			Sq += exp(1i * q * (k + 0.0)) * Sz;
		}
		Sz.~Mat();
		double SpinFactor = 0;
		// planting here omp parallel produces data race - check, because speed boost = 8x, while below = 3x 
		for(int w=0; w < nloop; w++){
			omega = omega_vec[w];
			SpinFactor = 0;	
#pragma omp parallel for shared(omega, q, partition_function, temperature, resultSF, zero, Sq)\
reduction(+: SpinFactor) schedule(static) // approx 3x faster & correct results
			for (int n = 0; n < N; n++) {
				for (int m = n; m < N; m++) { // m = n, because only for E_m > E_n we have constribution - eigenvalues are sorted and omega>0
					double omega_mn = eigenvalues(m) - eigenvalues(n);
					if ( (omega_mn < omega + domega/2.) && (omega_mn >= omega - domega/2.) ) {
						//double matrix_element = abs( cdot(cx_vec(eigenvectors.col(n), zero), Sq * eigenvectors.col(m)) );
						//matrix_element *= matrix_element;
						cx_double matrix_element = 0;
						for (int o = 0; o < N; o++) matrix_element += eigenvectors.col(n)(o) * Sq(o, o) * eigenvectors.col(m)(o);
						SpinFactor += exp(-eigenvalues(n) / temperature) * abs(matrix_element) * abs(matrix_element) / partition_function;
					}
				}
			}
			resultSF[ID][w] = to_string(q) + "\t\t" + to_string(omega_vec[w]) + "\t\t" + to_string(SpinFactor);
		}
		std::cout << q << endl;
		Sq.~cx_mat();
	}
	for (int j = 0; j < L + 1; j++) {
		for (int k = 0; k < nloop; k++)
			SpinFactorFile << resultSF[j][k] << endl;
			//SpinFactorFile << q << "\t\t" << omega_vec[k] << "\t\t" << SF[j][k] << endl;
	}
	resultSF.~vector();
	SF.~vector();
	zero.~vec();
	omega_vec.~vector();
	SpinFactorFile.close(); //file2.close();
}

void Heisenberg1D::spinStructureFactor_T_0() {
	double q, omega = 0; // eigenvalues(0) - eigenvalues(N - 1);
	double domega = 0.05;

	vec zero(N, fill::zeros);
	vector<double> omega_vec;
	while (omega <= pi) {
		omega_vec.push_back(omega);
		omega += domega;
	}
	int nloop = omega_vec.size();

	stringstream s2;
	s2 << fixed << setprecision(0) << Delta;
	string del = s2.str();
	ofstream SpinFactorFile;
	SpinFactorFile.open("SSF_T=0_Delta=" + del + ".txt");
	//ofstream SpinFactorFile("SSF_Delta=" + del + ".txt");
	//ofstream file2("SSF-snap_T=" + T + "_Delta=" + del + ".txt");

	vector<vector<string>> resultSF(L + 1);
	for (int ID = 0; ID <= L; ID++) { //sum over q
		resultSF[ID] = vector<string>(nloop);
		q = 2 * pi / L * ID;
		cx_mat Sq(mat(N, N, fill::zeros), mat(N, N, fill::zeros));
		mat Sz(N, N, fill::zeros); //z-spin matrix for k^th site : S_k^z
		for (int k = 0; k < L; k++) {//Calculate  S_q^z
#pragma omp parallel for shared(Sz) schedule(static) // approx 3x faster & correct results
			for (int p = 0; p < N; p++) {
				vector<int> vect = int_to_binary(mapping_inv(p), L); //j^th base_vector
				Sz(p, p) = (double)vect[k] - 0.5;
			}
			//psi_m += exp(1i * q * (k + 0.0)) * Sz * eigenVec.col(m);
			Sq += exp(1i * q * (k + 0.0)) * Sz;
		}
		Sz.~Mat();
		double SpinFactor = 0;
		// planting here omp parallel produces data race - check, because speed boost = 8x, while below = 3x 
		for (int w = 0; w < nloop; w++) {
			omega = omega_vec[w];
			SpinFactor = 0;
#pragma omp parallel for shared(omega, q, resultSF, Sq) reduction(+: SpinFactor) schedule(static) // approx 3x faster & correct results
			for (int m = 0; m < N; m++) { // m = n, because only for E_m > E_n we have constribution - eigenvalues are sorted and omega>0
				double omega_mn = eigenvalues(m) - eigenvalues(0);
				if ((omega_mn < omega + domega / 2.) && (omega_mn >= omega - domega / 2.)) {
					cx_double matrix_element = 0;
					for (int o = 0; o < N; o++) matrix_element += eigenvectors.col(0)(o) * Sq(o, o) * eigenvectors.col(m)(o);
					SpinFactor += abs(matrix_element) * abs(matrix_element);
				}
			}
			resultSF[ID][w] = to_string(q) + "\t\t" + to_string(omega_vec[w]) + "\t\t" + to_string(SpinFactor);
		}
		std::cout << q << endl;
		Sq.~cx_mat();
	}

	for (int j = 0; j < L + 1; j++) {
		for (int k = 0; k < nloop; k++) {
			SpinFactorFile << resultSF[j][k] << endl;
			//file2 << resultSF[j][k] << endl;
		}
		//file2 << "\n\n" << endl;
	}
	resultSF.~vector();
	zero.~vec();
	omega_vec.~vector();
	SpinFactorFile.close(); //file2.close();
}

double Heisenberg1D::get_average_gap() {
	double h = 0;
	for (int i = 0; i < N - 1; i++) h += eigenvalues(i + 1) - eigenvalues(i);
	return h / (N - 1.0);
}

void Heisenberg1D::SSF_magnetic_field(double h) {

	double q, omega = 0, time = 0;
	double domega = 0.1;
	double dt = domega;
	vector<double> omega_vec;
	while (omega <= pi) {
		omega_vec.push_back(omega);
		omega += domega;
	}
	int nloop = omega_vec.size();
	vector<double> time_vec;
	while (time <= 5*L) {
		time_vec.push_back(time);
		time += dt;
	}
	int no_iterations = time_vec.size();

	stringstream s1, s2;
	s1 << fixed << setprecision(4) << 0 << "_h=" << h;
	s2 << fixed << setprecision(0) << Delta;
	string del = s2.str(); string T = s1.str();
	//ofstream SpinFactorFile("SSF_T=" + T + "_Delta=" + del + ".txt");
	ofstream SpinFactorFile("L=10test.txt");

	vector<vector<string>> resultSF(L + 1);
	for (int ID = 0; ID <= L; ID++) { //sum over q
		q = 2 * pi / L * ID;
		Magnetic_perturbance(q, h); // create perturbation
		// New temporary eigenval;ues, eigenvectors
			vec eigenVal(N);
			mat eigenVec(N, N);
		eig_sym(eigenVal, eigenVec, H + H_magnetic); // diagonalization with perturbation
		cx_vec wavefunction = cx_vec(eigenVec.col(0), vec(N, fill::zeros)); //initial state - ground state
			eigenVal.~vec(); eigenVec.~Mat(); 
			wavefunction.zeros();
			int idx = mapping(binary_to_int(first_vec));
			wavefunction(idx) = 1;
		generate_coefficients_vector(wavefunction); //generates the coefficients:  <wavefunction|n>, where |n> are eigenevectors of the unperturbed hamiltonian
		vector<cx_double> S_qt(no_iterations);
		for (int m = 0; m < no_iterations; m++) { //loop over time
			//Evolution
				wavefunction.zeros();
				for (int k = 0; k < N; k++) {//here problems with shared coeff and eigenvector, defining mutex in loop slowers down the program execution
					wavefunction += exp(-1i * time_vec[m] * eigenvalues(k)) * coeff(k) * eigenvectors.col(k);
				}
		//--------------FFT(q){<Sz(t)>}------------------
			double re = 0, imag = 0;
//#pragma omp parallel for shared(wavefunction, q) reduction(+:real, imag) 
				/*for (int j = 0; j < N; j++) {
					vector<int> vect = int_to_binary(mapping_inv(j), L); //jth base_vector
					//Sz = Sz + (vect[k] - 0.5) * abs(wavefunction(j)) * abs(wavefunction(j));
					for (int k = 0; k < L; k++) {
						re += (vect[k] - 0.5) * abs(wavefunction(j)) * abs(wavefunction(j)) * cos(q * static_cast<double>(k));
						imag += (vect[k] - 0.5) * abs(wavefunction(j)) * abs(wavefunction(j)) * sin(q * static_cast<double>(k));
					}
					vect.~vector();
					//Sq += exp(1i * q * static_cast<double>(k)) * Sz;
				}*/
				cx_double Sq = 0;
				for (int k = 0; k < L; k++) {
					double Sz = 0;
					for (int j = 0; j < N; j++) {
						vector<int> vect = int_to_binary(mapping_inv(j), L); //jth base_vector
						Sz = Sz + (vect[k] - 0.5) * abs(wavefunction(j)) * abs(wavefunction(j)); //output should be real, but C++ forbids such typecast
						vect.~vector();
					}
					SpinFactorFile << k + 1 << "\t\t" << time_vec[m] << "\t\t" << Sz << endl;
					Sq += exp(1i * q * static_cast<double>(k)) * Sz;
				}
				S_qt[m] = re + 1i * imag; //vector of S(q,t) for different t
		}
		//------------------FFT(omega){ S(q,t) }--------------
		//resultSF[ID] = vector<string>(nloop);
			for (int w = 0; w < nloop; w++) {
				 /*double S_qw_re = 0, S_qw_Im = 0;
#pragma omp parallel for reduction(+: S_qw_re, S_qw_Im)
					S_qw_re += real(exp(1i * omega_vec[w] * time_vec[t]) * S_qt[t] * dt);
					S_qw_Im += imag(exp(1i * omega_vec[w] * time_vec[t]) * S_qt[t] * dt);
					S_qw = S_qw_re + 1i * S_qw_Im;*/
				cx_double S_qw = 0;
				for (int t = 0; t < no_iterations; t++)
					S_qw += exp(1i * omega_vec[w] * time_vec[t]) * S_qt[t] * dt;
				//SpinFactorFile << q << "\t\t" << omega_vec[w] << "\t\t" << abs(S_qw) << "\t\t" << real(S_qw) << "\t\t" << imag(S_qw) << endl;
				//resultSF[ID][w] = to_string(q) + "\t\t" + to_string(omega_vec[w]) + "\t\t" + to_string(real(S_qw));
			}
		std::cout << "q = " << q << endl;
	}
	/*for (int j = 0; j < L / 2 + 1; j++) {
		for (int k = 0; k < nloop; k++)
			SpinFactorFile << resultSF[j][k] << endl;
		//SpinFactorFile << q << "\t\t" << omega_vec[k] << "\t\t" << SF[j][k] << endl;
	}*/
	resultSF.~vector();
	omega_vec.~vector(); time_vec.~vector();
	SpinFactorFile.close(); //file2.close();
}



//----------------------------------------------------------------------------------------------
//----------------------------------------LANCZOS ALGORITHM AND MORE------------------------------
//----------------------------------------------------------------------------------------------

// Using Lanczos power iteration method - diagonalization - already with matrix build starting from random vector
void Heisenberg1D::Lanczos_Diagonalization(int m) {
	this->Krylov_space = mat(N, m);
	this->H_Lanczos = mat(m, m, fill::zeros);

	//	Generate random first vector
	srand(time(NULL));
	for (int j = 0; j < N; j++) {
		Krylov_space.col(0)(j) = static_cast<double>(rand()) / (RAND_MAX + 0.0) - 0.5;
		//fi_j_1(j) = static_cast<double>(rand()) / (RAND_MAX + 0.0) - 0.5;
	}
	double beta = dot(Krylov_space.col(0), Krylov_space.col(0));
	Krylov_space.col(0) = Krylov_space.col(0) / sqrt(beta); //normalized fi_0
	vec tmp = H * Krylov_space.col(0);

	double alfa = dot(Krylov_space.col(0), tmp);
	tmp = tmp - alfa * Krylov_space.col(0);
	H_Lanczos(0, 0) = alfa;
	for (int j = 1; j < m; j++) {
		beta = sqrt(dot(tmp, tmp));
		Krylov_space.col(j) = tmp / beta;

		tmp = H * Krylov_space.col(j);
		alfa = dot(Krylov_space.col(j), tmp);
		tmp = tmp - alfa * Krylov_space.col(j) - beta * Krylov_space.col(j - 1);

		H_Lanczos(j, j) = alfa;
		H_Lanczos(j, j - 1) = beta;
		H_Lanczos(j - 1, j) = beta;

	}

	this->Lanczos_eigenVal = vec(m);
	mat eigenVec(m, m);
	eig_sym(Lanczos_eigenVal, eigenVec, H_Lanczos); // no eigenvectors for now
	//cout << Lanczos_eigenVal.t();
	this->Lanczos_eigenVec = mat(N, m); // m eigenvectors of lanczos iteration
	Lanczos_eigenVec = Krylov_space * eigenVec;
	eigenVec.~Mat();

	//Lanczos_gs.push_back(Lanczos_eigenVal(0));
	tmp.~vec();
}
//Print eigenvalues to stream
void Heisenberg1D::Lanczos_convergence(int m, stringstream& ss) {
	for (int k = 0; k < Lanczos_eigenVal.size(); k++) {
		ss << m << "\t" << Lanczos_eigenVal(k) << "\n";
	}
}
void Heisenberg1D::Lanczos_gs_energy() {
	int idx = 0;
	ofstream file("Relative_lanczos_Delta=" + to_string((int)Delta) + ".txt");
	for (int k = 0; k < Lanczos_gs.size(); k++) {
		file << k + 1 << "\t\t" << (Lanczos_gs[k] - eigenvalues(0)) / eigenvalues(0) << endl;
		if ((fabs(Lanczos_gs[k] - eigenvalues(0)) <= 1e-18) && idx < 1) {
			cout << "L = " << L << "\tM = " << k << endl;
			idx++;
		}
	}
	file.close();
}
//-------------------------------------------------------

//Lanczos Spin Structure Factor
void Heisenberg1D::Build_Lanczos_Hamil(cx_vec initial_vec, int m) {
	this->Krylov_space = mat(N, m);
	this->H_Lanczos = mat(m, m, fill::zeros);

	Krylov_space.col(0) = real(initial_vec);

	double beta = dot(Krylov_space.col(0), Krylov_space.col(0));
	Krylov_space.col(0) = Krylov_space.col(0) / sqrt(beta); //normalized fi_0
	vec tmp = H * Krylov_space.col(0);

	double alfa = dot(Krylov_space.col(0), tmp);
	tmp = tmp - alfa * Krylov_space.col(0);
	H_Lanczos(0, 0) = alfa;
	for (int j = 1; j < m; j++) {
		beta = sqrt(dot(tmp, tmp));
		Krylov_space.col(j) = tmp / beta;

		tmp = H * Krylov_space.col(j);
		alfa = dot(Krylov_space.col(j), tmp);
		tmp = tmp - alfa * Krylov_space.col(j) - beta * Krylov_space.col(j - 1);

		H_Lanczos(j, j) = alfa;
		H_Lanczos(j, j - 1) = beta;
		H_Lanczos(j - 1, j) = beta;

	}
	tmp.~vec();
}
void Heisenberg1D::SSF_Lanczos() {
	double q, omega = 0; // eigenvalues(0) - eigenvalues(N - 1);
	double domega = 0.01;
	int M = Lanczos_eigenVal.size();

	stringstream s1;
	s1 << fixed << setprecision(0) << Delta;
	string del = s1.str();
	ofstream SpinFactorFile;
	SpinFactorFile.open("SSF_Lanczos_M=" + to_string(Lanczos_eigenVal.size()) + "_Delta=" + del + ".txt");

	for (int ID = 0; ID <= L; ID++) { //sum over q
		q = 2 * pi / L * ID;
		cx_mat Sq(mat(N, N, fill::zeros), mat(N, N, fill::zeros));
		mat Sz(N, N, fill::zeros); //z-spin matrix for k^th site : S_k^z
		for (int k = 0; k < L; k++) {//Calculate  S_q^z
#pragma omp parallel for shared(Sz) schedule(static)
			for (int p = 0; p < N; p++) {
				vector<int> vect = int_to_binary(mapping_inv(p), L); //j^th base_vector
				Sz(p, p) = (double)vect[k] - 0.5;
			}
			Sq += exp(1i * q * (k + 0.0)) * Sz;
		}
		Sz.~Mat();

		//Second Lanczos procedure
			vec ground_state = Lanczos_eigenVec.col(0);

			cx_vec SqSq_GS(N, fill::zeros);		// = Sq.t() * Sq * ground_state;
#pragma omp parallel for shared(SqSq_GS, Sq, ground_state) schedule(static)
			for (int s = 0; s < N; s++) SqSq_GS(s) = conj(Sq(s, s)) * Sq(s, s) * ground_state(s);
			
			cx_double alfa = cdot(cx_vec(ground_state, vec(N, fill::zeros)), SqSq_GS);
			Build_Lanczos_Hamil(Sq * ground_state / sqrt(alfa), M);
		//------------------------
		double SpinFactor = 0;
		for (double omega = 0; omega < pi; omega += 0.01) { // Calculate S(q,w)
			cx_double z = omega + 1i * 0.01 + Lanczos_eigenVal(0);
			cx_double Continous_Fraction = z - H_Lanczos(M-1, M-1);
			cx_double SqSq = 0;
			for (int m = M - 2; m >= 0; m--) { // m = n, because only for E_m > E_n we have constribution - eigenvalues are sorted and omega>0
				Continous_Fraction = z - H_Lanczos(m, m) - H_Lanczos(m, m + 1) * H_Lanczos(m, m + 1) / Continous_Fraction;
			}
			SpinFactor = -1 / pi * imag(alfa / Continous_Fraction);
			SpinFactorFile << q << "\t\t" << omega << "\t\t" << SpinFactor << endl;
		}
		std::cout << q << endl;
		Sq.~cx_mat();
		SqSq_GS.~cx_vec();
		ground_state.~vec();
	}
	SpinFactorFile.close();
}

//Lanczos Static Quantitites
void Heisenberg1D::Lanczos_Heat_capacity_Randomized(int M) {
	double dT = 0.01, T = dT;
	int R = 20;
	ofstream file("Heat_Capacity_Lanczos_M=" + to_string(M) + "_R=" + to_string(R) + ".txt");

	srand(time(NULL));
	while (T <= 5) {

		double Partition_Function = 0;
		double overlap = 0;
		double energy_av = 0, energy2_av = 0;

		for (int r = 1; r <= R; r++) { // sum over random vectors - different Lanczos procedures
			vec random_vector(N);
			for (int j = 0; j < N; j++)
				random_vector(j) = static_cast<double>(rand()) / (RAND_MAX + 0.0) - 0.5;
			
			Build_Lanczos_Hamil(cx_vec(random_vector, vec(N,fill::zeros)), M); //Lanczos procedure with |r> - random vector
			
			this->Lanczos_eigenVal = vec(M);
			mat eigenVec(M, M);
			eig_sym(Lanczos_eigenVal, eigenVec, H_Lanczos); // no eigenvectors for now
			this->Lanczos_eigenVec = mat(N, M); // m eigenvectors of lanczos iteration
			Lanczos_eigenVec = Krylov_space * eigenVec;
			
			eigenVec.~Mat();

			for (int m = 0; m < M; m++) {
				for (int j = 0; j < Lanczos_eigenVal.size(); j++) {
					overlap = abs(dot(random_vector, Lanczos_eigenVec.col(j)));
					overlap *= overlap;
					Partition_Function += overlap * std::exp(-Lanczos_eigenVal(j) / T); //partition function(T)
					energy_av += overlap*Lanczos_eigenVal(j) * std::exp(-Lanczos_eigenVal(j) / T); //average energy(T)
					energy2_av += overlap*Lanczos_eigenVal(j) * Lanczos_eigenVal(j) * std::exp(-Lanczos_eigenVal(j) / T); //average energy^2(T)
				}
			}
		}
		Partition_Function *= (double)N / (double)R;
		energy_av = energy_av / Partition_Function * (double)N / (double)R;
		energy2_av = energy2_av / Partition_Function *(double)N / (double)R;
		double heat_capacity_lancz = (energy2_av - energy_av * energy_av) / T / T / (L + 0.0);

		file << T << "\t\t" << heat_capacity_lancz << endl;
		T += dT;
	}
	file.close();
}



void Heisenberg1D::spinStructureFactor_FFT(double temperature) {
	double time = 0;
	double dt = 0.01;
	ofstream file("SpinFactor_FFT.txt");

	vector<double> time_vec;
	while (time <= 2 * L) {
		time_vec.push_back(time);
		time += dt;
	}
	int nloop = time_vec.size();

	mat I = eye(N, N);
	cx_mat temp_mat(H - I*eigenvalues(0), mat(N,N,fill::zeros));
	cx_vec psi_0(eigenvectors.col(0), vec(N, fill::zeros));
	psi_0 = exp(-eigenvalues(0) / temperature / 2.0) * psi_0;
	vector<vector<cx_double>> correlator(L);
	for (int r = 0; r < L; r++) { //corelation through whole chain
		correlator[r] = vector<cx_double>(nloop);
#pragma omp parallel for shared(temp_mat, time_vec, psi_0, nloop, correlator)
		for(int t = 0; t < nloop; t++) {
			cx_mat Sz0(N,N,fill::zeros), Szr(N,N,fill::zeros);
			cx_mat U(N, N, fill::zeros);
			U = expmat(-1i * time_vec[t] * temp_mat); //unitary evolution
			for (int j = 0; j < N; j++) {
				vector<int> vect = int_to_binary(mapping_inv(j), L);
				Sz0(j, j) = (double)vect[0] - 0.5;
				Szr(j, j) = (double)vect[r] - 0.5;
				vect.~vector();
			}
			correlator[r][t] = cdot(psi_0, Szr * U * Sz0 * psi_0); //C(r,t) = <psi_0|Sz^r * exp(-i*(H-E_0)*t ) * Sz^0|psi_0> - interaction picture
			//file << r << "\t\t" << time_vec[t] << "\t\t" << real(correlator[r][t]) << endl;
			Sz0.~Mat(); Szr.~Mat();
			U.~Mat();
		}
		std::cout << r << endl;
	}
	psi_0.~cx_vec(); I.~Mat(); temp_mat.~cx_mat();
	vector<vector<cx_double>> shit(nloop);
	ofstream file2("temp.txt");
//#pragma omp parallel for shared(temp_mat, time_vec, nloop, correlator, shit, file2)
	for (int t = 0; t < nloop; t++) {
		shit[t] = vector<cx_double>(L + 1);
		for (int k = 0; k <= L; k++) {
			shit[t][k] = 0;
			double q = 2 * pi / L * k;
			for (int r = 0; r < L; r++) {
				shit[t][k] += exp(1i * q * (r + 0.0)) * correlator[r][t];
			}
		}
		file2 << time_vec[t] << "\t\t" << real(shit[t][L/4 + 1]) << endl; // pi/2 plot
	}
	file2.close();
	/*for (int k = 0; k <= L; k++)
		for (int t = 0; t < nloop; t++)
			file << 2 * pi / L * k << "\t\t" << time_vec[t] << "\t\t" << real(shit[t][k]) << endl;*/

	shit.~vector();
	std::cout << "Time to FFT dis shit" << endl;
	//FFT of correlatorvector<double> omega_vec;
		double omega = 0, domega = dt;
		vector<double> omega_vec;
		while (omega <= eigenvalues(N - 1) - eigenvalues(0)) {
			omega_vec.push_back(omega);
			omega += domega;
		}
		nloop = omega_vec.size();

		vector<vector<cx_double>> SpinFactor(L + 1);
		for (int j = 0; j < L + 1; j++) {
			double q = 2 * pi / L * j;
			SpinFactor[j] = vector<cx_double>(nloop);
#pragma omp parallel for shared(omega_vec, nloop, time_vec, correlator, SpinFactor) //  check if HERE worsk parallel !!!!
			for (int w = 0; w < nloop; w++) {
				SpinFactor[j][w] = 0;
				for (int k = 0; k < L; k++) {
					cx_double temp = (exp(1i * time_vec[0] * omega_vec[w]) * correlator[k][0] + exp(1i * time_vec[time_vec.size() - 1] * omega_vec[w]) * correlator[k][time_vec.size()-1]) / 2.0 * dt;
					for (int t = 1; t < time_vec.size() - 1; t++) {
						temp += exp(1i * time_vec[t] * omega_vec[w]) * correlator[k][t] * dt;
					}
					SpinFactor[j][w] += exp(1i * q * (k + 0.0)) * temp;
				}

				file << q << "\t\t" << omega_vec[w] << "\t\t" << abs(SpinFactor[j][w]) << endl;
			}
			std::cout << j << endl;
		}
	//-----------------------------
	file.close();
	time_vec.~vector();
	omega_vec.~vector();
	correlator.~vector();
	SpinFactor.~vector();
}
//----------------------------------------------------------------------------------------------
//----------------------------------------Rest methods & functions------------------------------
//----------------------------------------------------------------------------------------------

//Getting private fields for usage outside the class
mat Heisenberg1D::get_hamil() {
	return this->H;
}
//----------------------------------------------------

//Factorial!!! - screw tgamma function, is shit
long long int factorial(int n){
	if (n > 1) return n * factorial(n - 1);
	else return 1;
}

// Binomial: n po k = n!/( k!*(n-k)! )
long long int Binomial(int n, int k) {
	return factorial(n) / factorial(k) / factorial(n - k);
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

//---------------------------------Main program to compute result-------------------------------
void Compute_subblocks(int L, double J, double Delta, string del) {	

	Heisenberg1D Object(L, L/2, J, Delta);
	Object.Hamiltonian();

	//Diagonalization
		Object.Diagonalization();
		std::cout << "\tHamiltonian L=" << L << " diagonalized!!" << endl;
	//----------------
	ofstream energies("Eigenvalues_Lanczos" + del + ".txt");
	stringstream ss;
	for (int m = 0; m <= 200; m += 1) {
		if (m >= 1) {
			Object.Lanczos_Diagonalization(m);
			Object.Lanczos_convergence(m, ss);
		}
	}
	energies << ss.str();
	energies.close();


	Object.~Heisenberg1D();
}

//---------------------------------Main program to compute result-------------------------------
void Compute_total(int L, double J, double Delta, string del) {

	Heisenberg1D Object(L, J, Delta);
	Object.Hamiltonian();
	
	//Diagonalization
		Object.Diagonalization();
		std::cout << "\tHamiltonian L=" << L << " diagonalized!!" << endl;
	//---------------

		ofstream savefile("HeatCapacity_ED_L="+to_string(L)+"_Delta=" + del + ".txt");
		Object.Heat_Capacity(0.01, 0.01, 5, savefile);
		savefile.close();
		/*ofstream energies("Eigenvalues_Lanczos_Delta=" + del + ".txt");
		stringstream ss;
		for (int m = 0; m <= 150; m += 1) {
			if (m >= 1) {
				Object.Lanczos_Diagonalization(m);
				Object.Lanczos_convergence(m, ss);
			}
		}
		energies << ss.str();
		energies.close();
		Object.Lanczos_gs_energy();*/

		int m = 100;
		Object.Lanczos_Diagonalization(m);
		Object.Lanczos_Heat_capacity_Randomized(m);
		//Object.SSF_Lanczos();

		/*for (int dim = 6; dim <= 12; dim++) {
			Heisenberg1D Hamil(dim, 1, 0);
			Hamil.Hamiltonian();
			Hamil.Diagonalization();
			int m_end;
			if (dim <= 8) m_end = pow(2, dim);
			else m_end = 250; // Binomial(dim, dim / 2);
			for (int m = 0; m <= m_end; m += 2) {
				if (m >= 1) {
					Hamil.Lanczos_Diagonalization(m);
					//Object.Lanczos_convergence(m, ss);
				}
			}
			Hamil.Lanczos_gs_energy();
			Hamil.~Heisenberg1D();
		}*/

	Object.~Heisenberg1D();
}






//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//---------------------------------------------Old definitions----------------------------------

//Time evolution, magentization
		/*Object.Magnetization_evolution(del);
		double domain_wall_temperature = Object.get_domain_wall_state_temp(0.001, 3 * L, del);
		std::cout << domain_wall_temperature << endl;
		Object.Overlap(5 * L, domain_wall_temperature, del);*/

		//double h = Object.get_average_gap();
		//cout << "Average energy gap:\t" << h << endl;
		/*double temperature = 0.0; //min 1./160.
			std::cout << "\n Temperature :" << temperature << endl << endl;
			//Object.spinStructureFactor_lin_responce(temperature);*/
			//Object.SSF_magnetic_field(h);
			//Object.spinStructureFactor_T_0();

		//Spin sctructure Factor
			//double temperature = 0.0625;
			//Object.spinStructureFactor_lin_responce(temperature);

	//Energy Gap, heat_capacity etc.
		/*
		ofstream fileC, fileGapeven, fileGapodd, fileEigenvalues;
		fileC.open("C_L=" + to_string(L) + "_Delta=" + del + ".txt");
		fileGapeven.open("Egap_even__Delta=" + del + ".txt"); //energy gap for even chain length
		fileGapodd.open("Egap_odd__Delta=" + del + ".txt"); //energy gap for odd chain length
		fileEigenvalues.open("Eigenvalues_Delta=" + del + ".txt", ios_base::app);
		//my_mutex.lock(); //Blocking for other threads
		Object.print_eigenvalues(fileEigenvalues);
		//my_mutex.unlock(); //unlocking
		if (clear_H == false) {
			std::cout << "Eigenvalues correct?:\n";
			Col<int> check = Object.check_eigenvalues();
			std::cout << check;
			std::cout << "\nOh yeah";
		}

		double Gap = Object.get_energy_gap();
		double invL = 1. / L;

		//my_mutex.lock();//Blocking for other threads
		if (L % 2 == 0) {
			fileGapeven << fixed << setprecision(16) << invL << ", " << Gap << endl;
		}
		else {
			fileGapodd << fixed << setprecision(16) << invL << ", " << Gap << endl;
		}
		//my_mutex.unlock(); //unlocking

		Object.Heat_Capacity(0.01, 0.01, 4., fileC); // dT, T_0, T_end, file to save data

		fileC.close();
		fileEigenvalues.close();
		fileGapeven.close();
		fileGapodd.close();*/


//Print hamiltonian
	/*//Object.print_Hamiltonian_total(file2);
	//create a locale that has the ctype category  copied from the "en_US.UTF-8"
		locale loc = locale(locale(), "en_US.UTF8", locale::ctype);
		wofstream file2;
		file2.open("Hamiltonian_total_L=" + to_string(L) + "_Delta=" + del + ".txt");
		file2.imbue(loc); // adding the locale to the stream
	//-------------------------------------------------------
	file2.close();*/
//----------

// Generating total Hamiltonian
/*void Heisenberg1D::Hamiltonian_total() {
	int s_i, s_j; //j=i+1
	int idx; //indices equivalent to spin-filp due to kinetic term
	for (int k = 0; k < N; k++) {
		base_vector = int_to_binary(k, L);
		for (int j = 0; j < L - 1; j++) {
			s_i = base_vector[j];
			s_j = base_vector[j + 1];
			H(k, k) += Delta * (s_i - 0.5) * (s_j - 0.5); // "potential" term - diagonal part, spin s_z_i = eigenvector[i] - 1/2
			if (s_i == 0 && s_j == 1) { // S_i^+ S_i+1^-
				/*base_vector(j) = 1; //spin filp
				base_vector(j + 1) = 0;
				idx = binary_to_int(base_vector); //getting index of the resulting eigenvector
				base_vector(j) = 0; //in order to preserve eigenvector for next iteration
				base_vector(j + 1) = 1;
				idx = k - std::pow(2, j); //  we now what changed, therefore the decimal representation has to change by +-2^j
				H(idx, k) = J / 2.;
			}
			if (s_i == 1 && s_j == 0) { // S_i^- S_i+1^+
				/*base_vector(j) = 0; //spin filp (inverse to the previous)
				base_vector(j + 1) = 1;
				idx = binary_to_int(base_vector);
				base_vector(j) = 1; //in order to preserve eigenvector for next iteration
				base_vector(j + 1) = 0;
				idx = k + std::pow(2, j);
				H(idx, k) = J / 2.;
			}
		}
	}
}*/
//----------------------------------------------------

//Calculates the magnetization on each site over time evolution -> output to file
/*void Heisenberg1D::Magnetization_evolution_total(string del) {

	double dt = 0.05;
	double time = 0;

	cx_vec wavefunction(N, fill::zeros);
	int idx = binary_to_int(first_vec);
	wavefunction(idx) = 1; //first base_vector is 0011 = first_vec
	ofstream file_evolution("Unitary evolution_Delta=" + del + ".txt");
	double S_z;
	while (time <= 2 * L) { // +dt because new wavefunction is at loop end
		for (int k = 0; k < L; k++) {
			S_z = 0;
			for (int j = 0; j < N; j++) {
				vector<int> vect = int_to_binary(mapping_inv(j), L); //jth base_vector
				S_z = S_z + (vect[k] - 0.5) * abs(wavefunction(j)); //otput should be real, but C++ forbids such typecast
				vect.~vector();
			}
			file_evolution << k << "\t" << time << "\t" << S_z << endl;
		}
		time += dt;
		//wavefunction = Time_evolution_stationary(wavefunction, time);
		wavefunction = Time_evolution(wavefunction, dt);
	}
	wavefunction.clear();
	file_evolution.close();
}*/
//-----------------------------------------------------------------------------