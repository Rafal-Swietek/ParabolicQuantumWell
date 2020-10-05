#include "Hamiltionian.h"
double m_0 = 9.10938356; //electron mass in [kg] // 10^-31
double hbar = 6.582119569; //planck constant in [eV*s] //10^-16
double e_c = 1.60211; // elementary electric charge //10^-19
std::ofstream logfile("log.txt"); //logfile

// Contructor - destructor
Hamiltonian::Hamiltonian(double omega, double QW_pos, double x_1, double x_2, double dx){
	this->N = static_cast<int>((x_2 - x_1) / dx + 1);
	this->x_1 = x_1; 
	this->x_2 = x_2; 
	this->dx = dx;
	// Tridiagonal symmetric matrix:
		this->H = std::vector<std::vector<double>>(2);
		H[0] = std::vector<double>(N); //diagonal part
		H[1] = std::vector<double>(N); //off-diagonal
	//this->eigenVal = std::vector<double>(N);
	this->omega = omega;
	this->QW_position = QW_pos;
	logfile << "Welcome to log!\n\t Below you'll see any bug and progress of your programm:\n" << std::endl;
}
Hamiltonian::~Hamiltonian(){
	H.~vector();
	eigenVal.~vector();
}


// Building Hamiltonian matrix for parabolic potential
void Hamiltonian::Build_Hamiltonian(std::vector<double> potential) {

	logfile << "Start Hamiltonian Build" << std::endl;
	logfile << "H = \n diagonal\toff-diagonal" << std::endl;
	int power = std::pow(10,2); // 10^-31 * (10^-9)^2 / 10^(-32) / 10^(-19)
	this->alfa = (2 * m_0 * dx * dx) / (hbar * hbar) * power / e_c;

	//For finding minimal/maximal value by Gerschgorin Cirle Theorem
		double tmp = 2. + 1. / 2. * m_0 * omega * omega * (x_1 - QW_position) * (x_1 - QW_position) * alfa - 1.; // H[0][0] - 1
		double temp2 = tmp + 2.; // H[0][0] + 1

	for (int j = 0; j < N - 1; j++) {
		H[1][j] = -1.;
		H[0][j] = 2. + potential[j] * alfa;
		if (H[0][j] - 2. < tmp) //minimal eigenvalue estimation
			tmp = H[0][j] - 2.;
		if (H[0][j] + 2. > temp2) //maximal eigenvalue estimation
			temp2 = H[0][j] + 2.;
		logfile << H[0][j] << "\t\t" << H[1][j] << std::endl;
	}
	H[1][N - 1] = 0;
	H[0][N - 1] = 2 + potential[N - 1] * alfa;
	logfile << H[0][N - 1] << std::endl;
	if (H[0][N - 1] - 1. < tmp)
		tmp = H[N - 1][N - 1] - 1.;
	if (H[0][N - 1] + 1. > temp2)
		temp2 = H[N - 1][N - 1] + 1.;
	E_min = tmp; E_max = temp2;
	logfile << "Hamiltonian Built successfully!\n" << std::endl;
	logfile << "Eigenvalue extrema:\t" << E_min << "-" << E_max << std::endl;
}

// Diagonalizing matrix with Sturm sequence and Ehrlich-Albert root finding
void Hamiltonian::Strum_Ehrlich_Diagonalization(int num_of_eigenVal){
	double a = E_min, b = E_max;
	double tol = dx;
	int ID = 0; // some index
	std::vector<double> solution(num_of_eigenVal);
	for (int k = 0; k < num_of_eigenVal / 2; k++) { //evenly distributed initial guess using x=0 mirror symmetry
		if (k < num_of_eigenVal / 2)
			solution[k] = b - (b - a) / (num_of_eigenVal / 2. - 1) * k;
		else if (k = num_of_eigenVal / 2)
			solution[k] = a;
		else
			solution[k] = a + (b - a) / (num_of_eigenVal / 2. - 1) * k;
	}
	/*for (int k = num_of_eigenVal / 2; k < num_of_eigenVal; k++) {
		solution[k] = omega * (k - num_of_eigenVal / 2. + 1./2.) + 0.5 * pow(-1,k+1);
		solution[k - num_of_eigenVal / 2] = solution[k];
	}*/
	// Calculate roots of characteristic equation with Ehrlich-Albert algorithm
	logfile << "Calculate eigenvalues" << std::endl;
	while (true) {
		std::vector<double> temp(num_of_eigenVal);
		temp = solution;
//#pragma omp parallel for shared(solution)
		for (int i = 0; i < num_of_eigenVal; i++) {
			double pol_k, pol_der_k;
			pol_k = Sturm_sequence(solution[i], N);
			pol_der_k = Sturm_sequence_derivative(solution[i], N);
			double sum = 0;
			//Calculate the sum
				for (int j = 0; j < num_of_eigenVal; j++)
					if( j != i )
						sum += 1 / (solution[i] - solution[j]);
			solution[i] = solution[i] - (pol_k/pol_der_k) / (1 - pol_k / pol_der_k * sum);
		}
		// Check convergence
		double modul1 = 0, modul2 = 0;
		for (int p = 0; p < num_of_eigenVal; p++) {
			modul1 += solution[p] * solution[p];
			modul2 += temp[p] * temp[p];
		}
		if (std::sqrt(fabs(modul1 - modul2)) < tol) break;
		temp.~vector();
		ID++;
		logfile << "Step:\t" << ID << std::endl;
	}
	this->eigenVal = std::vector<double>(num_of_eigenVal);
	for (int k = 0; k < num_of_eigenVal; k++)
		eigenVal[k] = solution[k] / alfa; // store calculated roots, which are eigenvalues // E = \alfa*eigenVal
	solution.~vector();
	logfile << "Eigenvalues calculated!" << std::endl;

	//sort(eigenVal.begin(), eigenVal.end()); // sort eigenvalues in ascending order

	// Bisection
	/*while (i < num) {
		c = (b - a) / 2.0;
		if (fabs(Sturm_sequence(H, c, N)) < tol) {
			eigenVal.push_back(c);
			break;
		}
		if (Sturm_sequence(H, a, N) * Sturm_sequence(H, c, N) < 0) b = c;
		else if (Sturm_sequence(H, c, N) * Sturm_sequence(H, b, N) < 0) a = c;
		i++;
	}*/
}

// Finding given number of eigenvalues using Sturm sequence
void Hamiltonian::Sturm_Diagonalization(int num_of_eigenVal) {
	this->eigenVal = std::vector<double>(0);
	logfile << "Calculate eigenvalues" << std::endl;
	for (int k = 0; k < num_of_eigenVal; k++) {
	repeat:
		//for (int k = 0; k < N; k++) {
		std::vector<double>interval = find_set_containing_eigenVal(k); // a,b,poly(a), poly(b)
		if (interval[2] == 0)
			interval[1] = interval[0];

		while (fabs(interval[1] - interval[0]) > dx) {
			double mid_point = (interval[0] + interval[1]) / 2.;
			std::vector<double> temp = Sturm_sequence_sign_change(mid_point);
			if (temp[0] == 0) {
				interval[0] = mid_point; interval[1] = mid_point;
			}
			if (temp[0] * interval[3] < 0) {
				interval[0] = mid_point;
				interval[2] = temp[0];
			}
			else {
				interval[1] = mid_point;
				interval[3] = temp[0];
			}
		}
		double E = (interval[1] - interval[0]) / 2. / alfa;
		/*bool result = false;
		for (int m = 0; m < eigenVal.size(); m++) {
			result = result || (E == eigenVal[m]);
		}
		if (result == true)
			goto repeat;
		else*/
		eigenVal.push_back(E);
	}
	logfile << "Eigenvalues calculated! E = :" << std::endl;
	sort(eigenVal.begin(), eigenVal.end()); // sort eigenvalues in ascending order
	//eigenVal.erase(std::unique(eigenVal.begin(), eigenVal.end()), eigenVal.end()); //delete duplicates from vector, i.e. remove degeneracy
	for (int k = 0; k < eigenVal.size(); k++) {
		logfile << eigenVal[k] << std::endl;
	}
}

std::vector<double> Hamiltonian::find_set_containing_eigenVal(int k) {
	std::vector<double> interval(4);
	double tol = dx;
	double b = E_max, a = E_min;
	std::vector<double> pol_k_min = Sturm_sequence_sign_change(a);
	std::vector<double> pol_k_max = Sturm_sequence_sign_change(b);
	while ((pol_k_max[1] > k) || (pol_k_min[1] < k - 1)) {
		double mid_point = (a + b) / 2.;
		std::vector<double> mid_vec = Sturm_sequence_sign_change(mid_point);
		if (mid_vec[1] >= k) {
			pol_k_max = mid_vec;
			b = mid_point;
		}
		else {
			pol_k_min = mid_vec;
			a = mid_point;
		}	
	}
	interval[0] = a; interval[1] = b; 
	interval[2] = pol_k_min[0]; interval[3] = pol_k_max[0];
	return interval;
}

std::vector<double> Hamiltonian::Sturm_sequence_sign_change(double x) {
	std::vector<double> poly_value_sign(2);
	int sign = 0;
	double pn1 = 1, pn = H[0][0] - x;
	if ( (pn1 * pn < 0) || pn == 0) sign = 1;
	for (int j = 1; j < N; j++) {
		double pn2 = pn1; pn1 = pn;
		pn = (H[0][j] - x) * pn1 - H[1][j] * H[1][j] * pn2;
		if ((pn1 * pn < 0) || pn == 0) sign++;
	}
	poly_value_sign[0] = pn;
	poly_value_sign[1] = sign;
	return poly_value_sign;
}

double Hamiltonian::Sturm_sequence(double x, int k) {
	if (k == 0) return 1;
	else if (k == 1) return H[0][0] - x;
	else return (H[0][k - 1] - x) * Sturm_sequence(x, k - 1) - H[1][k - 1] * H[1][k - 1] * Sturm_sequence(x, k - 2);
}
double Hamiltonian::Sturm_sequence_derivative(double x, int k) {
	if (k == 0) return 0;
	else if (k == 1) return -1;
	else
		return (H[0][k - 1] - x) * Sturm_sequence_derivative(x, k - 1)\
		- H[1][k - 1] * H[1][k - 1] * Sturm_sequence_derivative(x, k - 2) - Sturm_sequence(x, k - 1);
}


std::vector<double> Hamiltonian::Calculate_eigenvector_INVERSE_POWER(int n){
	std::vector<double> wavefunction(N);
	// matrix for gauss elimination - 
	//	matrix = H - E_n*eye(N)
		std::vector<double> a(N), b(N), c(N); //  - tridiagonal matrix => 3 vectors
													  
	//	Generate random first vector
		double norm = 0;
		srand(time(NULL));
		for (int j = 0; j < N; j++){
			wavefunction[j] = static_cast<double>(rand()) / (RAND_MAX + 0.0);
			norm += std::pow(wavefunction[j], 2);
			if (j < N - 1) 
				c[j] = H[1][j];
			if (j > 0)
				a[j] = H[1][j];
			b[j] = H[0][j] - eigenVal[n];
		}
		for (int k = 0; k < wavefunction.size(); k++) {
			//wavefunction[k] /= norm;
		}
	double difference = DBL_MAX; 
	double dot_prod;
	double mu = eigenVal[n];
	double E = mu + 0.5;
	// Inverse power method termineated, when || z_k - z_k_1 || < tolerance - z being the eigenvector
	logfile << "Find given eigenvector" << std::endl;
	int iter = 0;
	while (true) {
		std::vector<double> solution = wavefunction; // to edit wavefunction - d vector in Thomas algoritm
		std::vector<double> wafefunction_prev = wavefunction; // previous iteration result
		double E_prev = E;
		// Recreate the a,b,c vectors
			
		//Thomas algorithm for tridiagonal system of equation (E-H)*wavefunction = wavefunction_prev
			std::vector<double> a_tmp = a, b_tmp = b, c_tmp = c;
			for (int k = 1; k < N; k++) {
				double w = a_tmp[k] / b_tmp[k - 1];
				b_tmp[k] = b_tmp[k] - w * c_tmp[k - 1];
				solution[k] = solution[k] - w * solution[k - 1];
			}
			// Back substitution
				wavefunction[N - 1] = solution[N - 1] / b_tmp[N - 1];
				norm = std::pow(wavefunction[N - 1], 2);
				for (int k = N - 2; k >= 0; k--) {
					wavefunction[k] = (solution[k] - c_tmp[k] * wavefunction[k + 1]) / b_tmp[k];
					norm += std::pow(wavefunction[k], 2);
				}
			dot_prod = 0;
			//Wavefunction normalisation and energy calculation as: E=q^k * H * q^k - with q^k being wavefunction in k^th step
			E = wavefunction[0] * (wavefunction[0] * H[0][0] + wavefunction[1] * H[1][0]);
			for (int k = 0; k < wavefunction.size(); k++) {
				wavefunction[k] /= sqrt(norm);
				dot_prod += wafefunction_prev[k] * wavefunction[k];
				if ((k < N - 1) && (k > 0))
					E += wavefunction[k] * (wavefunction[k - 1] * H[1][k] + wavefunction[k]*H[0][k] + wavefunction[k + 1] * H[1][k]);
			}
			E += wavefunction[N-1] * (wavefunction[N - 2] * H[1][N-1] + wavefunction[N-1] * H[0][N-1]);
			logfile << E << std::endl;
		//E = mu + 1. / dot_prod;
		wafefunction_prev.~vector();
		solution.~vector();
		a_tmp.~vector(); b_tmp.~vector(); c_tmp.~vector();
		iter++;
		if (fabs(E - E_prev) < 1e-8) break;
		//if( fabs(dot_prod - 1.0) < 1e-12) break;
	}

	for (int k = 0; k < wavefunction.size(); k++) {
		wavefunction[k] /= (hbar * omega * 0.1);
		wavefunction[k] += eigenVal[n];
		logfile << wavefunction[k] << std::endl;
	}
	logfile << "Eigenvector found in " << iter << "iterations!!" << std::endl;
	a.~vector(); b.~vector(); c.~vector();
	return wavefunction;
}

std::vector<double> Hamiltonian::Calculate_eigenvector(int n) {
	logfile << "Find given eigenvector" << std::endl;

	double E_n = eigenVal[n] * (2 * m_0 * dx * dx) / (hbar * hbar) * pow(10, 2) / e_c;
	std::vector<double> wavefunction(N), delta_minus(N), delta_plus(N);
	wavefunction[0] = 0; wavefunction[N - 1] = 0;

	delta_minus[0] = 1.0 / (H[0][0] - E_n);
	delta_plus[N - 1] = 1.0 / (H[0][N - 1] - E_n);

	delta_minus[1] = 1.0 / (H[0][1] - E_n - H[1][0] * H[1][0] * delta_minus[0]);
	for (int k = 1; k <= N - 1; k++) {
		delta_minus[k] = 1.0 / (H[0][k] - E_n - H[1][k - 1] * H[1][k - 1] * delta_minus[k - 1]);
		int k2 = N - k - 1;
		delta_plus[k2] = 1.0 / (H[0][k2] - E_n - H[1][k2] * H[1][k2] * delta_plus[k2 + 1]);
	}

	double minimum = DBL_MAX;
	int index = 0;
	for (int k = 1; k < N - 1; k++) {
		double tmp = fabs(H[1][k - 1] * H[1][k - 1] * delta_minus[k - 1] + (H[0][k] - E_n) + H[1][k] * H[1][k] * delta_plus[k + 1]);
		if (tmp <= minimum) {
			minimum = tmp; 
			index = k;
		}
	}
	index = N / 2 + 1;
	logfile << "Index:" << index << " for energy = " << eigenVal[n] << std::endl;
	// Create wavefunction
	wavefunction[index] = 1;
	double norm = 1;
	for (int k = index - 1; k >= 0; k--) {
		wavefunction[k] = -delta_minus[k] * H[1][k] * wavefunction[k + 1];
		norm += wavefunction[k] * wavefunction[k];
	}
	for (int k = index + 1; k < N; k++) {
		wavefunction[k] = -delta_plus[k] * H[1][k - 1] * wavefunction[k - 1];
		norm += wavefunction[k] * wavefunction[k];
	}
	wavefunction = Analitycal_eigenvector(n);
	double maximal = 0;
	for (int k = 0; k < wavefunction.size(); k++) {
		//wavefunction[k] /= sqrt(norm);
		if ((wavefunction[k]) > maximal) {
			maximal = (wavefunction[k]);
			index = k;
		}
	}
	index -= static_cast<int>(0.2 / dx + 1) - 1;
	wavefunction = std::vector<double>(N);
	logfile << "Index:" << index << " for energy = " << eigenVal[n] << std::endl;
	// Create wavefunction
	wavefunction[index] = 0.1;
	norm = 1;
	for (int k = index - 1; k >= 0; k--) {
		wavefunction[k] = -delta_minus[k] * H[1][k] * wavefunction[k + 1];
		norm += wavefunction[k] * wavefunction[k];
	}
	for (int k = index + 1; k < N; k++) {
		wavefunction[k] = -delta_plus[k] * H[1][k - 1] * wavefunction[k - 1];
		norm += wavefunction[k] * wavefunction[k];
	}

	maximal = 0;
	for (int k = 0; k < N; k++) {
		wavefunction[k] /= sqrt(norm);
		if (fabs(wavefunction[k]) > maximal) maximal = fabs(wavefunction[k]);
	}

	// Normalize by energy Gap and shift to given eigenvalue
		for (int k = 0; k < N; k++) {
			wavefunction[k] /= maximal;
			if (hbar * omega * 0.1 < 1)
				wavefunction[k] *= (0.5 * hbar * omega * 0.1);
			else
				wavefunction[k] /= (0.5 * hbar * omega * 0.1);
			wavefunction[k] += eigenVal[n];
		}

	logfile << "Eigenvector found!!" << std::endl;
	return wavefunction;
}


std::vector<double> Hamiltonian::Calculate_eigenvector_Numerov(std::vector<double> potential, int n){

	logfile << "Find given eigenvector" << std::endl;
	this->alfa = (2 * m_0 * dx * dx) / (hbar * hbar) * pow(10, 2) / e_c;

	double E_n = eigenVal[n];
	std::vector<double> wavefunction(N), temp(N);
	/*wavefunction[0] = 0;
	wavefunction[1] = std::pow(-1, n);
	double norm = fabs(wavefunction[1] * wavefunction[1]);
	for (int n = 1; n < N - 1; n++) {
		double f_n = 1 + 1.0 / 12.0 * alfa * (E_n - potential[n]);
		double f_n_1 = 1 + 1.0 / 12.0 * alfa * (E_n - potential[n - 1]);
		double f_n_plus = 1 + 1.0 / 12.0 * alfa * (E_n - potential[n + 1]);
		wavefunction[n + 1] = ((12. - 10. * f_n) * wavefunction[n] - f_n_1 * wavefunction[n - 1]) / f_n_plus;
		norm += fabs(wavefunction[n] * wavefunction[n]);
	}*/
	// From left side Numerov
	double E = hbar * omega * 0.1 * (n + 0.5);
	double x1 = QW_position - sqrt(2 * E / (m_0 * omega * omega));
	int break_point = static_cast<int>( (x1 - x_1)/dx + 1.);
	wavefunction[0] = 1;
	wavefunction[1] = std::pow(-1, n);
	double norm = fabs(wavefunction[1] * wavefunction[1]);
	double max1 = -DBL_MAX, max2 = -DBL_MAX;
	for (int n = 1; n <= break_point; n++) {
		double f_n = 1 + 1.0 / 12.0 * alfa * (E_n - potential[n]);
		double f_n_1 = 1 + 1.0 / 12.0 * alfa * (E_n - potential[n - 1]);
		double f_n_plus = 1 + 1.0 / 12.0 * alfa * (E_n - potential[n + 1]);
		wavefunction[n + 1] = ((12. - 10. * f_n) * wavefunction[n] - f_n_1 * wavefunction[n - 1]) / f_n_plus;
		norm += fabs(wavefunction[n] * wavefunction[n]);
		if (fabs(wavefunction[n + 1]) > max1)
			max1 = fabs(wavefunction[n + 1]);
	}
	// From right side Numerov
	wavefunction[N - 1] = 0; 
	wavefunction[N - 2] = 0.1;
	norm += fabs(wavefunction[N - 2] * wavefunction[N - 2]);
	for (int n = N - 2; n > break_point; n--) {
		double f_n = 1 + 1.0 / 12.0 * alfa * (E_n - potential[n]);
		double f_n_1 = 1 + 1.0 / 12.0 * alfa * (E_n - potential[n - 1]);
		double f_n_plus = 1 + 1.0 / 12.0 * alfa * (E_n - potential[n + 1]);
		wavefunction[n - 1] = ((12. - 10. * f_n) * wavefunction[n] - f_n_plus * wavefunction[n + 1]) / f_n_1;
		norm += fabs(wavefunction[n] * wavefunction[n]);
		if (fabs(wavefunction[n - 1]) > max2) 
			max2 = fabs(wavefunction[n - 1]);
	}
	for (int n = N - 1; n >= break_point; n--) {
		wavefunction[n] *= fabs(max1 / max2);
	}

	double maximal = 0;
	for (int k = 0; k < N; k++) {
		wavefunction[k] += temp[k];
		wavefunction[k] /= sqrt(norm);
		if (wavefunction[k] > maximal) maximal = wavefunction[k];
	}

	for (int k = 0; k < wavefunction.size(); k++) {
		wavefunction[k] /= maximal * (hbar * omega * 0.1);
		wavefunction[k] += eigenVal[n];
		logfile << wavefunction[k] << std::endl;
	}

	return wavefunction;
}

std::vector<double> Hamiltonian::Analitycal_eigenvector(int n){
	std::vector<double> wavefunction(0);
	int power = sqrt(pow(10, 1));
	double B = sqrt(m_0 * omega / hbar / e_c) * power;
	double A = sqrt(factorial(n) / pow(2, n)) * sqrt(sqrt(m_0 * omega / hbar / 3.141592653));
	double norm = 0;
	double maximum = 0;
	if (n % 2 == 0) {
		double x = x_1 - 0.2;
		while(x<= x_2 + 0.2 + dx){
			double temp = 0;
			for (int l = 0; l <= n / 2; l++)
				temp += pow(-1, n / 2 - l) / factorial(2 * l) / factorial(n / 2 - l) * pow(2 * B * x, 2 * l);
			temp *= A * exp(-B * B * x * x / 2);
			if (temp > maximum) maximum = temp;
			wavefunction.push_back(temp);
			norm +=temp * temp;
			x += dx;
		}
	}
	else {
		double x = x_1 - 0.2;
		while (x <= x_2 + 0.2 + dx) {
			double temp = 0;
			for (int l = 0; l <= (n - 1) / 2; l++)
				temp += pow(-1, (n - 1) / 2 - l) / factorial(2 * l + 1) / factorial((n - 1) / 2 - l) * pow(2 * B * x, 2 * l + 1);
			temp *= A * exp(-B * B * x * x / 2);
			if (temp > maximum) maximum = temp;
			wavefunction.push_back(temp);
			norm += temp * temp;
			x += dx;
		}
	}
	for (int k = 0; k < wavefunction.size(); k++) {
		wavefunction[k] /= maximum;
		if(hbar * omega * 0.1 < 1)
			wavefunction[k] *= (0.5*hbar*omega*0.1);
		else
			wavefunction[k] /= (0.5 * hbar * omega * 0.1);
		wavefunction[k] += hbar * omega * 0.1 * (n + 0.5);
	}
	return wavefunction;
}

//Factorial!!! - screw tgamma function, is shit
long long int factorial(int n) {
	if (n > 1) return n * factorial(n - 1);
	else return 1;
}