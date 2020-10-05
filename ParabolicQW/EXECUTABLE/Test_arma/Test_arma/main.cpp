#include <iostream>
#include <armadillo>
double m_0 = 9.10938356; //electron mass in [kg] // 10^-31
double hbar = 6.582119514; //planck constant in [eV*s] //10^-16

using namespace arma;
using namespace std;

int main(int argc, const char* argv) {
	int power = 1; // 10^-31 * (10^15)^2 * 10^(-18) / 10^(-19) = 1
	double m_0 = 9.10938356; //electron mass in [kg] // 10^-31
	double e_c = 1.60211; // elementary electric charge //10^-19

	double omega = 1; //in pHz 10^15
	double x_0 = 0, x_1 = -5, x_2 = 5, dx = 0.0001; // in nm
	int num_of_eigenval;
	ifstream file2("parameters.dat");
	file2 >> x_0 >> x_1 >> x_2 >> dx >> omega >> num_of_eigenval;
	file2.close();

	// Build Quantum Well potential with same x-grid as Hamil:
	std::vector<double> potential(static_cast<int>((x_2 - x_1) / dx + 1));
	for (int k = 0; k < potential.size(); k++) {
		double x = x_1 + k * dx;
		potential[k] = 1. / 2. / e_c * m_0 * power * omega * omega * (x - x_0) * (x - x_0);
	}

	int N = potential.size();
	power = pow(10, 2); // 10^-31 * (10^-8)^2 / 10^(-32) / 10^(-19)
	double alfa = (2 * m_0 * dx * dx) / (hbar * hbar) * power / e_c;

	//For finding minimal/maximal value by Gerschgorin Cirle Theorem
	double tmp = 2. + 1. / 2. * m_0 * omega * omega * (x_1 - x_0) * (x_1 - x_0) * alfa - 1.; // H[0][0] - 1
	double temp2 = tmp + 2.; // H[0][0] + 1

	mat H(N, N, fill::zeros);
	for (int j = 0; j < N - 1; j++) {
		H(j, j + 1) = -1.;
		H(j + 1, j) = -1.;
		H(j, j) = 2. + potential[j] * alfa;
	}
	H(N - 1, N - 1) = 2. + potential[N - 1] * alfa;
	ofstream file("eigenVal.dat"), file3("eigenVec.dat");
	vec eigenVal(N);
	mat eigenVec(N, N);
	eig_sym(eigenVal, eigenVec, H);
	eigenVal = eigenVal / alfa;
	for (int j = 0; j < num_of_eigenval; j++) {
		file << eigenVal(j) << endl;
		for (int k = 0; k < N; k++)
			file3 << eigenVec.col(j)(k) << endl;
	}
	file.close();
	file3.close();

	H.~Mat(); eigenVal.~vec(); eigenVec.~Mat();

	return 0;

}