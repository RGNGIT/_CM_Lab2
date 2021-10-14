#include <iostream>
#include <conio.h>
#include <string>
#include "gausslib.h"

std::string * SLAU2String(double ** dptrMatrix, double * ptrEqual) {
	std::string * stringify = new std::string[N];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			stringify[i] += dptrMatrix[i][j] >= 0 ? ((j > 0 ? " + " : "")
				+ std::to_string(dptrMatrix[i][j])) : (" - "
			    + std::to_string(dptrMatrix[i][j] * -1));
			stringify[i] += "x" + std::to_string(j + 1);
		}
		stringify[i] += " = " + std::to_string(ptrEqual[i]);
	}
	return stringify;
}

double ** dptrMatrix = new double * [N] {
	new double[N] {-1, 0.28, -0.17, 0.06},
	new double[N] {0.52, -1, 0.12, 0.17},
	new double[N] {0.17, -0.18, -0.79, 0},
	new double[N] {0.11, 0.22, 0.03, -0.95}
};
double * ptrEqual = new double[N] {
	-21, 117, 0.81, -0.72
};

int main(void) {
	setlocale(LC_ALL, "");
	std::cout << "СЛАУ:" << std::endl;
	for (int i = 0; i < N; i++) {
		std::cout << SLAU2String(copyMatrix(dptrMatrix), copyVector(ptrEqual))[i] << std::endl;
	}
	std::cout << std::endl << "Корни:" << std::endl;
	for (int i = 0; i < N; i++) {
		std::cout << "x" << (i + 1) << " == " << gaussSolve(copyMatrix(dptrMatrix), copyVector(ptrEqual))[i] << std::endl;
	}
	std::cout << std::endl << "Определитель: " << gaussDet(copyMatrix(dptrMatrix)) << std::endl;
	std::cout << std::endl << "Обратная матрица: ";
	for (int i = 0; i < N; i++) {
		std::cout << std::endl;
		for (int j = 0; j < N; j++) {
			std::cout << gaussReverse(copyMatrix(dptrMatrix))[i][j] << "|";
		}
	}
	_getch();
	return 0;
}