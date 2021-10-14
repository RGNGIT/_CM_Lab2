#include <cmath>
#define N 4

double ** dptrSolvedMatrix;

double ** copyMatrix(double ** dptrMatrix) {
	double ** clone = new double * [N];
	for (int i = 0; i < N; i++) {
		clone[i] = new double[N];
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			clone[i][j] = dptrMatrix[i][j];
		}
	}
	return clone;
}

double * copyVector(double * ptrVector) {
	double * clone = new double[N];
	for (int i = 0; i < N; i++) {
		clone[i] = ptrVector[i];
	}
	return clone;
}

double * gaussSolve(double ** dptrMatrix, double * ptrEqual) {
	double * x = new double[N];
	double R;
	for (int q = 0; q < N; q++) {
		R = 1 / dptrMatrix[q][q];
		dptrMatrix[q][q] = 1;
		for (int j = q + 1; j < N; j++) {
			dptrMatrix[q][j] *= R;
		}
		ptrEqual[q] *= R;
		for (int k = q + 1; k < N; k++) {
			R = dptrMatrix[k][q];
			dptrMatrix[k][q] = 0;
			for (int j = q + 1; j < N; j++) {
				dptrMatrix[k][j] = dptrMatrix[k][j] - dptrMatrix[q][j] * R;
			}
			ptrEqual[k] = ptrEqual[k] - ptrEqual[q] * R;
		}
	}
	dptrSolvedMatrix = dptrMatrix;
	for (int q = N - 1; q >= 0; q--) {
		R = ptrEqual[q];
		for (int j = q + 1; j < N; j++) {
			R -= dptrMatrix[q][j] * x[j];
		}
		x[q] = R;
	}
	return x;
}

double gaussDet(double ** dptrMatrix) {
	double tmp, det = 1;
	for (int k = 0; k < N - 1; k++) {
		for (int i = k + 1; i < N; i++) {
			tmp = -dptrMatrix[i][k] / dptrMatrix[k][k];
			for (int j = 0; j < N; j++) {
				dptrMatrix[i][j] += dptrMatrix[k][j] * tmp;
			}
		}
	}
	for (int i = 0; i < N; i++) {
		det *= dptrMatrix[i][i];
	}
	return det;
}

double ** gaussReverse(double ** dptrMatrix) {
	double ** res = new double * [N];
	double * y = new double[N];
	double * gauss;
	int i, j, k;
	for (int i = 0; i < N; i++) {
		res[i] = new double[N];
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				y[j] = 1;
			}
			else {
				y[j] = 0;
			}
		}
		gauss = gaussSolve(copyMatrix(dptrMatrix), y);
		for (int k = 0; k < N; k++) {
			res[k][i] = gauss[k];
		}
	}
	return res;
}


