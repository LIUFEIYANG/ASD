/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: Processer.cpp
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#include "Processer.h"
#include <time.h>
#include <iostream>
#include <math.h>
#include <limits>
#include <map>
using namespace std;

#ifndef BITDEPTH_12_TEST
#define BITDEPTH_12_TEST 0
#endif
#if BITDEPTH_12_TEST
const double threshold = 4094;
const double offset    = 4096;
#else
const double threshold = 32767;
const double offset    = 65536;
#endif

#define USE_GSL 0
#if USE_GSL
#include <gsl/gsl_eigen.h>
#endif

bool Processer::xcomparefunction(DistanceReference a, DistanceReference b)
{
	return a.distance < b.distance;
}

MatrixXd Processer::MatrixDRound(MatrixXd matrix)
{
	MatrixXd res = MatrixXd::Zero(matrix.rows(), matrix.cols());
	for (auto row = 0; row < matrix.rows(); row++)
	{
		for (auto col = 0; col < matrix.cols(); col++)
		{
			if (matrix(row, col) < 0)
			{
				res(row, col) = 0 - (double)((int)(0.5 - matrix(row, col)));
			}
			else
			{
				res(row, col) = (double)((int)(matrix(row, col) + 0.5));
			}
		}
	}
	return res;
}

MatrixXd Processer::Cov(MatrixXd matrix,int degree)
{
	MatrixXd downSamMatrix;
	if (degree == 100)
	{
		downSamMatrix=matrix;
	}
	else
	{
		size_t numofpoint=(size_t)((double)degree/100*matrix.rows());
		downSamMatrix=MatrixXd::Zero(numofpoint,matrix.cols());
		srand((unsigned)time(NULL));
		for (size_t index = 0; index < numofpoint; index++)
		{
			downSamMatrix.row(index) += matrix.row(rand() % matrix.rows());
		}
	}
	MatrixXd tmp = ZeroMean(downSamMatrix);
	MatrixXd cov = tmp.transpose()*tmp;
	cov /= (tmp.rows() - 1);
	return cov;
}

MatrixXd Processer::ZeroMean(MatrixXd matrix)
{
	MatrixXi mean = matrix.colwise().mean().cast<int>();
	MatrixXd res = MatrixXd::Zero(matrix.rows(), matrix.cols());
	for (auto col = 0; col < matrix.cols(); col++)
	{
		res.col(col) = matrix.col(col).array() - mean(0, col);
	}
	return res;
}

MatrixXd Processer::Cor(MatrixXd cov)
{
	MatrixXd Cor = MatrixXd::Zero(cov.rows(), cov.cols());
	for (auto row = 0; row < cov.rows(); row++)
	{
		for (auto col = 0; col < cov.cols(); col++)
		{
			Cor(row, col) = cov(row, col) / (sqrt(cov(row, row) + 0.0001f)*sqrt(cov(col, col) + 0.0001f));
		}
	}
	for (auto index = 0; index < cov.cols(); index++)
	{
		Cor(index, index) = 1;
	}
	return Cor;
}

void Processer::GetPLUSD(MatrixXd A, MatrixXd &P, MatrixXd &L, MatrixXd &U, MatrixXd &S, MatrixXd &D)
{
	assert(A.rows() == A.cols());

	int n = A.rows();
	L.setIdentity();
	P.setIdentity();
	S.setIdentity();
	for (auto i = 0; i < (n - 1); i++)
	{
		MatrixXd tempP = MatrixXd::Identity(n, n);
		MatrixXd tempL = MatrixXd::Identity(n, n);
		MatrixXd tempS = MatrixXd::Identity(n, n);

		int Rowmin = i;
		int Colmin = i + 1;
		double tempMin = numeric_limits<double>::max();
		for (auto j = i; j < n; j++)
		{
			for (auto k = i + 1; k < n; k++)
			{
				if(fabs(A(j, k)) < 0.000001f)
				{
					continue;
				}
				else
				{
					double Sc = fabs((A(j, i) - 1) / A(j, k));
					if (Sc < tempMin)
					{
						tempMin = Sc;
						Rowmin = j;
						Colmin = k;
					}
				}
			}
		}
		if (Rowmin != i)
		{
			tempP(i, i) = tempP(Rowmin, Rowmin) = 0;
			tempP(i, Rowmin) = tempP(Rowmin, i) = 1;
		}

		MatrixXd B = tempP*A;

		double temp_s = (B(i, i) - 1) / B(i, Colmin);
		for (int j = i + 1; j < n; j++)
		{
			tempL(j, i) = temp_s * B(j, Colmin) - B(j, i);
		}

		tempS(Colmin, i) = -temp_s;
		L = tempL * tempP * L;
		P = tempP * P;
		S *= tempS;

		A = tempL*B*tempS;
		for (auto j=i+1; j < n; j++)
		{
			A(j, i) = 0;
		}
	}
	L = L * P.transpose();
	D.setIdentity();
	D(n - 1, n - 1) = A.determinant() > 0 ? 1 : -1;
	A(n - 1, n - 1) = 1;
	U.setZero();
	U += A;
	return;
}

MatrixXd Processer::invProductLower(MatrixXd y, MatrixXd Lower)
{
	MatrixXd res = MatrixXd::Zero(y.rows(), y.cols());
	res.col(y.cols() - 1) += y.col(y.cols() - 1) / Lower(y.cols() - 1, y.cols() - 1);
	for (auto col = y.cols() - 2; col >= 0; col--)
	{
		res.col(col) += MatrixDRound(y.col(col) - res.rightCols(y.cols() - 1 - col)*Lower.col(col).bottomRows(y.cols() - 1 - col));
	}
	return res;
}

MatrixXd Processer::invProductUpper(MatrixXd y, MatrixXd Upper)
{
	MatrixXd res = MatrixXd::Zero(y.rows(), y.cols());
	res.col(0) += y.col(0) / Upper(0, 0);
	for (auto col = 1; col < y.cols(); col++)
	{
		res.col(col) += MatrixDRound(y.col(col) - res.leftCols(col)*Upper.col(col).topRows(col));
	}
	return res;
}

MatrixXd Processer::RKLT(MatrixXd matrix, MatrixXd &cov, int degree) {
	if (matrix.cols() == 1)
	{
		cov = MatrixXd::Zero(1, 1);
		return matrix;
	}
	else
	{
		cov = Cov(matrix, degree);
		MatrixXd *PLUSD = new MatrixXd[5];
		for (auto index = 0; index < 5; index++)
		{
			PLUSD[index] = MatrixXd::Zero(cov.rows(), cov.cols());
		}
#if USE_GSL
		int n = cov.rows();
		gsl_matrix_view m = gsl_matrix_view_array(cov.data(), n, n);
		gsl_vector *eval = gsl_vector_alloc(n);
		gsl_matrix *evec = gsl_matrix_alloc(n, n);
		gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(n);
		gsl_eigen_symmv(&m.matrix, eval, evec, w);
		gsl_eigen_symmv_free(w);
		gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

		MatrixXd eigenvectors(n, n);
		for (auto i =0; i < n; i++)
		{
			for (auto j=0; j < n; j++)
			{
				eigenvectors(i ,j) = gsl_matrix_get(evec, i, j);
			}
		}
		GetPLUSD(eigenvectors, PLUSD[0], PLUSD[1], PLUSD[2], PLUSD[3], PLUSD[4]);/*U=PLAS*/
		gsl_vector_free(eval);
		gsl_matrix_free(evec);
#else
		SelfAdjointEigenSolver<MatrixXd> eigenSolver(cov);
		GetPLUSD(eigenSolver.eigenvectors(), PLUSD[0], PLUSD[1], PLUSD[2], PLUSD[3], PLUSD[4]);/*U=PLAS*/
#endif
		MatrixXd invL = PLUSD[1].inverse();
		MatrixXd invP = PLUSD[0].inverse();
		MatrixXd invS = PLUSD[3].inverse();
		for (auto tmp_row = 0; tmp_row < cov.rows(); tmp_row++)
		{
			for (auto tmp_col = tmp_row + 1; tmp_col < cov.cols(); tmp_col++)
			{
				invL(tmp_row, tmp_col) = 0;
				invS(tmp_row, tmp_col) = 0;
			}
		}

		MatrixXd KLTimage = matrix*invP;
		KLTimage *= invL;
		KLTimage = MatrixDRound(KLTimage);
		KLTimage *= PLUSD[4]; // D
		KLTimage = MatrixDRound(KLTimage);
		KLTimage *= PLUSD[2]; // U
		KLTimage = MatrixDRound(KLTimage);
		KLTimage *= invS;
		KLTimage = MatrixDRound(KLTimage);
		return KLTimage;
	}
}

MatrixXd Processer::invRKLT(MatrixXd matrix, MatrixXd cov)
{
	if (matrix.cols() == 1)
	{
		return matrix;
	}
	else
	{
		MatrixXd *PLUSD = new MatrixXd[5];
		for (auto index = 0; index < 5; index++)
		{
			PLUSD[index] = MatrixXd::Zero(cov.rows(), cov.cols());
		}
#if USE_GSL
		int n = cov.rows();
		gsl_matrix_view m = gsl_matrix_view_array(cov.data(), n, n);
		gsl_vector *eval = gsl_vector_alloc(n);
		gsl_matrix *evec = gsl_matrix_alloc(n, n);
		gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(n);
		gsl_eigen_symmv(&m.matrix, eval, evec, w);
		gsl_eigen_symmv_free(w);
		gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

		MatrixXd eigenvectors(n, n);
		for (auto i =0; i < n; i++)
		{
			for (auto j=0; j < n; j++)
			{
				eigenvectors(i ,j) = gsl_matrix_get(evec, i, j);
			}
		}
		GetPLUSD(eigenvectors, PLUSD[0], PLUSD[1], PLUSD[2], PLUSD[3], PLUSD[4]);/*U=PLAS*/
		gsl_vector_free(eval);
		gsl_matrix_free(evec);
#else
		SelfAdjointEigenSolver<MatrixXd> eigenSolver(cov);
		GetPLUSD(eigenSolver.eigenvectors(), PLUSD[0], PLUSD[1], PLUSD[2], PLUSD[3], PLUSD[4]);
#endif
		MatrixXd invL = PLUSD[1].inverse();
		MatrixXd invS = PLUSD[3].inverse();
		for (auto tmp_row = 0; tmp_row < cov.rows(); tmp_row++)
		{
			for (auto tmp_col = tmp_row + 1; tmp_col < cov.cols(); tmp_col++)
			{
				invL(tmp_row, tmp_col) = 0;
				invS(tmp_row, tmp_col) = 0;
			}
		}

		MatrixXd Ori = invProductLower(matrix, invS);

		Ori = invProductUpper(Ori, PLUSD[2]);
		Ori *= PLUSD[4].inverse();
		Ori = MatrixDRound(Ori);
		Ori = invProductLower(Ori, invL);
		Ori *= PLUSD[0];
		Ori = MatrixDRound(Ori);
		return Ori;
	}
}

MatrixXd Processer::GetDistanceFromCorrelation(MatrixXd Correlation, distanceMode mode)
{
	MatrixXd res = MatrixXd::Zero(Correlation.rows(), Correlation.cols());
	const float b = 0.00001f;
	for (auto row = 0; row < Correlation.rows(); row++)
	{
		for (auto col = 0; col < Correlation.cols(); col++)
		{
			res(row, col) = (1 / (fabs(Correlation(row, col)) + b)) - 1;
			res(row, col) = res(row, col) > 0 ? res(row, col) : 0;
		}
	}
	return res;
}

void Processer::CheckIsNeedClustered(MatrixXd image, VectorXi &flags)
{
	map<double,int> probability;
	for (auto i = 0; i < image.cols(); i++)
	{
		probability.clear();
		for (auto j = 0; j < image.rows(); j++)
		{
			if (probability.find(image(j, i)) != probability.end())
			{
				probability[image(j, i)]++;
			}
			else
			{
				probability[image(j, i)] =1;
			}
		}
		int max_count=0;
		for(auto it=probability.begin();it!=probability.end();it++)
		{
			if(it->second>max_count)
			{
				max_count=it->second;
			}
		}
		if((double(max_count)/image.rows())>0.1)
		{
			flags(i)=0;
		}
	}
	return;
}

int Processer::ExtraReferenceBand(VectorXd distance, VectorXi isAbleToBeReference)
{
	assert(distance.size() == isAbleToBeReference.size());
	int length = distance.size();
	DistanceReference * distanceReference = new DistanceReference[length];
	for (auto index = 0; index < length; index++)
	{
		distanceReference[index].indexorder=index;
		distanceReference[index].distance=distance(index);
		distanceReference[index].isAbleReference=(bool)isAbleToBeReference(index);
	}
	sort(distanceReference,distanceReference+length,xcomparefunction);
	for (auto index = 0; index < length; index++)
	{
		if (distanceReference[index].isAbleReference)
		{
			int actualOrder = distanceReference[index].indexorder;
			delete[] distanceReference;
			distanceReference = nullptr;
			return actualOrder;
		}
	}
	return -1;
}

MatrixXd Processer::invPolynomialPredict(MatrixXd coeff, MatrixXd residual, MatrixXd reference)
{
	int order = coeff.rows() - 1;
	MatrixXd x = MatrixXd::Zero(residual.rows(), order + 1);
	x.col(0) += MatrixXd::Ones(residual.rows(), 1);
	for (auto index = 0; index < order; index++)
	{
		x.col(index + 1) += x.col(index).cwiseProduct(reference);
	}
	return (x*coeff).cast<int>().cast<double>() + residual;
}

int Processer::FindNewOrder(VectorXi colindex, int index)
{
	for (auto indextmp = 0; indextmp < colindex.size(); indextmp++)
	{
		if (colindex(indextmp) == index)
		{
			return indextmp;
		}
	}
	return -1;
}

void Processer::PreProcess(double *data, int length)
{
	for (auto index = 0; index < length; index++)
	{
		if (data[index] > threshold)
		{
			data[index] -= offset;
		}
	}
}

void Processer::PreProcess(MatrixXd &matrix)
{
	for (auto row = 0; row<matrix.rows(); row++)
	{
		for (auto col = 0; col<matrix.cols(); col++)
		{
			if (matrix(row, col) > threshold)
			{
				matrix(row, col) -= offset;
			}
		}
	}
}

void Processer::invPreProcess(double *data, int length)
{
	for (auto index = 0; index < length; index++)
	{
		if (data[index] < 0)
		{
			data[index] += offset;
		}
	}
}

void Processer::invPreProcess(MatrixXd &matrix)
{
	for (auto row = 0; row < matrix.rows(); row++)
	{
		for (auto col = 0; col < matrix.cols(); col++)
		{
			if (matrix(row, col) < 0)
			{
				matrix(row, col) += offset;
			}
		}
	}
}
