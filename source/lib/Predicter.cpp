/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: Predicter.cpp
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#include "Predicter.h"
#include "Common.h"
#include "Statistic.h"
#include <iostream>
#include <limits>
using namespace std;

MatrixXd* Predicter::SpectralPredict(MatrixXd cur,MatrixXd reference, PredictMode mode ,double threshold){
	MatrixXd *res = nullptr;
	double error = numeric_limits<double>::max();
	if (mode == Polynomial)
	{
		res = xPolynomialPredict(cur, reference, error);
	}
	else if (mode == DPCM)
	{
		res = xDPCMPredict(cur, reference, error);
	}
	else if (mode == Auto)
	{
		double errortmp;
		MatrixXd *restmp;
		restmp = xPolynomialPredict(cur, reference, errortmp);
		if (errortmp < error)
		{
			res = restmp;
			error = errortmp;
		}
		restmp = xDPCMPredict(cur, reference, errortmp);
		if (errortmp < error) {
			res = restmp;
			error = errortmp;
		}
	}
	double original_range = (cur.maxCoeff() - cur.minCoeff());
	Statistic state;
	double original_entropy=state.cal_entropy(cur);
	double current_entropy=state.cal_entropy(res[1]);
	if (error > threshold || ((res[1].maxCoeff() - res[1].minCoeff()) > original_range) || original_entropy<current_entropy)
	{
		res[0] = MatrixXd(1, 1);
		res[0](0, 0) = 0;
		res[1] = cur;
	}
	return res;
}

MatrixXd* Predicter::xPolynomialPredict(MatrixXd cur,MatrixXd reference,double &error,int order)
{
	assert(cur.cols() == 1 && reference.cols() == 1 && cur.rows() == reference.rows());
	MatrixXd *res = nullptr;
	if (order == -1)
	{
		MatrixXd *restmp;
		error = numeric_limits<double>::max();
		double errortmp;
		for (auto ordertmp = 1; ordertmp < MAX_ORDER; ordertmp++)
		{
			restmp = xPolynomialFitting(cur, reference, errortmp, ordertmp);
			if(errortmp<error)
			{
				res = restmp;
				error = errortmp;
			}
		}
	}
	else
	{
		res = xPolynomialFitting(cur, reference, error, order);
	}
	return res;
}

MatrixXd* Predicter::xPolynomialFitting(MatrixXd cur, MatrixXd reference, double &error, int order){
	assert(cur.cols() == 1 && reference.cols() == 1 && cur.rows() == reference.rows());
	MatrixXd *res = new MatrixXd[2];
	MatrixXd x = MatrixXd::Zero(cur.rows(), order + 1);
	x.col(0) += MatrixXd::Ones(cur.rows(), 1);
	for (auto index = 0; index < order; index++)
	{
		x.col(index + 1) += x.col(index).cwiseProduct(reference);
	}

	res[0] = (x.transpose()*x).inverse()*x.transpose()*cur;
	res[1] = cur - (x*res[0]).cast<int>().cast<double>();
	error = res[1].cwiseAbs2().sum() / res[1].rows();
	return res;
}

MatrixXd* Predicter::xDPCMPredict(MatrixXd cur,MatrixXd reference,double &error){
	assert(cur.cols() == 1 && reference.cols() == 1 && cur.rows() == reference.rows());
	MatrixXd *res = new MatrixXd[2];
	res[0] = MatrixXd::Zero(2, 1);
	res[0](0, 0) = 0; res[0](1, 0) = 1;
	res[1] = cur - reference;
	error = res[1].cwiseAbs2().sum() / res[1].rows();
	return res;
}
