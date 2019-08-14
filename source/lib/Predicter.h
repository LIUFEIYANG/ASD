/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: Predicter.h
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#ifndef PREDICTER_H
#define PREDICTER_H

#include "Common.h"
#include <Eigen/Eigen>
using namespace Eigen;

enum PredictMode
{
	Polynomial,
	DPCM,
	Auto,
	AllPredictMode
};

class Predicter
{
private:
	static MatrixXd* xPolynomialPredict(MatrixXd cur, MatrixXd reference, double &error, int order = -1);
	static MatrixXd* xPolynomialFitting(MatrixXd cur, MatrixXd reference, double &error, int order);
	static MatrixXd* xDPCMPredict      (MatrixXd cur, MatrixXd reference, double &error);

public:
	static MatrixXd* SpectralPredict(MatrixXd cur,MatrixXd reference,PredictMode mode,double threshold);
};

#endif
