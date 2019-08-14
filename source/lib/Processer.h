/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: Processer.h
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#ifndef PROCESSER_H
#define PROCESSER_H

#include <Eigen/Eigen>
using namespace Eigen;

struct DistanceReference
{
	int    indexorder;
	double distance;
	bool   isAbleReference;
};

enum distanceMode
{
	distancemode1,
	distancemode2,
	distancemode3,
	distancemode4,
	AllDistanceMode
};

class Processer
{
private:
	static bool     xcomparefunction(DistanceReference a, DistanceReference b);

public:
	static MatrixXd MatrixDRound(MatrixXd matrix);                                                               // round elements of matrix
	static MatrixXd Cov(MatrixXd matrix, int degree);                                                            // calculate covariance by col major
	static MatrixXd Cor(MatrixXd cov);                                                                           // get correlations form covariances
	static MatrixXd ZeroMean(MatrixXd matrix);                                                                   // zero meaning a matrix by col major
	static void     GetPLUSD(MatrixXd A, MatrixXd &P, MatrixXd &L, MatrixXd &U, MatrixXd &S, MatrixXd &D);       // decompose A to corresponding P*L*U*D*S
	static MatrixXd invProductLower(MatrixXd y, MatrixXd Lower);                                                 // inverse product of lower triangle ERM
	static MatrixXd invProductUpper(MatrixXd y, MatrixXd Upper);                                                 // inverse product of upper triangle ERM
	static MatrixXd RKLT(MatrixXd matrix, MatrixXd &cov, int degree);                                            // reversible K-L transform
	static MatrixXd invRKLT(MatrixXd matrix, MatrixXd cov);                                                      // inverse reversible K-L transform
	static MatrixXd GetDistanceFromCorrelation(MatrixXd Correlation,distanceMode mode);                          // derive distance from correlation

	static int      ExtraReferenceBand(VectorXd distance, VectorXi isAbleToBeReference);                         // obtain the reference band
	static MatrixXd invPolynomialPredict(MatrixXd coeff, MatrixXd residual, MatrixXd reference);                 // inverse polynomial prediction
	static int      FindNewOrder(VectorXi colindex, int index);                                                  // search the new order of a band
	static void     PreProcess(double *data, int length);                                                        // preprocess
	static void     PreProcess(MatrixXd &matrix);                                                                // preprocess
	static void     invPreProcess(double *data, int length);                                                     // inverse preprocess
	static void     invPreProcess(MatrixXd &matrix);                                                             // inverse preprocess
	static void     CheckIsNeedClustered(MatrixXd image, VectorXi &flag);                                        // check the transform is useful for each band
};

#endif
