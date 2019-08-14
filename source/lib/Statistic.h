/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: Statistic.h
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#ifndef STATISTIC_H
#define STATISTIC_H

#include "Eigen/Eigen"
using namespace Eigen;

class Statistic{
public:
	static double cal_entropy(MatrixXd image);
	static double cal_spatial_correlation(MatrixXd image,int offset_x=1,int offset_y=1);
};

#endif
