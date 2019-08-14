/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: Statistic.cpp
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#include "Statistic.h"
#include "Processer.h"
#include <map>
#include <cmath>
using namespace std;

double Statistic::cal_entropy(MatrixXd image) {
	map<double, int> counter;
	for (auto i = 0; i < image.rows(); i++)
	{
		for (auto j = 0; j < image.cols(); j++)
		{
			if(counter.find(image(i,j)) == counter.end())
			{
				counter[image(i, j)] = 1;
			}
			else
			{
				counter[image(i, j)] = counter[image(i, j)] + 1;
			}
		}
	}

	double entropy = 0;
	for (auto it = counter.begin(); it != counter.end(); it++)
	{
		double probability = (double)it->second / (image.rows()*image.cols());
		entropy += probability * log(1 / probability) / log(2);
	}
	return entropy;
}

double Statistic::cal_spatial_correlation(MatrixXd image, int offset_x, int offset_y) {
	int count = (image.rows() - offset_y)*(image.cols() - offset_x);
	MatrixXd seq(count, 2);
	int iter = 0;
	for (auto i = 0; i < (image.rows() - offset_y); i++)
	{
		for (auto j = 0; j < (image.cols() - offset_x); j++)
		{
			seq(iter, 0) = image(i, j);
			seq(iter, 1) = image(i + offset_y, j + offset_x);
			iter++;
		}
	}
	MatrixXd cov = Processer::Cov(seq, 100);
	return cov(1, 0) / (sqrt(cov(0, 0))*sqrt(cov(1, 1)));
}
