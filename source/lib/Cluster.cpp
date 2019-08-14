/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: Cluster.cpp
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#include "Cluster.h"
#include "Common.h"
#include <iostream>
#include <math.h>
#include <cassert>
#include <algorithm>
#include <vector>
using namespace std;

bool Cluster::xcomparefunction(IndexAndNum a, IndexAndNum b)
{
	return a.num > b.num;
}

bool Cluster::xcomparefunction2(IndexAndNum a, IndexAndNum b)
{
	return a.gamma > b.gamma;
}

bool Cluster::xcomparevectors(vector<int> a, vector<int> b)
{
	if (a.size() != b.size())
	{
		return false;
	}

	if (a.size() == 0)
	{
		return false;
	}

	for (auto i = 0; i<a.size(); i++)
	{
		if (a[i] != b[i])
		{
			return false;
		}
	}
	return true;
}

VectorXd Cluster::xGetRho(MatrixXd distance, double threshold)
{
	VectorXd res = VectorXd::Zero(distance.cols());
	for (auto i = 0; i<distance.rows(); i++)
	{
		for (auto j = 0; j < distance.cols(); j++)
		{
			if( i != j )
			{
				res(i) += exp((-1)*(distance(i, j)*distance(i, j)) / (threshold*threshold));
			}
		}
	}
	return res;
}

double Cluster::xGetEntropy(VectorXd rho)
{
	double entropy = 0;
	double sum = 0;
	for (auto i = 0; i < rho.size(); i++)
	{
		sum += rho(i);
	}
	for (auto i = 0; i < rho.size(); i++)
	{
		double ratio = rho(i) / sum;
		entropy -= ratio * log(ratio); // / log(2)
	}
	return entropy;
}

double Cluster::xGetDc(MatrixXd distance, double precision)
{
	// ternary search
	double l = distance.maxCoeff();
	double r = distance.minCoeff();
	for (auto row = 0; row<(distance.rows() - 1); row++)
	{
		for (auto col = row + 1; col<distance.cols(); col++)
		{
			if (distance(row, col) < l)
			{
				l = distance(row, col);
			}
			if (distance(row, col) > r)
			{
				r = distance(row, col);
			}
		}
	}

	//double l = 0, r = distance.maxCoeff();

	while ((r - l) > precision) {
		double mid1 = (2 * l + r) / 2;
		double mid2 = (l + 2 * r) / 3;
		VectorXd rho1   = xGetRho(distance, mid1);
		double entropy1 = xGetEntropy(rho1);
		VectorXd rho2   = xGetRho(distance, mid2);
		double entropy2 = xGetEntropy(rho2);
		if (entropy1 < entropy2)
		{
			r = mid2;
		}
		else
		{
			l = mid1;
		}
	}
	return 3 * l / sqrt(2);
}

MatrixXd Cluster::xGetDelta(MatrixXd distance, VectorXd rho)
{
	MatrixXd res = MatrixXd::Constant(2, distance.cols(),-1);
	res.row(0) = distance.colwise().maxCoeff();
	for (auto i = 0; i < distance.cols(); i++)
	{
		for (auto j = 0; j < distance.cols(); j++)
		{
			double curDis = distance(i, j);
			if ((i != j) && (rho(i) <= rho(j)) && (res(0, i) >= curDis))
			{
				if (res(1, j) != i)
				{
					res(0, i) = curDis;
					res(1, i) = j;
				}
			}
		}
	}
	return res;
}

VectorXd Cluster::xGetGamma(VectorXd delta, VectorXd rho)
{
	double minDelta = delta.minCoeff();
	double maxDelta = delta.maxCoeff();
	double minRho = rho.minCoeff();
	double maxRho = rho.maxCoeff();
	VectorXd normDelta = (delta) / (maxDelta);
	VectorXd normRho = (rho) / (maxRho);
	for (auto i = 0; i < normRho.size(); i++)
	{
		normRho(i) = pow(normRho(i), 0.05);
	}
	return normDelta.cwiseProduct(normRho);
}

MatrixXi Cluster::xGetClusterCenter(VectorXd gamma, MatrixXd cor, double EPS, int numofCluster, MatrixXd deltaandpos)
{
	MatrixXd gammaPair = MatrixXd::Zero(2, gamma.size());
	for (auto i = 0; i<gamma.size(); i++)
	{
		gammaPair(0, i) = gamma(i);
		gammaPair(1, i) = i;
	}
	MatrixXd gammaSort = MatrixXd::Zero(2, gamma.size());
	for (auto i = 0; i<gamma.size(); i++)
	{
		int position;
		gamma.maxCoeff(&position);
		gammaSort.col(i) += gammaPair.col(position);
		gamma(position) = -1;
	}
	numofCluster = L_min(gamma.size(), numofCluster);
	if (numofCluster <= 0)
	{
		int length;
		double threshold = 0;
		double factor=1;
		for (length = gamma.size(); (length>1) && (gammaSort(0, length - 1)<EPS); length--);

		// geometry mean value as threshold
		//for (auto i = 0; i < length; i++) {
		//  threshold += log(gammaSort(0, i));
		//}
		//threshold /= length;
		//threshold = factor*exp(threshold);

		// unequal weight mean as threshold
		int tot;
		for (tot = 0; length > 1; tot++, length >>= 1)
		{
			double cur = 0;
			for (int j = 0; j < length; j++)
			{
				cur += log(gammaSort(0, j));
			}
			threshold += cur / length;
		}
		threshold = factor*exp(threshold / tot);

		// cout << "threshold : " << threshold << endl;
		numofCluster = 0;
		for (auto i = 0; (i < gamma.size()) && (gammaSort(0, i) >= threshold); i++)
		{
			numofCluster++;
		}
		// block : check whether the minimal coefficient > 0.8 or maximal coefficient < 0.8
		{
			VectorXd cortmp(cor.rows()*(cor.cols() - 1) / 2);
			int tmpindex = 0;
			for (auto i = 0; i < cor.rows(); i++)
			{
				for (auto j = i + 1; j < cor.cols(); j++)
				{
					cortmp(tmpindex) = fabs(cor(i, j));
					tmpindex++;
				}
			}
			if (cortmp.minCoeff() >= 0.8)
			{
				numofCluster = 1;
			}
			else if (cortmp.maxCoeff() < 0.8)
			{
				numofCluster = gamma.size();
			}
		}
	}

	MatrixXi res = MatrixXi::Constant(2, gamma.size(), -1);
	res.row(1).setZero();
	for (auto i = 0; i < numofCluster; i++)
	{
		res(0, gammaSort(1, i)) = (int)gammaSort(1, i);
		res(1, gammaSort(1, i)) = 1;
	}
	return res;
}

void Cluster::xGetCluster(MatrixXi &ClusterAndFlag, MatrixXd deltaAndPos)
{
	UnionSet unionset(ClusterAndFlag.cols());
	for (auto i = 0; i < ClusterAndFlag.cols(); i++)
	{
		if (ClusterAndFlag(1, i) != 1)
		{
			unionset.merge(i, (int)deltaAndPos(1, i));
		}
	}
	for (auto i = 0; i < ClusterAndFlag.cols(); i++)
	{
		if (ClusterAndFlag(1, i) == 1)
		{
			ClusterAndFlag(0, unionset.getFather(i)) = i;
		}
	}
	for (auto i = 0; i < ClusterAndFlag.cols(); i++)
	{
		ClusterAndFlag(0, i) = ClusterAndFlag(0, unionset.getFather(i));
	}
}

VectorXi Cluster::xGetClusterCore(VectorXi Cluster, MatrixXd distance, VectorXd rho, double dc)
{
	VectorXd maxRhoBorderRegion = VectorXd::Zero(distance.cols());
	for (auto i = 0; i < distance.cols(); i++)
	{
		if (Cluster(i) == -1)
		{
			continue;
		}
		bool flag = false;
		for (auto j = 0; j < distance.cols(); j++)
		{
			flag |= (Cluster(i) != Cluster(j)) && (distance(i, j) <= dc);
		}
		if (flag)
		{
			maxRhoBorderRegion(Cluster(i)) = L_max(maxRhoBorderRegion(Cluster(i)), rho(i));
		}
	}
	VectorXi clusterCore = VectorXi::Zero(distance.cols());
	for (auto i = 0; i < distance.cols(); i++)
	{
		if (Cluster(i) != -1)
		{
			clusterCore(i) = (rho(i) > maxRhoBorderRegion(Cluster(i))) ? 1 : 0;
		}
	}
	return clusterCore;
}

void Cluster::xChangeCluster(VectorXi &Cluster)
{
	VectorXi newNum = VectorXi::Constant(Cluster.size(),-1);
	int tot = 0;
	for (auto i = 0; i < Cluster.size(); i++)
	{
		if (Cluster(i) == i)
		{
			newNum(i) = tot++;
		}
	}
	for (auto i = 0; i < Cluster.size(); i++)
	{
		if (Cluster(i) != -1)
		{
			Cluster(i) = newNum(Cluster(i));
		}
	}
}

VectorXi Cluster::DensityPeakCluster(MatrixXd source, MatrixXd distance, MatrixXd cor, VectorXd &gammavalue, int numofCluster)
{
	double   EPS = 0; // = 0; // if set as 0, all the gamma values are used to calculate the threshold
	// cout << "distance : " << endl << distance << endl;
	double   dc_value = xGetDc(distance, 1e-4);
	// cout << "range of distance : " << distance.minCoeff() << " - " << distance.maxCoeff() << endl;
	// cout << "dc value : " << dc_value << endl;
	VectorXd rho = xGetRho(distance,dc_value);
	// cout << "density : " << (RowVectorXd)rho << endl;
	MatrixXd delta = xGetDelta(distance, rho);
	// cout << "delta distance and position  : " << endl << delta << endl;
	VectorXd gamma = xGetGamma((VectorXd)delta.row(0), rho);
	gammavalue = gamma;
	MatrixXi clusterAndFlag = xGetClusterCenter(gamma, cor, EPS, numofCluster, delta);
	//cout << clusterAndFlag << endl;
	xGetCluster(clusterAndFlag, delta);
	//cout << clusterAndFlag << endl;
	VectorXi cluster = clusterAndFlag.row(0);
	//if (numofCluster <= 0) {
	//  VectorXi clusterCore = xGetClusterCore(cluster, distance, rho, dc_value);
	//  for (auto i = 0; i < clusterCore.size(); i++) {
	//    if (clusterCore(i))
	//      cluster(i) = i;
	//  }
	//}

	//map<double, int> probability;
	//for (auto i = 0; i < source.cols(); i++) {
	//  probability.clear();
	//  for (auto j = 0; j < source.rows(); j++) {
	//    if (probability.find(source(j, i)) != probability.end()) {
	//      probability[source(j, i)]++;
	//    }
	//    else {
	//      probability[source(j, i)] = 1;
	//    }
	//  }
	//  int max_count = 0;
	//  for (auto it = probability.begin(); it != probability.end(); it++) {
	//    if (it->second>max_count) {
	//      max_count = it->second;
	//    }
	//  }
	//  if ((double(max_count) / source.rows())>0.5) {
	//    cluster(i) = i;
	//  }
	//}

updateclusterslabel:
	// cout << "clusters : " << endl << (RowVectorXi)cluster << endl;
	map<int, int> clusterstate;
	clusterstate.clear();
	for (auto i = 0; i < cluster.size(); i++)
	{
		if (clusterstate.find(cluster(i)) != clusterstate.end())
		{
			clusterstate[cluster(i)]++;
		}
		else
		{
			clusterstate[cluster(i)] = 1;
		}
	}
	vector<int> *tmpstate=new vector<int>[clusterstate.size()];
	int index = 0;
	for (auto i = clusterstate.begin(); i != clusterstate.end(); i++, index++)
	{
		for (auto j = 0; j < cluster.size(); j++)
		{
			if (cluster(j) == i->first)
			{
				tmpstate[index].push_back(j);
			}
		}
	}
	index = 0;
	for (auto i = clusterstate.begin(); i != clusterstate.end(); i++, index++)
	{
		if (i->second != 1)
		{
			// cout << " index : " << i->first << endl;
			MatrixXd cortmp = MatrixXd::Ones(i->second, i->second);
			for (auto row = 0; row < i->second; row++)
			{
				for (auto col = 0; col < i->second; col++)
				{
					if (row != col)
					{
						cortmp(row, col) = fabs(cor(tmpstate[index][row], tmpstate[index][col]));
					}
				}
			}
			// cout << cortmp << endl;
			MatrixXd weight = cortmp.colwise().sum();
			for (auto tmpj = 0; tmpj < i->second; tmpj++)
			{
				if (tmpstate[index][tmpj]==i->first)
				{
					weight(tmpj) = i->second + 1;
				}
			}
			// cout << weight << endl;
			int minimalcoeffrow,minimalcoeffcol;
			weight.minCoeff(&minimalcoeffrow, &minimalcoeffcol);
			double minimalcorrelation = cortmp.col(minimalcoeffcol).minCoeff();
			// cout << "minimal coefficient : " << minimalcorrelation << endl;
			if (minimalcorrelation < 0.75)
			{
				cluster(tmpstate[index][minimalcoeffcol]) = tmpstate[index][minimalcoeffcol];
				delete[] tmpstate;
				tmpstate = nullptr;
				goto updateclusterslabel;
			}
		}
	}
	xChangeCluster(cluster);
	return cluster;
}

VectorXi Cluster::GetClusterGroupState(VectorXi Cluster)
{
	VectorXi res(Cluster.maxCoeff() + 1);
	res.setZero();
	for (int index = 0; index < Cluster.size(); index++)
	{
		res(Cluster(index))++;
	}
	return res;
}

VectorXi Cluster::ReOrderClusterByNum(VectorXi Cluster, VectorXd gamma)
{
	VectorXi groupState = GetClusterGroupState(Cluster);
	int length = groupState.size();
	IndexAndNum * indexandnum=new IndexAndNum[length];
	int tmpcount=0;
	for (auto index = 0; index < length; index++)
	{
		indexandnum[index].index=index;
		indexandnum[index].num=groupState(index);
		for (auto indextmp = 0; indextmp < Cluster.size(); indextmp++)
		{
			if (Cluster(indextmp) == index)
			{
				indexandnum[index].gamma = gamma(indextmp);
				break;
			}
		}
		if(groupState(index)==1)
		{
			tmpcount++;
		}
	}
	sort(indexandnum,indexandnum+length,xcomparefunction);
	sort(indexandnum+length-tmpcount,indexandnum+length,xcomparefunction2);
	VectorXi res(length);
	for (auto index = 0; index < length; index++)
	{
		res(index) = indexandnum[index].index;
	}
	return res;
}

void Cluster::UpdateCluster(VectorXi &Cluster, VectorXi isNeedTransform)
{
	for (auto index = 0; index < isNeedTransform.size(); index++)
	{
		if (!isNeedTransform(index))
		{
			Cluster(index)=Cluster.maxCoeff();
		}
	}
}

UnionSet::UnionSet(int size)
{
	father = new int[size];
	for (auto index = 0; index<size; index++)
	{
		father[index] = index;
	}
}

UnionSet::~UnionSet()
{
	if (father)
	{
		delete[] father;
	}
	father = nullptr;
}

int UnionSet::getFather(int x)
{
	if (x == father[x])
	{
		return x;
	}
	return father[x] = getFather(father[x]);
}

void UnionSet::merge(int x, int y) {
	int fatherx = getFather(x);
	int fathery = getFather(y);
	if (fatherx == fathery)
	{
		return;
	}
	father[fatherx] = fathery;
	return;
}
