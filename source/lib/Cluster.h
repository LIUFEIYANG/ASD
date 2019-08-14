/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: Cluster.h
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <Eigen/Eigen>
using namespace Eigen;
using namespace std;

struct IndexAndNum
{
	int    index;
	int    num;
	double gamma;
};

enum clusterMode
{
	DensityPeakCluster,
	AffinityPropagationCluster,
	AllClusterMode
};

class Cluster
{
private:
	static VectorXd xGetRho(MatrixXd distance, double threshold);
	static double   xGetEntropy(VectorXd rho);
	static double   xGetDc(MatrixXd distance, double precision);
	static MatrixXd xGetDelta(MatrixXd distance, VectorXd rho);
	static VectorXd xGetGamma(VectorXd delta, VectorXd rho);
	static void     xGetCluster(MatrixXi &ClusterAndFlag, MatrixXd deltaAndPos);
	static VectorXi xGetClusterCore(VectorXi Cluster, MatrixXd distance, VectorXd rho, double dc);
	static void     xChangeCluster(VectorXi &Cluster);

	static bool     xcomparefunction(IndexAndNum a, IndexAndNum b);
	static bool     xcomparefunction2(IndexAndNum a, IndexAndNum b);
	static bool     xcomparevectors(vector<int> a, vector<int> b);

public:
	static VectorXi DensityPeakCluster(MatrixXd source, MatrixXd distance, MatrixXd cor, VectorXd &gammavalue, int numofCluster = 0);

	static MatrixXi xGetClusterCenter(VectorXd gamma, MatrixXd cor, double EPS, int numofCluster, MatrixXd deltaandpos);
	static VectorXi GetClusterGroupState(VectorXi Cluster);
	static VectorXi ReOrderClusterByNum(VectorXi Cluster,VectorXd gamma);
	static void     UpdateCluster(VectorXi &Cluster, VectorXi isNeedTransform);
};


class UnionSet
{
private:
	int *father;

public:
	UnionSet(int size);
	~UnionSet();

	int getFather(int x);
	void merge(int x, int y);
};

#endif
