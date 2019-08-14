#include "../lib/SyntaxWriter.h"
#include "../lib/Processer.h"
#include "../lib/Cluster.h"
#include "../lib/ImageReader.h"
#include "../lib/ImageWriter.h"
#include "../lib/Predicter.h"
#include "../lib/Common.h"
#include "../lib/cmdline.h"
#include <cassert>
#include <fstream>

void setup_parameter(cmdline::parser &parameters)
{
	// mandatory : required
	parameters.add<string>("input"          , 'i', "input raw image file"                               ,  true, "");
	parameters.add<string>("output"         , 'o', "output transformed file"                            ,  true, "");
	parameters.add<string>("auxi"           , 'a', "auxilary data to recover original image"            ,  true, "");
	parameters.add<string>("inputsize"      , 's', "size of input raw image file"                       ,  true, "");
	parameters.add<string>("outputsize"     , 'f', "size of output raw image file"                      ,  true, "");
	// optional : with default value
	parameters.add<int>   ("degree"         , 'd', "percentage of points to calculate covariance matrix", false, 100, cmdline::range(1,100));
	parameters.add<int>   ("numberofcluster", 'c', "number of cluster"                                  , false,   0, cmdline::range(0,1000));
		// 0 means adaptive, positive value means corresponding clusters
	parameters.add<int>   ("distance"       , 'r', "functions to calculate distance from correlation"   , false,   1, cmdline::range(0, (int)distanceMode::AllDistanceMode - 1));
	parameters.add<int>   ("cluster"        , 'w', "clustering algorithm decision"                      , false,   0, cmdline::range(0, (int)clusterMode::AllClusterMode - 1));
	// switch
	parameters.add        ("preprocess"     , '\0', "using prepocess");
	parameters.add        ("checklossless"  , '\0', "check the losslessness of transformation");
	parameters.add        ("useplse"        , '\0', "whether use PLSE");
	parameters.add        ("userklt"        , '\0', "whether use RKLT");
	parameters.add        ("verbose"        , '\0', "whether display information");
	return;
}

image_geometry parse_size(string size)
{
	/* geometry format should be
		width:height:band:depth:sign:endianess
		sign : 1 means signed, 0 means unsigned
		endianess : 1 means big endian, 0 means little endian
	*/
	image_geometry geo;
	const char * delim = ":"; // delimter
	vector<string> tmp;
	char *p = nullptr;
	p = strtok(const_cast<char*>(size.c_str()),delim);
	while(p)
	{
		tmp.push_back(string(p));
		p = strtok(NULL,delim);
	}
	assert(tmp.size()==6);
	geo.width     = atoi(tmp[0].c_str());
	geo.height    = atoi(tmp[1].c_str());
	geo.band      = atoi(tmp[2].c_str());
	geo.depth     = atoi(tmp[3].c_str());
	geo.sign      = (bool)atoi(tmp[4].c_str());
	geo.endianess = (bool)atoi(tmp[5].c_str());
	return geo;
}

int main(int argc, char* argv[])
{
	// print the input command line
	printCmdLine(argc, argv);

	// parse command line
	cmdline::parser parameters;
	setup_parameter(parameters);
	parameters.parse_check(argc,argv);

	// setup parameters from command line parser
	string infn              = parameters.get<string>("input");
	string oufn              = parameters.get<string>("output");
	string aufn              = parameters.get<string>("auxi");
	int    degree            = parameters.get<int>("degree");
	int    countofclusters   = parameters.get<int>("numberofcluster");
	clusterMode  clsMode     = clusterMode(parameters.get<int>("cluster"));
	distanceMode disMode     = distanceMode(parameters.get<int>("distance"));
	image_geometry igeometry = parse_size(parameters.get<string>("inputsize"));
	image_geometry ogeometry = parse_size(parameters.get<string>("outputsize"));
	bool isPreprocess        = parameters.exist("preprocess") ? true : false;
	bool isCheckLossless     = parameters.exist("checklossless") ? true : false;
	bool isCheckTransformed  = true;
	bool isUsePLSE           = parameters.exist("useplse") ? true : false;
	bool isUseRKLT           = parameters.exist("userklt") ? true : false;
	bool isVerbose           = parameters.exist("verbose") ? true : false;

	// force checking the input and output geometry parameters
	assert( igeometry.width  == ogeometry.width &&
			igeometry.height == ogeometry.height &&
			igeometry.band   == ogeometry.band);

	// print parsed parameters
	if (isVerbose)
	{
		cout << "inputfilename      : " << infn << endl;
		cout << "outputfilename     : " << oufn << endl;
		cout << "auxilaryfilename   : " << aufn << endl;
		cout << "input geometry       " << endl;
		cout << "-size              : " << igeometry.height << "*" << igeometry.width << "*" << igeometry.band << endl;
		cout << "-depth             : " << igeometry.depth << endl;
		cout << "-sign              : " << igeometry.sign << endl;
		cout << "-endianess         : " << (igeometry.endianess ? "big endian" : "little endian") << endl;
		cout << "output geometry      " << endl;
		cout << "-size              : " << ogeometry.height << "*" << ogeometry.width << "*" << ogeometry.band << endl;
		cout << "-depth             : " << ogeometry.depth << endl;
		cout << "-sign              : " << ogeometry.sign << endl;
		cout << "-endianess         : " << (ogeometry.endianess ? "big endian" : "little endian") << endl;
		cout << "transform parameters " << endl;
		cout << "-simplifydegree    : " << degree << endl;
		cout << "-distancefunction  : " << disMode << endl;
		cout << "-preprocess        : " << (isPreprocess ? "yes" : "no") << endl;
		cout << "-checklossless     : " << (isCheckLossless ? "yes" : "no") << endl;
		cout << "-use PLSE          : " << (isUsePLSE ? "yes" : "no") << endl;
		cout << "-use RKLT          : " << (isUseRKLT ? "yes" : "no") << endl;
		cout << "-countofclusters   : " << countofclusters << endl << endl;
	}


	if (isUsePLSE && !isUseRKLT)
	{
		cout << "non-sense spectral decorrelation configuration!" << endl;
		return ERROR_NONSENSE_CFG;
	}

	int lengthPerBand = igeometry.width*igeometry.height; // w * h
	int lengthTotal = igeometry.band*lengthPerBand;       // w * h * b

	Cluster cls;
	Predicter pec;

	// reading image data
	ImageReader imageReader(infn);
	double *dataDouble = imageReader.ReadRAWdata(igeometry);

	if (igeometry.sign)
	{
		isPreprocess = false;
	}
	if (isPreprocess)
	{
		// Processer::PreProcess(dataDouble, lengthTotal);
		cout << "pre-process should not be enabled now!" << endl;
		exit(0);
	}

	Map<MatrixXd> MatrixDoubleOriNoDemean(dataDouble, lengthPerBand, igeometry.band); // original data
	if (isVerbose)
	{
		cout << "data range         : " << MatrixDoubleOriNoDemean.minCoeff() << " - " << MatrixDoubleOriNoDemean.maxCoeff() << endl;
	}

	MatrixXi meanOri = MatrixDoubleOriNoDemean.colwise().mean().cast<int>();  // colwise mean vector
	MatrixXd MatrixDoubleOri = Processer::ZeroMean(MatrixDoubleOriNoDemean);  // de-meaned data

	VectorXi OriColIndex       = VectorXi::Zero(igeometry.band);
	VectorXi isNeedClustered = VectorXi::Ones(igeometry.band);

	Processer::CheckIsNeedClustered(MatrixDoubleOri, isNeedClustered);
	if(isVerbose)
	{
		cout << "bands to be clustered : " << (RowVectorXi)isNeedClustered << endl;
	}


	MatrixXd MatrixDoubleOriReOrder = MatrixXd::Zero(MatrixDoubleOri.rows(), MatrixDoubleOri.cols());

	/* block 1 : re-order the image
		bands should be clustered are put forward while bands should not be clustered are put backward
	*/
	{
		int p_f = 0;
		int p_b = isNeedClustered.sum();
		for (auto i = 0; i < isNeedClustered.size(); i++)
		{
			if(isNeedClustered(i) == 1)
			{
				MatrixDoubleOriReOrder.col(p_f) += MatrixDoubleOri.col(i);
				OriColIndex(p_f) = i;
				p_f++;
			}
			else
			{
				MatrixDoubleOriReOrder.col(p_b) += MatrixDoubleOri.col(i);
				OriColIndex(p_b) = i;
				p_b++;
			}
		}
	}

	VectorXd gamma   = VectorXd::Zero(igeometry.band);
	VectorXi cluster = VectorXi::Zero(igeometry.band);

	/* block 2 : cluster the bands need to be clustered
	*/
	{
		MatrixXd cov = Processer::Cov(MatrixDoubleOriReOrder.leftCols(isNeedClustered.sum()), degree);
		MatrixXd cor = Processer::Cor(cov);
		MatrixXd distance = Processer::GetDistanceFromCorrelation(cor, disMode);
		VectorXd gammatmp;
		VectorXi clustertmp;
		if (clsMode == DensityPeakCluster)
		{
			clustertmp = cls.DensityPeakCluster(MatrixDoubleOriReOrder.leftCols(isNeedClustered.sum()), distance, cor, gammatmp, countofclusters);
		}
		cluster.head(isNeedClustered.sum()) += clustertmp;
		gamma.head(isNeedClustered.sum()) += gammatmp;

		for (auto i = isNeedClustered.sum(); i < isNeedClustered.size(); i++)
		{
			cluster(i) = cluster.maxCoeff() + 1;
		}
	}

	if (isVerbose){
		cout << "original column index : " << (RowVectorXi)OriColIndex << endl;
		cout << "clusters : " << (RowVectorXi)cluster << endl;
	}

	/* block 3 : adjust the cluster and gamma to its original index
	*/
	{
		VectorXi clustertmp(cluster.size());
		VectorXd gammatmp(gamma.size());
		for (auto i = 0; i < OriColIndex.size(); i++)
		{
			clustertmp(OriColIndex(i)) = cluster(i);
			gammatmp(OriColIndex(i)) = gamma(i);
			OriColIndex(i) = i;
		}
		cluster.setZero();
		cluster += clustertmp;
		gamma.setZero();
		gamma += gammatmp;
	}

	// dealing with the image according to the clusters
	MatrixXd cov = Processer::Cov(MatrixDoubleOri, degree);
	MatrixXd cor = Processer::Cor(cov);
	MatrixXd distance = Processer::GetDistanceFromCorrelation(cor, disMode);
	VectorXi groupstate = cls.GetClusterGroupState(cluster);
	VectorXi groupReorder = cls.ReOrderClusterByNum(cluster, gamma);

	if(isVerbose)
	{
		cout << "cov : " << endl << cov << endl << endl;
		cout << "cor : " << endl << cor << endl << endl;
		cout << "cluster            : " << (RowVectorXi)cluster << endl;
		cout << "group status       : " << (RowVectorXi)groupstate << endl;
		cout << "group rank         : " << (RowVectorXi)groupReorder << endl;
		cout << "gamma value        : " << (RowVectorXd)gamma << endl;
	}

	/* ReferenceList is used to help extract the band able to be refercenced */
	VectorXi ReferenceList = VectorXi::Zero(igeometry.band);

	/* covList is used to store the cov matrices
		for bands to RKLT, it stores the cov matrix
		for bands to PLSE, it stores the coefficient vector
	*/
	MatrixXd *covList = new MatrixXd[groupstate.size()];


	/* out is the used to store the transformed data */
	MatrixXd out = MatrixXd::Zero(lengthPerBand, igeometry.band);

	int numofoutband = 0;
	for (auto index = 0; index < groupstate.size(); index++)
	{
		// use both PLSE and RKLT
		if (isUsePLSE && isUseRKLT)
		{
			MatrixXd tmpIn = MatrixXd::Zero(lengthPerBand, groupstate(groupReorder(index)));
			if (isVerbose)
			{
				cout << "==================================================================================" << endl;
				cout << "dealing with " << groupstate(groupReorder(index)) << " band : ";
			}

			// obtain the corresponding data
			int j = 0;
			for (auto i = 0; i < cluster.size(); i++)
			{
				if (cluster(i) == groupReorder(index))
				{
					if (isVerbose)
					{
						cout << i << " ";
					}
					tmpIn.col(j) += MatrixDoubleOri.col(i);
					OriColIndex(numofoutband + j) = i;
					if (groupstate(groupReorder(index)) > 1)
					{
						ReferenceList(i) = 1;
					}
					j++;
				}
			}
			if (isVerbose)
			{
				cout << endl << "Reference list " << (RowVectorXi)ReferenceList << endl;
			}

			MatrixXd tmpOut;
			if (groupstate(groupReorder(index)) == 1)
			{
				int indextmp = Processer::ExtraReferenceBand(distance.col(OriColIndex(numofoutband)), ReferenceList);
				if (indextmp != -1)
				{
					/*  tmpres[0] contain the coefficients vector
					    tmpres[1] contain the residuals
					*/
					MatrixXd *tmpres = pec.SpectralPredict(tmpIn, MatrixDoubleOri.col(indextmp), Auto, cov(OriColIndex(numofoutband), OriColIndex(numofoutband)));
					if (tmpres[0].rows() != 1 || tmpres[0](0, 0) != 0)
					{
						covList[groupReorder(index)] = MatrixXd::Zero(tmpres[0].rows() + 1, 1);
						covList[groupReorder(index)].topRows(tmpres[0].rows()) += tmpres[0];
						covList[groupReorder(index)](tmpres[0].rows(), 0) = indextmp;
						tmpOut = tmpres[1];
						if (isVerbose)
						{
							cout << "reference index : " << indextmp << endl;
							cout << "coefficience    : " << tmpres[0].transpose() << endl;
						}
					}
					else
					{
						tmpOut = tmpres[1];
						covList[groupReorder(index)] = MatrixXd::Zero(1, 1);
						covList[groupReorder(index)](0, 0) = -1;
						if (isVerbose)
						{
							cout << "reference index : " << -1 << endl;
						}
					}
				}
				else
				{
					tmpOut = tmpIn;
					covList[groupReorder(index)] = MatrixXd::Zero(1, 1);
					covList[groupReorder(index)](0, 0) = indextmp;
					if (isVerbose)
					{
						cout << "reference index : " << indextmp << endl;
					}
				}
				ReferenceList(OriColIndex(numofoutband)) = 1;
			}
			else
			{
				tmpOut = Processer::RKLT(tmpIn, covList[groupReorder(index)], degree);
			}
			for (auto i = 0; i < tmpOut.cols(); i++, numofoutband++)
			{
				out.col(numofoutband) += tmpOut.col(i);
			}
		}

		// use only RKLT
		else if (!isUsePLSE && isUseRKLT)
		{
			MatrixXd tmpIn = MatrixXd::Zero(lengthPerBand, groupstate(groupReorder(index)));
			if (isVerbose)
			{
				cout << "==================================================================================" << endl;
				cout << "dealing with " << groupstate(groupReorder(index)) << " band : ";
			}

			int j = 0;
			for (auto i = 0; i < cluster.size(); i++)
			{
				if (cluster(i) == groupReorder(index))
				{
					if (isVerbose)
					{
						cout << i << " ";
					}
					tmpIn.col(j) += MatrixDoubleOri.col(i);
					OriColIndex(numofoutband + j) = i;
					j++;
				}
			}
			if (isVerbose)
			{
				cout << endl;
			}

			MatrixXd tmpOut;
			if (groupstate(groupReorder(index)) == 1)
			{
				tmpOut = tmpIn;
			}
			else
			{
				tmpOut = Processer::RKLT(tmpIn, covList[groupReorder(index)], degree);
			}
			for (auto i = 0; i < tmpOut.cols(); i++, numofoutband++)
			{
				out.col(numofoutband) += tmpOut.col(i);
			}
		}

		// don't use PLSE and RKLT either, just re-order the image as the cluster results
		else
		{
			MatrixXd tmpIn = MatrixXd::Zero(lengthPerBand, groupstate(groupReorder(index)));
			if (isVerbose)
			{
				cout << "==================================================================================" << endl;
				cout << "dealing with " << groupstate(groupReorder(index)) << " band : ";
			}

			int j = 0;
			for (auto i = 0; i < cluster.size(); i++)
			{
				if (cluster(i) == groupReorder(index))
				{
					if (isVerbose)
					{
						cout << i << " ";
					}

					tmpIn.col(j) += MatrixDoubleOriNoDemean.col(i);
					OriColIndex(numofoutband + j) = i;
					j++;
				}
			}
			if (isVerbose)
			{
				cout << endl;
			}

			MatrixXd tmpOut = tmpIn;
			for (auto i = 0; i < tmpOut.cols(); i++, numofoutband++)
			{
				out.col(numofoutband) += tmpOut.col(i);
			}
		}
	}
	if (isVerbose)
	{
		cout << "col index : " << (RowVectorXi)OriColIndex << endl;
	}

	/* write decorrelated data */
	ImageWriter imageWriter(oufn);
	imageWriter.WriteRAWdata(out.data(), ogeometry);

	/* write necessary information for reconstructing into side file */
	SyntaxWriter sw(aufn);
	sw.Write<uint16_t>(ogeometry.width);                        // width      : unsigned 16 bit (1~65536)
	sw.Write<uint16_t>(ogeometry.height);                       // height     : unsigned 16 bit (1~65536)
	sw.Write<uint16_t>(ogeometry.band);                         // band       : unsigned 16 bit (1~65536)
	sw.Write<uint8_t>(ogeometry.depth);                         // detph      : unsigned  8 bit (2^(1~65536))
	sw.Write<uint8_t>(ogeometry.sign ? 1 : 0);                  // sign       : unsigned  8 bit ( 1 bit is enough )
	sw.Write<uint8_t>(ogeometry.endianess ? 1 : 0);             // endianess  : unsigned  8 bit ( 1 bit is enough )
	sw.Write<uint8_t>(isPreprocess ? 1 : 0);                    // preprocess : unsigned  8 bit ( 1 bit is enough )
	sw.Write<uint8_t>(isUsePLSE ? 1 : 0);                       // usePLSE    : unsigned  8 bit ( 1 bit is enough )
	sw.Write<uint8_t>(isUseRKLT ? 1 : 0);                       // useRKLT    : unsigned  8 bit ( 1 bit is enough )
	if (isUsePLSE || isUseRKLT)
	{
		VectorXi meanTmp(igeometry.band);
		for (auto index = 0; index < igeometry.band; index++)
		{
			meanTmp(index) = meanOri(0, index);
		}
		sw.WriteArray<VectorXi, int32_t>(meanTmp);              // mean       :   signed 32 bit * band
	}
	sw.WriteArray<VectorXi, uint16_t>(cluster);                 // cluster    : unsigned 16 bit * band
	sw.WriteArray<VectorXd, double>(gamma);                     // gamma      :   double 64 bit * band
	sw.WriteArray<VectorXi, uint16_t>(OriColIndex);             // oricolindex: unsigned 16 bit * band
	if (isUsePLSE&&isUseRKLT)
	{
		for (auto index = 0; index < groupstate.size(); index++)
		{
			if (groupstate(index) != 1)
			{
				sw.WriteSymmetricMatrix<MatrixXd, double>(covList[index]); // cov matrix
			}
			else
			{
				sw.Write<uint16_t>(covList[index].rows());
				if (covList[index].rows()>1)
				{
					sw.WriteArray<VectorXd, double>(covList[index].topRows(covList[index].rows() - 1));
					sw.Write<uint16_t>(covList[index](covList[index].rows() - 1, 0));
				}
			}
		}
	}
	else if (!isUsePLSE && isUseRKLT)
	{
		for (auto index = 0; index < groupstate.size(); index++)
		{
			if (groupstate(index) != 1)
			{
				sw.WriteSymmetricMatrix<MatrixXd, double>(covList[index]);
			}
		}
	}
	sw.SyntaxWriterClose();

	if (isCheckLossless)
	{
		int pixelperframe = igeometry.width*igeometry.height;
		MatrixXd recon = out;
		int startCol = 0;
		for (auto index = 0; index < groupstate.size(); index++)
		{
			if (isUsePLSE && isUseRKLT)
			{
				MatrixXd recontmp;
				if (groupstate(groupReorder(index)) == 1)
				{
					int actualOrder = covList[groupReorder(index)].rows() - 1;
					if (actualOrder != 0)
					{
						recontmp = Processer::invPolynomialPredict(covList[groupReorder(index)].topRows(actualOrder), recon.col(startCol), recon.col(Processer::FindNewOrder(OriColIndex, (int)covList[groupReorder(index)](actualOrder))));
						recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))).setZero();
						recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))) += recontmp;
						startCol += groupstate(groupReorder(index));
					}
					else
					{
						startCol++;
					}
				}
				else
				{
					recontmp = Processer::invRKLT(recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))), covList[groupReorder(index)]);
					recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))).setZero();
					recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))) += recontmp;
					startCol += groupstate(groupReorder(index));
				}

			}
			else if (!isUsePLSE&&isUseRKLT)
			{
				MatrixXd recontmp;
				if (groupstate(groupReorder(index)) == 1)
				{
					startCol++;
				}
				else
				{
					recontmp = Processer::invRKLT(recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))), covList[groupReorder(index)]);
					recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))).setZero();
					recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))) += recontmp;
					startCol += groupstate(groupReorder(index));
				}
			}
		}

		MatrixXd reconReOrder = MatrixXd::Zero(pixelperframe, igeometry.band);
		if (isUsePLSE || isUseRKLT)
		{
			for (int col = 0; col < igeometry.band; col++)
			{
				reconReOrder.col(OriColIndex(col)) += (recon.col(col).array() + meanOri(0, OriColIndex(col))).matrix();
			}
		}
		else
		{
			for (int col = 0; col < igeometry.band; col++)
			{
				reconReOrder.col(OriColIndex(col)) += recon.col(col).matrix();
			}
		}
		MatrixXd resi = MatrixDoubleOriNoDemean - reconReOrder;
		if (resi.cwiseAbs().sum())
		{
			cout << "it is lossy"    << endl;
		}
		else
		{
			cout << "it is lossless" << endl;
		}
	}
	return SUCCESSIVE_EXEC;
}
