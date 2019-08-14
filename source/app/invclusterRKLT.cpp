#include "../lib/Processer.h"
#include "../lib/SyntaxReader.h"
#include "../lib/Cluster.h"
#include "../lib/ImageReader.h"
#include "../lib/ImageWriter.h"
#include "../lib/Predicter.h"
#include "../lib/Common.h"
#include "../lib/cmdline.h"
#include <fstream>
using namespace std;

void setup_parameter(cmdline::parser &parameters){
	// mandatory
	parameters.add<string>("input"     , 'i', "input raw image file"                   , true, "");
	parameters.add<string>("output"    , 'o', "output transformed file"                , true, "");
	parameters.add<string>("auxi"      , 'a', "auxilary data to recover original image", true, "");
	parameters.add<string>("outputsize", 'f', "size of output raw image file"          , true, "");
	return;
}

image_geometry parse_size(string size) {
	image_geometry geo;
	const char * delim = ":";
	vector<string> tmp;
	char *p = nullptr;
	p = strtok(const_cast<char*>(size.c_str()), delim);
	while (p) {
		tmp.push_back(string(p));
		p = strtok(NULL, delim);
	}
	assert(tmp.size() == 6);
	geo.width = atoi(tmp[0].c_str());
	geo.height = atoi(tmp[1].c_str());
	geo.band = atoi(tmp[2].c_str());
	geo.depth = atoi(tmp[3].c_str());
	geo.sign = (bool)atoi(tmp[4].c_str());
	geo.endianess = (bool)atoi(tmp[5].c_str());
	return geo;
}

int main(int argc, char* argv[])
{
	// parse command line
	printCmdLine(argc, argv);
	cmdline::parser parameters;
	setup_parameter(parameters);
	parameters.parse_check(argc,argv);

	Cluster cls;
	//setup parameters from cmdline parser
	string infn              = parameters.get<string>("input");
	string oufn              = parameters.get<string>("output");
	string aufn              = parameters.get<string>("auxi");
	image_geometry ogeometry = parse_size(parameters.get<string>("outputsize"));

	// read auxilary information for backward transform
	SyntaxReader sr(aufn);
	image_geometry igeometry;
	igeometry.width     = sr.Read<uint16_t>();
	igeometry.height    = sr.Read<uint16_t>();
	igeometry.band      = sr.Read<uint16_t>();
	igeometry.depth     = sr.Read<uint8_t>();
	igeometry.sign      = (sr.Read<uint8_t>() == 1);
	igeometry.endianess = (sr.Read<uint8_t>() == 1);
	assert(igeometry.width == ogeometry.width && igeometry.height == ogeometry.height && igeometry.band == ogeometry.band);

	bool isPreprocess       = (sr.Read<uint8_t>() == 1);
	bool isUsePLSE          = (sr.Read<uint8_t>() == 1);
	bool isUseRKLT          = (sr.Read<uint8_t>() == 1);
	VectorXi meanOri;
	if (isUsePLSE || isUseRKLT)
		meanOri = sr.ReadArray<VectorXi, int32_t>(igeometry.band);
	VectorXi cluster      = sr.ReadArray<VectorXi,uint16_t>(igeometry.band);
	VectorXd gamma        = sr.ReadArray<VectorXd,double>(igeometry.band);
	VectorXi OriColIndex  = sr.ReadArray<VectorXi,uint16_t>(igeometry.band);
	VectorXi groupstate   = cls.GetClusterGroupState(cluster);
	VectorXi groupReorder = cls.ReOrderClusterByNum(cluster,gamma); 
	MatrixXd *covList     = new MatrixXd[groupstate.size()];
	if (isUsePLSE && isUseRKLT) {
		for (auto i = 0; i<groupstate.size(); i++) {
			if (groupstate(i) != 1) {
				covList[i] = sr.ReadSymmetricMatrix<MatrixXd, double>(groupstate(i), groupstate(i));
			}
			else {
				int order = sr.Read<uint16_t>();
				covList[i] = MatrixXd::Zero(order, 1);
				if (order>1) {
					covList[i].topRows(order - 1) += (MatrixXd)sr.ReadArray<VectorXd, double>(order - 1);
					covList[i](order - 1, 0) = sr.Read<uint16_t>();
				} else {
					covList[i](order - 1, 0) = -1;
				}
			}
		}
	} else if (!isUsePLSE && isUseRKLT) {
		for (auto i = 0; i < groupstate.size(); i++) {
			if (groupstate(i) != 1)
				covList[i] = sr.ReadSymmetricMatrix<MatrixXd, double>(groupstate(i), groupstate(i));
		}
	}
	sr.SyntaxReaderClose();

	cout << "inputfilename   :" << infn << endl;
	cout << "outputfilename  :" << oufn << endl;
	cout << "auxilaryfilename:" << aufn << endl;
	cout << "input geometry   " << endl;
	cout << "-size           :" << igeometry.height << "*" << igeometry.width << "*" << igeometry.band << endl;
	cout << "-depth          :" << igeometry.depth << endl;
	cout << "-sign           :" << igeometry.sign << endl;
	cout << "-endianess      :" << (igeometry.endianess ? "big endian" : "little endian") << endl;
	cout << "output geometry  " << endl;
	cout << "-size           :" << ogeometry.height << "*" << ogeometry.width << "*" << ogeometry.band << endl;
	cout << "-depth          :" << ogeometry.depth << endl;
	cout << "-sign           :" << ogeometry.sign << endl;
	cout << "-endianess      :" << (ogeometry.endianess ? "big endian" : "little endian") << endl;

	cout << "preprocess      :" << (isPreprocess ? "yes" : "no") << endl;
	cout << "use PLSE        :" << (isUsePLSE ? "yes" : "no") << endl;
	cout << "use RKLT        :" << (isUseRKLT ? "yes" : "no") << endl;

	cout << "cluster         : " << (RowVectorXi)cluster << endl;
	cout << "orderindex      : " << (RowVectorXi)OriColIndex << endl;
	cout << "groupstate      : " << (RowVectorXi)groupstate << endl;
	cout << "group rank      : " << (RowVectorXi)groupReorder << endl;

	int lengthPerBand = igeometry.width*igeometry.height;
	//MatrixXd recon = out;
	ImageReader imageReader(infn);
	double *dataDouble = imageReader.ReadRAWdata(igeometry);
	Map<MatrixXd> recon(dataDouble,lengthPerBand, igeometry.band);

	int pixelperframe = igeometry.width*igeometry.height;
	int startCol = 0;
	for (auto index = 0; index<groupstate.size(); index++) {
		if (isUsePLSE && isUseRKLT) {
			MatrixXd recontmp;
			if (groupstate(groupReorder(index)) == 1) {
				int actualOrder = covList[groupReorder(index)].rows() - 1;
				if (actualOrder != 0) {
					recontmp = Processer::invPolynomialPredict(covList[groupReorder(index)].topRows(actualOrder), recon.col(startCol), recon.col(Processer::FindNewOrder(OriColIndex, (int)covList[groupReorder(index)](actualOrder))));
					recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))).setZero();
					recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))) += recontmp;
					startCol += groupstate(groupReorder(index));
				} else {
					startCol++;
				}
			} else {
				recontmp = Processer::invRKLT(recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))), covList[groupReorder(index)]);
				recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))).setZero();
				recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))) += recontmp;
				startCol += groupstate(groupReorder(index));
			}

		} else if (!isUsePLSE&&isUseRKLT) {
			MatrixXd recontmp;
			if (groupstate(groupReorder(index)) == 1) {
				startCol++;
			}
			else {
				recontmp = Processer::invRKLT(recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))), covList[groupReorder(index)]);
				recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))).setZero();
				recon.block(0, startCol, pixelperframe, groupstate(groupReorder(index))) += recontmp;
				startCol += groupstate(groupReorder(index));
			}
		}
	}

	MatrixXd reconReOrder = MatrixXd::Zero( lengthPerBand, igeometry.band);
	if (isUsePLSE || isUseRKLT) {
		for (int col = 0; col < igeometry.band; col++)
			reconReOrder.col(OriColIndex(col)) += (recon.col(col).array() + meanOri(OriColIndex(col))).matrix();
	} else {
		for (int col = 0; col < igeometry.band; col++)
			reconReOrder.col(OriColIndex(col)) += recon.col(col).matrix();
	}
	if (isPreprocess)
		Processer::invPreProcess(reconReOrder);
	ImageWriter imageWriter(oufn);
	imageWriter.WriteRAWdata(reconReOrder.data(),ogeometry);
	return SUCCESSIVE_EXEC;
}
