/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: ImageReader.h
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#ifndef IMAGE_READER_H
#define IMAGE_READER_H

#include <string>
#include "Common.h"
using namespace std;

class ImageReader
{
private:
	string datafilename;

public:
	ImageReader(string filename)
	{
		datafilename = filename;
	};
	~ImageReader(){};

	double * ReadRAWdata(image_geometry geo);


private:
	void xTransform2double(double * dst, void * src, image_geometry geo);
	template<typename T> void EndianConvert(T * src, int32_t length = 1);
};

#endif

