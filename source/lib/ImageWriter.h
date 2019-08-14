/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: ImageWriter.h
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#ifndef IMAGE_WRITER_H
#define IMAGE_WRITER_H

#include <string>
#include "Common.h"
using namespace std;

class ImageWriter
{
private:
	string datafilename;

public:
	ImageWriter(string filename)
	{
		datafilename=filename;
	}
	~ImageWriter(){};

	void WriteRAWdata(double * data, image_geometry geo);

private:
	template<typename T> T * xTransformFromDouble(double * src, size_t length);
	template<typename T> void EndianConvert(T * src, int32_t length = 1);
};

#endif
