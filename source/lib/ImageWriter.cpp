/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: ImageWriter.cpp
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#include "ImageWriter.h"
#include "Common.h"
#include <iostream>
#include <fstream>
#include <assert.h>
#include <cstdint>
using namespace std;

template<typename T> void ImageWriter::EndianConvert(T * src, int32_t length)
{
	int step = sizeof(T);
	char * start = (char*)src;
	if (step == 1)
	{
		return;
	}
	else if (step == 2)
	{
		for (auto i = 0; i < length; i++)
		{
			swap(start[0], start[1]);
			start += step;
		}
	}
	else if (step == 4)
	{
		for (auto i = 0; i < length; i++)
		{
			swap(start[0], start[3]);
			swap(start[1], start[2]);
			start += step;
		}
	}
	else if (step == 8)
	{
		for (auto i = 0; i < length; i++)
		{
			swap(start[0], start[7]);
			swap(start[1], start[6]);
			swap(start[2], start[5]);
			swap(start[3], start[4]);
			start += step;
		}
	}
	return;
}

void ImageWriter::WriteRAWdata(double * data, image_geometry geo)
{
	size_t length = (size_t)geo.width*geo.height*geo.band;
	size_t lengthByByte = length * (geo.depth >> 3);

	void *output=nullptr;
	if (geo.depth == 8)
	{
		if (geo.sign)
		{
			output=xTransformFromDouble<int8_t>(data,length);
			if (geo.endianess)
			{
				EndianConvert((int8_t *)output, length);
			}
		}
		else
		{
			output=xTransformFromDouble<uint8_t>(data,length);
			if (geo.endianess)
			{
				EndianConvert((uint8_t *)output, length);
			}
		}
	}
	else if(geo.depth == 16)
	{
		if (geo.sign)
		{
			output=xTransformFromDouble<int16_t>(data,length);
			if (geo.endianess)
			{
				EndianConvert((int16_t *)output, length);
			}
		}
		else
		{
			output=xTransformFromDouble<uint16_t>(data,length);
			if (geo.endianess)
			{
				EndianConvert((uint16_t *)output, length);
			}
		}
	}
	else if (geo.depth == 32)
	{
		if (geo.sign)
		{
			output = xTransformFromDouble<int32_t>(data, length);
			if (geo.endianess)
			{
				EndianConvert((int32_t *)output, length);
			}
		}
		else {
			output = xTransformFromDouble<uint32_t>(data, length);
			if (geo.endianess)
			{
				EndianConvert((uint32_t *)output, length);
			}
		}
	}

	fstream outfile(datafilename, ios::out | ios::binary);
	if (!outfile.is_open())
	{
		cout << "Fail to open output file : " << datafilename << endl;
		exit(ERROR_FILE_OPEN);
	}
	outfile.write((char*)output, lengthByByte);
	outfile.close();
}

template<typename T>
T* ImageWriter::xTransformFromDouble(double * src, size_t length)
{
	T* dst=nullptr;
	try
	{
		dst = new T[length];
	}
	catch (const bad_alloc& e)
	{
		cout << e.what() << endl;
		cout << "fail to malloc [" << length << "] memory" << endl;
		exit(ERROR_MEM_MALLOC);
	}
	for (size_t index = 0; index<length; index++)
	{
		dst[index] = (T)src[index];
	}
	return dst;
}
