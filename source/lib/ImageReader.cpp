/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: ImageReader.cpp
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#include "ImageReader.h"
#include "Common.h"
#include <iostream>
#include <fstream>
#include <cassert>
using namespace std;

template<typename T> void ImageReader::EndianConvert(T * src, int32_t length)
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

double * ImageReader::ReadRAWdata(image_geometry geo)
{
	size_t length = (size_t)geo.width*geo.height*geo.band;
	size_t lengthByByte = length * (geo.depth >> 3);
	uint8_t * data = nullptr;
	try
	{
		data = new uint8_t[lengthByByte];
	}
	catch (const bad_alloc& e)
	{
		cout << e.what() << endl;
		cout << "fail to malloc [" << lengthByByte << "] byte memory" << endl;
		exit(ERROR_MEM_MALLOC);
	}

	fstream infile(datafilename,ios::in|ios::binary);
	if (!infile.is_open())
	{
		cout << "Fail to open input file : " << datafilename << endl;
		exit(ERROR_FILE_OPEN);
	}
	infile.read((char *)data, lengthByByte);
	if (lengthByByte!=infile.gcount())
	{
		cout << "file is too small" << endl;
		exit(ERROR_FILE_SIZE);
	}
	infile.close();

	double *dataDouble = nullptr;
	try
	{
		dataDouble = new double[length];
	}
	catch (const bad_alloc& e)
	{
		cout << e.what() << endl;
		cout << "fail to malloc [" << length << "] double memory" << endl;
		exit(ERROR_MEM_MALLOC);
	}
	xTransform2double(dataDouble,data,geo);
	delete[] data;
	data = nullptr;
	return dataDouble;
}

void ImageReader::xTransform2double(double * dst, void * src, image_geometry geo)
{
	int32_t length = geo.width*geo.height*geo.band;
	if (geo.depth == 8)
	{
		if (geo.sign)
		{
			/*signed 8 bit-depth [-127,128] */
			int8_t *tmp = (int8_t*)src;
			for (size_t index = 0; index<length; index++)
			{
				if (geo.endianess) EndianConvert(tmp+index);
				dst[index] = (double)tmp[index];
			}
		}
		else
		{
			/*unsigned 8 bit-depth [0,255] */
			uint8_t *tmp = (uint8_t*)src;
			for (size_t index = 0; index<length; index++)
			{
				if (geo.endianess) EndianConvert(tmp + index);
				dst[index] = (double)tmp[index];
			}
		}
	}
	else if (geo.depth == 16)
	{
		if (geo.sign)
		{
			/*signed 16 bit-depth [-32767,32768] */
			int16_t *tmp = (int16_t*)src;
			for (size_t index = 0; index<length; index++)
			{
				if (geo.endianess) EndianConvert(tmp + index);
				dst[index] = (double)tmp[index];
			}
		}
		else
		{
			/*unsigned 16 bit-depth [0,65535] */
			uint16_t *tmp = (uint16_t*)src;
			for (size_t index = 0; index<length; index++)
			{
				if (geo.endianess) EndianConvert(tmp + index);
				dst[index] = (double)tmp[index];
			}
		}
	}
	else if (geo.depth == 32)
	{
		if (geo.sign)
		{
			/*signed 32 bit-depth */
			int32_t *tmp = (int32_t*)src;
			for (size_t index = 0; index<length; index++)
			{
				if (geo.endianess) EndianConvert(tmp + index);
				dst[index] = (double)tmp[index];
			}
		}
		else
		{
			/*unsigned 32 bit-depth */
			uint32_t *tmp = (uint32_t*)src;
			for (size_t index = 0; index<length; index++)
			{
				if (geo.endianess) EndianConvert(tmp + index);
				dst[index] = (double)tmp[index];
			}
		}
	}
}
