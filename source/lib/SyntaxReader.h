/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: SyntaxReader.h
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#ifndef SYNTAX_READER_H
#define SYNTAX_READER_H

#include <cstdint>
#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
using namespace Eigen;
using namespace std;

class SyntaxReader
{
private:
	fstream sideinfo;

	template<class T> T xRead()
	{
		T data;
		sideinfo.read((char *)&data, sizeof(T));
		return data;
	}

	template<class T> T * xRead(int length)
	{
		T *data = new T[length];
		sideinfo.read((char *)data, sizeof(T)*length);
		return data;
	}

public:
	SyntaxReader(string filename)
	{
		sideinfo.open(filename, ios::binary | ios::in);
		if (!sideinfo.is_open())
		{
			cout << "Fail to open file " << filename.c_str() << endl;
		}
	}

	~SyntaxReader()
	{
		if (sideinfo.is_open())
		{
			sideinfo.close();
		}
	}

	void SyntaxReaderClose()
	{
		if (sideinfo.is_open())
		{
			sideinfo.close();
		}
	}

	template<class T> T Read()
	{
		return xRead<T>();
	}

	template<class T1, class T2> T1 ReadArray(int length)
	{
		T2 * tmp = xRead<T2>(length);
		T1 res(length);
		for (auto i = 0; i < length; i++)
		{
			res(i) = tmp[i];
		}
		return res;
	}

	template<class T1, class T2> T1 ReadMatrix(int rows, int cols) {
		T2 * tmp = xRead<T2>(rows*cols);
		T1 res(rows, cols);
		int index = 0;
		for (auto j = 0; j < cols; j++)
		{
			for (auto i = 0; i < rows; i++)
			{
				res(i, j) = tmp[index];
				index++;
			}
		}
		return res;
	}

	template<class T1, class T2> T1 ReadTriMatrix(int rows, int cols) {
		int length = (rows*(rows + 1)) >> 1;
		T2 * tmp = xRead<T2>(length);
		T1 res(rows, cols);
		res.setZero();
		int index = 0;
		for (auto j = 0; j < cols; j++)
		{
			for (auto i = j; i < rows; i++)
			{
				res(i, j) = tmp[index];
				index++;
			}
		}
		return res;
	}

	template<class T1, class T2> T1 ReadSymmetricMatrix(int rows, int cols) {
		T1 res = ReadTriMatrix<T1, T2>(rows, cols);
		for (auto i = 0; i < rows; i++)
		{
			for (auto j = i + 1; j < cols; j++)
			{
				res(i, j) = res(j, i);
			}
		}
		return res;
	}
};

#endif
