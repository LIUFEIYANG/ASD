/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: SyntaxWriter.h
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#ifndef SYNTAX_WRITER_H
#define SYNTAX_WRITER_H

#include <cstdint>
#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
using namespace Eigen;
using namespace std;

class SyntaxWriter
{
private:
	fstream sideinfo;

	template<class T> void xWrite(T data)
	{
		sideinfo.write((char *)&data, sizeof(T));
	}

	template<class T> void xWrite(T * data, size_t length)
	{
		sideinfo.write((char *)data, sizeof(T)*length);
	}

public:
	SyntaxWriter(string filename)
	{
		sideinfo.open(filename, ios::binary | ios::out);
		if (!sideinfo.is_open())
		{
			cout << "Fail to open file " << filename.c_str() << endl;
		}
	}

	~SyntaxWriter()
	{
		if (sideinfo.is_open())
		{
			sideinfo.close();
		}
	}

	void SyntaxWriterClose()
	{
		if (sideinfo.is_open())
		{
			sideinfo.close();
		}
	}

	template<class T> void Write(T data)
	{
		xWrite(data);
	}

	template<class T1, class T2> void WriteArray(T1 data)
	{
		if (data.size())
		{
			T2 * tmp = new T2[data.size()];
			for (auto i = 0; i < data.size(); i++)
			{
				tmp[i]=(T2)data(i);
			}

			xWrite(tmp, data.size());
		}
	}

	template<class T1, class T2> void WriteMatrix(T1 data)
	{
		if(data.rows()&&data.cols())
		{
			T2 * tmp = new T2[data.rows()*data.cols()];
			int index = 0;
			for (auto j = 0; j < data.cols(); j++)
			{
				for (auto i = 0; i < data.rows(); i++)
				{
					tmp[index] = (T2)data(i, j);
					index++;
				}
			}
			xWrite(tmp, data.rows()*data.cols());
		}
	}

	template<class T1, class T2> void WriteTriMatrix(T1 data)
	{
		if (data.rows() && data.cols())
		{
			int length = (data.rows()*(data.rows() + 1)) >> 1;
			T2 * tmp = new T2[length];
			int index = 0;
			for (auto j = 0; j < data.cols(); j++)
			{
				for (auto i = j; i < data.rows(); i++)
				{
					tmp[index] = (T2)data(i, j);
					index++;
				}
			}
			xWrite(tmp, length);
		}
	}

	template<class T1, class T2> void WriteSymmetricMatrix(T1 data)
	{
		WriteTriMatrix<T1, T2>(data);
	}
};

#endif
