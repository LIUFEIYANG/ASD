/*================================================================
*   Copyright (C) 2017 IIP. All rights reserved.
*
*   filename: Common.h
*   creator : Feiyang Liu
*   date    : 2017-09-05
*   describe:
*
================================================================*/

#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <iostream>

// macro function
#define L_squa(x)                (x*x)
#define L_cube(x)                (x*x*x)
#define L_max(a,b)               (a>b?a:b)
#define L_min(a,b)               (a>b?b:a)
#define L_clip(a,bitdepth)       (L_min(L_max(a,0),(1<<bitdepth-1))
#define L_clip3(a,minval,maxval) (L_min(L_max(a,minval),maxval)

#define DEBUG_MODE 1
#if DEBUG_MODE
#define DEBUG_INFO(msg)          cout << msg << endl ;
#else
#define DEBUG_INFO(msg)
#endif

// error type
#define SUCCESSIVE_EXEC      0

#define ERROR_MEMORY_MALLOC -1
#define ERROR_FILE_OPEN     -2
#define ERROR_FILE_SIZE     -3
#define ERROR_PARA_CFG      -4
#define ERROR_NONSENSE_CFG  -5

#define ERROR_MEM_MALLOC   -11

// some parameters
#define MAX_ORDER           4

// common functions
inline void printCmdLine(int argv, char* argc[]){
	std::cout << "running : ";
	for (int index = 0; index<argv; index++){
		std::cout << argc[index] << " ";
	}
	std::cout << std::endl;
}

struct image_geometry
{
	int32_t width;
	int32_t height;
	int32_t band;
	int32_t depth;
	bool    sign;
	bool    endianess;
};

struct point3d
{
	int32_t x;
	int32_t y;
	int32_t z;
};

struct size3d
{
	int32_t x_length;
	int32_t y_length;
	int32_t z_length;
};

#endif
