// BasicOpenCLApplication.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Chrono.h"
#include "RingBuffer.h"
#include <math.h>
#include <process.h>
#include <time.h>
#include <immintrin.h>


#define float8 __m256
#define set8(x) _mm256_set1_ps(x)
#define nbThreads 4
#define ringSize (8192)
const unsigned long long minimumJobSize = 274;
__m256 multi = _mm256_set1_ps(2);
__m256 max = _mm256_set1_ps(255);
float8 adder = _mm256_set1_ps(0.5);
float8 start = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
float8 ranges1, ranges2, rangec1, rangec2, dime[2];



void SaveBMP(char *fname, unsigned char *image, int width, int height, int componentPerPixel=1, int reverseColor=0)
{
	FILE *destination;
    int i,j;
	int *pt;
	char name[512],hdr[0x36];
	unsigned char *imsource=new unsigned char [width*height*3];
	//int al=(ImageSize*3)%4;
	
	if (componentPerPixel==1)
		for (i=0;i<width*height*3;i++)
			imsource[i]=image[i/3];
	else 
		for (i=0;i<width*height*3;i++)
			imsource[i]=image[i];
	if (reverseColor)
		for (j=0;j<height;j++)
			for (i=0;i<width;i++)
			{
				unsigned char aux;
				aux=imsource[3*(i+width*j)];
				imsource[3*(i+width*j)]=imsource[3*(i+width*j)+2];
				imsource[3*(i+width*j)+2]=aux;
			}
	strcpy(name,fname);
	i=(int)strlen(name);
	if (!((i>4)&&(name[i-4]=='.')&&(name[i-3]=='b')&&(name[i-2]=='m')&&(name[i-1]=='p')))
	{
		name[i]='.';
		name[i+1]='b';
		name[i+2]='m';
		name[i+3]='p';
		name[i+4]=0;
	}
	if ((destination=fopen(name, "wb"))==NULL) 
		perror("erreur de creation de fichier\n");
    hdr[0]='B';
    hdr[1]='M';
	pt=(int *)(hdr+2);// file size
	*pt=0x36+width*height*3;
	pt=(int *)(hdr+6);//reserved
	*pt=0x0;
	pt=(int *)(hdr+10);// image address
	*pt=0x36;
	pt=(int *)(hdr+14);// size of [0E-35]
	*pt=0x28;
	pt=(int *)(hdr+0x12);// Image width
	*pt=width;
	pt=(int *)(hdr+0x16);// Image heigth
	*pt=height;
	pt=(int *)(hdr+0x1a);// color planes
	*pt=1;
	pt=(int *)(hdr+0x1c);// bit per pixel
	*pt=24;
	for (i=0x1E;i<0x36;i++) 
		hdr[i]=0;
	fwrite(hdr,0x36,1,destination);
	fwrite (imsource,width*height*3,1,destination);
    fclose(destination);
	delete[] imsource;
}

typedef struct { float real; float im; } complex;

complex add(complex a, complex b)
{
	complex res;
	res.real = a.real + b.real;
	res.im = a.im + b.im;
	return res;
}
complex sub(complex a, complex b)
{
	complex res;
	res.real = a.real - b.real;
	res.im = a.im - b.im;
	return res;
}
complex mul(complex a, complex b)
{
	complex res;
	res.real = a.real*b.real - a.im*b.im;
	res.im = a.real*b.im + a.im*b.real;
	return res;
}

float squaredNorm(complex c)
{
	return c.real*c.real + c.im*c.im;
}

int Iterate(complex c)
{
	const int max_iterations = 255;
	complex z,a,l;
	a.real = 0.91;
	a.im = 0.;

	z = c;
	l.real = 4;
	l.im = 0;
	int i = 1;
	while (i < max_iterations)
	{
		z = mul(l, mul(z,sub(a,z)));
		if (squaredNorm(z) > 128)
			break;
		i += 2;
	}
	return (min(i, max_iterations));
}


__m256i Iterategroup(float8 real, float8 img)
{
	const int max_iterations = 255;
	__m256i maxi = _mm256_set1_epi32(255);
	__m256 zr,zi, ar,ai, lr, li, rs, is, rm, im, z;
	ar = _mm256_set1_ps(0.91);
	ai = _mm256_set1_ps(0.);
	__m256i trc = _mm256_set1_epi32(1);
	__m256i fac = _mm256_set1_epi32(0);
	__m256i trp = _mm256_set1_epi32(1);
	__m256 comparison = _mm256_set1_ps(128);
	__m256i j = _mm256_set1_epi32(1);


	zr = real;
	zi = img;
	lr = _mm256_set1_ps(4);
	li = _mm256_set1_ps(0);
	int i = 1;
	__m256i tr, trn;
	while (i < max_iterations)
	{

		rs = _mm256_sub_ps(ar, zr); is = _mm256_sub_ps(ai,zi);
		rm = _mm256_sub_ps(_mm256_mul_ps(zr, rs), _mm256_mul_ps(zi, is));
		im = _mm256_add_ps(_mm256_mul_ps(zr, is), _mm256_mul_ps(zi, rs));
		zr = _mm256_sub_ps(_mm256_mul_ps(lr, rm), _mm256_mul_ps(li, im));
		zi = _mm256_add_ps(_mm256_mul_ps(lr, im), _mm256_mul_ps(li, rm));
		z = _mm256_add_ps(_mm256_mul_ps(zr, zr), _mm256_mul_ps(zi, zi)); // z = mul(l, mul(z,sub(a,z))) from the old version
		tr = _mm256_max_epi32(_mm256_add_epi32(trc, _mm256_cvtps_epi32(_mm256_cmp_ps(z, comparison, _CMP_GT_OQ))),fac);//decides for each element if it should continue iterating 
		trn = _mm256_and_si256(_mm256_max_epi32(_mm256_cvtps_epi32(_mm256_div_ps(z, z)), fac) , _mm256_and_si256(trp, tr));// if it stops iterating it shouldn't start reiterating (1 for continue, 0 for stop)
		trp = tr;
		j = _mm256_add_epi32(_mm256_add_epi32(trn,trn), j);
		if (_mm256_testz_si256(trn, trc) == 1) {
			break;
		}
		i += 2;
	}
	return _mm256_min_epi32(j,maxi);
}



void SimpleFractalDrawing(unsigned char *image, int dim[2],float range[2][2])
{
	Chrono c;
	for (int j=0;j<dim[1];j++)
		for (int i=0;i<dim[0];i++)
		{
			complex c;
			c.real=range[0][0]+(i+0.5)*(range[0][1]-range[0][0])/dim[0]; //Create x coordinates within the range [range[0][0] .. range[0][1]]
			c.im=range[1][0]+(j+0.5)*(range[1][1]-range[1][0])/dim[1]; //Create x coordinates within the range [range[1][0] .. range[1][1]] 
			float f = 2 * Iterate(c);
			if (f > 255.)
				f = 255.;
			image[j*dim[0]+i]=f; 
		}
	c.PrintElapsedTime_ms("time CPU (s): ");

}


void SimpleFractalDrawingSIMD(unsigned char* image, int dim[2], float range[2][2])
{
	Chrono c;
	float8 ranges1 = _mm256_set1_ps(range[0][1] - range[0][0]);
	float8 ranges2 = _mm256_set1_ps(range[1][1] - range[1][0]);
	float8 rangec1 = _mm256_set1_ps(range[0][0]);
	float8 rangec2 = _mm256_set1_ps(range[1][0]);
	float8 adder = _mm256_set1_ps(0.5);
	float8 dime[2] = { _mm256_set1_ps(dim[0]),_mm256_set1_ps(dim[1]) };
	float8 start = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
	for (int j = 0; j < dim[1]; j += 8)
		for (int i = 0; i < dim[0]; i += 8)
		{
			float8 loopi = _mm256_set1_ps(i);
			float8 loopj = _mm256_set1_ps(j);
			float8 real = _mm256_add_ps(_mm256_div_ps(_mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(start, loopi), adder), ranges1), dime[0]), rangec1);
			float8 img = _mm256_add_ps(_mm256_div_ps(_mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(start, loopj), adder), ranges2), dime[1]), rangec2);
			float8 f[8];
			f[0] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[0])))), max);
			f[1] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[1])))), max);
			f[2] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[2])))), max);
			f[3] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[3])))), max);
			f[4] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[4])))), max);
			f[5] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[5])))), max);
			f[6] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[6])))), max);
			f[7] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[7])))), max);
			for (int k = 0; k < 8; k++)
			{
				image[(j + k) * dim[0] + i] = f[k].m256_f32[0];
				image[(j + k) * dim[0] + i + 1] = f[k].m256_f32[1];
				image[(j + k) * dim[0] + i + 2] = f[k].m256_f32[2];
				image[(j + k) * dim[0] + i + 3] = f[k].m256_f32[3];
				image[(j + k) * dim[0] + i + 4] = f[k].m256_f32[4];
				image[(j + k) * dim[0] + i + 5] = f[k].m256_f32[5];
				image[(j + k) * dim[0] + i + 6] = f[k].m256_f32[6];
				image[(j + k) * dim[0] + i + 7] = f[k].m256_f32[7];
			}
		}
	c.PrintElapsedTime_ms("time CPU (s): ");
}



void ThreadCal(unsigned char* image, int dim[2], int i, int j) {
	float8 loopi = _mm256_set1_ps(i);
	float8 loopj = _mm256_set1_ps(j);
	float8 real = _mm256_add_ps(_mm256_div_ps(_mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(start, loopi), adder), ranges1), dime[0]), rangec1);
	float8 img = _mm256_add_ps(_mm256_div_ps(_mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(start, loopj), adder), ranges2), dime[1]), rangec2);
	float8 f[8];
	f[0] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[0])))), max);
	f[1] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[1])))), max);
	f[2] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[2])))), max);
	f[3] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[3])))), max);
	f[4] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[4])))), max);
	f[5] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[5])))), max);
	f[6] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[6])))), max);
	f[7] = _mm256_min_ps(_mm256_mul_ps(multi, _mm256_cvtepi32_ps(Iterategroup(real, _mm256_set1_ps(img.m256_f32[7])))), max);
	for (int k = 0; k < 8; k++)
	{
		image[(j + k) * dim[0] + i] = f[k].m256_f32[0];
		image[(j + k) * dim[0] + i + 1] = f[k].m256_f32[1];
		image[(j + k) * dim[0] + i + 2] = f[k].m256_f32[2];
		image[(j + k) * dim[0] + i + 3] = f[k].m256_f32[3];
		image[(j + k) * dim[0] + i + 4] = f[k].m256_f32[4];
		image[(j + k) * dim[0] + i + 5] = f[k].m256_f32[5];
		image[(j + k) * dim[0] + i + 6] = f[k].m256_f32[6];
		image[(j + k) * dim[0] + i + 7] = f[k].m256_f32[7];
	}
}

void ThreadJob(unsigned char *image, int dim[2],float range[2][2] ,int begin, int end, RingBuffer* rb)
{
	
	ranges1 = _mm256_set1_ps(range[0][1] - range[0][0]);
	ranges2 = _mm256_set1_ps(range[1][1] - range[1][0]);
	rangec1 = _mm256_set1_ps(range[0][0]);
	rangec2 = _mm256_set1_ps(range[1][0]);
	dime[0] = _mm256_set1_ps(dim[0]);
	dime[1] = _mm256_set1_ps(dim[1]);
	for (int j = begin ; j < end; j+=8){
		for (int i = 0; i < (dim[0]); i+=8)
		{
			if ((i + 8) < 273) {
				ThreadCal(image, dim, i, j); // the side dark section is done straight away to save time
			}
			else {
				int val[2] = { i,j };
				if (!rb->AddToRing(val)) {
					ThreadCal(image, dim, i, j);
				}
			}
		}
	}
	while (!rb->IsRingEmptyCheck())
	{
		int *v = rb->RemoveFromRing();
		if (v != NULL) {
			ThreadCal(image, dim, *v, *(v + 1));
		}
	}
}

void SimpleFractalDrawingSIMD_MT(unsigned char* image, int dim[2], float range[2][2])
{
	Chrono c;
	std::thread t[nbThreads-1];
	RingBuffer rb(ringSize);
	int split = dim[1] / nbThreads;
	for (int i = 0; i < nbThreads; i++)
	{
		if (i < nbThreads - 1){
			t[i] = std::thread(ThreadJob, image, dim, range, split * i, split * (i + 1), &rb);
		}
		else {
			ThreadJob( image, dim, range, split * i, split * (i + 1), &rb);
		}
	}
	for (int i = 0; i < nbThreads-1; i++) {
		t[i].join();
	}
	c.PrintElapsedTime_ms("time CPU (s): ");
}


int main(int argc, char* argv[])
{
	int dims[2]={1024,1024};
	float range[2][2] = { {-0.003,0.008},{-0.0002,0.0005} };
	float range2[2][2] = { {-1.4,0.6},{-1.1,1.3} };
	unsigned char *image=new unsigned char[dims[0]*dims[1]];
	SimpleFractalDrawing(image,dims,range); //largest 64bit prime
	SaveBMP("fractal.bmp", image, dims[0], dims[1]);
	for (int i = 0; i < dims[0] * dims[1]; i++)
		image[i] = 127; //resetting image to grey
	SimpleFractalDrawingSIMD(image, dims, range);
	SaveBMP("fractalSIMD.bmp", image, dims[0], dims[1]);
	SimpleFractalDrawingSIMD_MT(image, dims,range);
	SaveBMP("fractalSIMD_MT.bmp",image,dims[0],dims[1]);
	delete[] image;
	return 0; 
}

/*
CPU 1 thread : 123	ms
SIMD_1_Thread:	16	ms
SIMD_MT		:	7	ms


*/