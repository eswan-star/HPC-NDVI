#include <algorithm>
#include "papi.h"
#include <numeric>
#include <fstream>
#include <x86intrin.h>
#include <stdalign.h>
#include <cstring>
#include <cmath>
#include <iostream>

void get_Pk_kernel(float* __restrict__ img_data, const int &Nt, const int &Nr, const int &Nc);
void get_Pk_kernel_serial(float* __restrict__ img_data, const int &Nt, const int &Nr, const int &Nc);
void left_transform(float* __restrict__ img_data, const int &Nt, const int&Nr, const int &Nc);
void right_transform(float* __restrict__ img_data, const int &Nt, const int &Nr, const int &Nc);

