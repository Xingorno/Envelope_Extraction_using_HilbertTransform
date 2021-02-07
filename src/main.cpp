#include <iostream>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex> 
using namespace std;

//macros for real and imaginary parts
#define REAL 0
#define IMAG 1
//length of complex array
#define N 19

void hilbert(const double* in, fftw_complex* out)
{
    // copy the data to the complex array
    for (int i = 0; i < N; ++i) {
        out[i][REAL] = in[i];
        out[i][IMAG] = 0;
    }
    // creat a DFT plan and execute it
    fftw_plan plan = fftw_plan_dft_1d(N, out, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    // destroy a plan to prevent memory leak
    fftw_destroy_plan(plan);
    int hN = N >> 1; // half of the length (N/2)
    int numRem = hN; // the number of remaining elements
    // multiply the appropriate value by 2
    //(those should multiplied by 1 are left intact because they wouldn't change)
    for (int i = 1; i < hN; ++i) {
        out[i][REAL] *= 2;
        out[i][IMAG] *= 2;
    }
    // if the length is even, the number of the remaining elements decrease by 1
    if (N % 2 == 0)
        numRem--;
    else if (N > 1) {
        out[hN][REAL] *= 2;
        out[hN][IMAG] *= 2;
    }
    // set the remaining value to 0
    // (multiplying by 0 gives 0, so we don't care about the multiplicands)
    memset(&out[hN + 1][REAL], 0, numRem * sizeof(fftw_complex));
    // creat a IDFT plan and execute it
    plan = fftw_plan_dft_1d(N, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    // do some cleaning
    fftw_destroy_plan(plan);
    fftw_cleanup();
    // scale the IDFT output
    for (int i = 0; i < N; ++i) {
        out[i][REAL] /= N;
        out[i][IMAG] /= N;
    }




}
/* Displays complex numbers in the form a +/- bi. */
void displayComplex(fftw_complex* y)
{
    for (int i = 0; i < N; i++) {
        if (y[i][IMAG] < 0)
            cout << y[i][REAL] << "-" << abs(y[i][IMAG]) << "i" << endl;
        else
            cout << y[i][REAL] << "+" << y[i][IMAG] << "i" << endl;
    }
}


/* Test */
int main()
{
    
 
    // input array
    double x[N];
    // output array
    fftw_complex y[N];
    // fill the first of some numbers
    for (int i = 0; i < N; ++i) {
        x[i] = i + 1;  // i.e.{1 2 3 4 5 6 7 8}
    }
    // compute the FFT of x and store the result in y.
    hilbert(x, y);
    // display the result
    cout << "Hilbert =" << endl;
    displayComplex(y);
    double envelope[N];
    for (int i = 0; i < N; i++)
    {
        complex<double> mycomplex(y[i][REAL], y[i][IMAG]);
        envelope[i] = abs(mycomplex);
        cout << envelope[i]<<"\n" << endl;
    }
    
 
}