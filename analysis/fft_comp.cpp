#ifndef fft_test
#define fft_test


#include <iostream>
#include <complex>
#include <string>
#include <vector>
#include <fftw3.h>
#include <fstream>

using namespace std;

int main()
{
    /* First load in the data array that contains the data to be fourier transformed */   
    string filename = "fft_data.dat";
    ifstream dataloader(filename.c_str());
    int n;
    double *data_in = (double *)fftw_malloc(n*sizeof(double));    /* The size of the output array is n/2 + 1*/
    //complex<double>  *data_out = (complex<double> *)fftw_malloc((n/2 + 1)*sizeof(complex<double> ));
    complex<double> *data_out = (complex<double> *)fftw_malloc(n*sizeof(complex<double>));

    /* Now we can compute the actual discrete fourier transform with n samples */
    fftw_plan plan = fftw_plan_dft_r2c_1d(n, data_in, reinterpret_cast<fftw_complex*>(data_out), FFTW_MEASURE);


    if (dataloader.is_open())
    {
        dataloader >> n;
        cout << "The no of data points: " << n << endl;
        for (int i =0 ; i < n; i++)
        {
            dataloader >> data_in[i];
        }
        dataloader.close();
    }
    else
    {
        cout << "Failed to open this file!" << endl;
    }

    cout << "The first data point is " << data_in[0] << endl;
    cout << "The last data point is " << data_in[n-1] << endl;
    cout << "Data successfully loaded!!!" << endl;

    cout << "Starting to compute the DFT " << endl;
    cout << "Plan completed. Now for the execution!!!!" << endl;
    fftw_execute(plan);
    cout << "DFT successfully computed. Writing the results to file. " << endl;
    
    cout << "The first data point output is " << data_out[0] << endl;
    cout << "The last data point output is " << data_out[n-1] << endl;
    /* We can now write out the data file of the fourier transform to a file */
    string outputfilename = "fft_results.dat";
    ofstream datasaver(outputfilename.c_str());

    if (datasaver.is_open())
    {
        for (int i = 0; i < (n/2+1); i++)
            datasaver << real(data_out[i])<< " " << imag(data_out[i]) << endl;
        datasaver.close();
    }
    else
    {
        cout << "Failed to open this file! " << endl;
    }

    /* Let's do some clean up lastly */
    fftw_destroy_plan(plan);
    fftw_free(data_in);
    fftw_free(data_out);

    return 0;

}

#endif
