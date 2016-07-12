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
<<<<<<< HEAD
    /* Create the data arrays that will hold the data to be fourier transformed. Need to read out the number of data points
    in the array first */   
    string filename = "fft_data.dat";
    ifstream dataloader(filename.c_str());
    int n = 0;
    if (dataloader.is_open()) 
    {
        dataloader >> n;
        cout << "The no of data points: " << n << endl;
    }
    else
    {
        cout << "Failed to open this file!" << endl;
        return 1;
    }

    /* Since we are doing a real-to-complex fft, the size of the input array is n while that of the output is a (n/2+1)
      array of complex<double> values. */
    double *data_in = (double *)fftw_malloc(n*sizeof(double));    
    /* The size of the output array is n/2 + 1*/
    complex<double> *data_out = (complex<double> *)fftw_malloc((n/2+1)*sizeof(complex<double>));

    /* Before loading in the data, we have to initialize the fft plan. */
    fftw_plan plan = fftw_plan_dft_r2c_1d(n, data_in, reinterpret_cast<fftw_complex*>(data_out), FFTW_MEASURE);

    /* Load in all the data from the data file. 
        We have already read in the first entry which is the number of data points. */
    if (dataloader.is_open())
    {
=======
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
>>>>>>> 5279165ddd213bbebb1a44f483fca852f01b318a
        for (int i =0 ; i < n; i++)
        {
            dataloader >> data_in[i];
        }
        dataloader.close();
    }
    else
    {
        cout << "Failed to open this file!" << endl;
<<<<<<< HEAD
        return 1;
=======
>>>>>>> 5279165ddd213bbebb1a44f483fca852f01b318a
    }

    cout << "The first data point is " << data_in[0] << endl;
    cout << "The last data point is " << data_in[n-1] << endl;
    cout << "Data successfully loaded!!!" << endl;

    cout << "Starting to compute the DFT " << endl;
<<<<<<< HEAD
    /* Now we can compute the actual discrete fourier transform with n samples */
    fftw_execute(plan);

=======
    cout << "Plan completed. Now for the execution!!!!" << endl;
    fftw_execute(plan);
>>>>>>> 5279165ddd213bbebb1a44f483fca852f01b318a
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
