#include <iostream>
#include <vector>
#include <fftw3.h>
#include <fstream>
#include <cmath>

#define SAMPLES_RATIO 2
#define FS 250
#define LEARNING_RATE 0.001
#define M_PI 3.14159265359
using namespace std;


// Find Inverse of Fast Fourier Transform (IFFT)
void ifft(fftw_complex* in, fftw_complex* out, int N)
{
    // create an IDFT plan
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    // execute the plan
    fftw_execute(plan);
    // do some cleaning
    fftw_destroy_plan(plan);
    fftw_cleanup();
    // scale the output to obtain the exact inverse
    for (int i = 0; i < N; ++i) {
        out[i][0] /= N;
        out[i][1] /= N;
    }
}

/*
Produces FIR filter coefficients based on
fs - sampling rate
fc - cutoff frequency
*/
double* highpassDesign(float fs, float fc) {
    //Number of coefficients/taps
    int m;
    m = SAMPLES_RATIO * fs;
    //Declaring fft complex arrays
    fftw_complex* H;
    H = (fftw_complex*)fftw_malloc(m * sizeof(fftw_complex));
    fftw_complex* x;
    x = (fftw_complex*)fftw_malloc(m * sizeof(fftw_complex));

    //Declare filter gains
    for (int i = 0; i <= int(fc / fs * m); i++) {
        H[i][0] = 0;
        H[i][1] = 0;
    }
    for (int i = int(fc / fs * m) + 1; i < int(m - fc / fs * m); i++) {
        H[i][0] = 1;
        H[i][1] = 0;
    }
    for (int i = int(m - fc / fs * m); i < m; i++) {
        H[i][0] = 0;
        H[i][1] = 0; 
    }

    // Apply inverse fourier transform to get coefficients
    ifft(H, x, m);

    // Exchange halves of coefficients and normalize
    double* out = new double[m];
    for (int i = 0; i < m / 2; i++) {
        out[i] = x[m / 2 + i][0];
    }
    for (int i = m / 2; i < m; i++) {
        out[i] = x[i - m / 2][0];
    }

    return out;
}
/*
Buffer class to store latest values for FIR filter
Commented out parts are unfinished ring buffer implementation
*/
class Buffer {
public:
    int size;
    double* buffer;
    int counter;
    Buffer(int len) {
        this->size = len;
        this->buffer = new double[len];
        for (int i = 0; i < len; i++) {
            this->buffer[i] = 0;
        }
        this->counter = 0;
    };

    void append(double v) {
        //Would be better to use a ring buffer
        for (int i = 0; i < this->size - 1; i++) {
            this->buffer[i] = this->buffer[i + 1];
        }
        this->buffer[this->size - 1] = v;
       // this->buffer[this->counter] = v;
        //cout << this->buffer[this->counter] << endl;
        //this->counter++;
        //if (this->counter >= this->size) { this->counter = 0; }
    }
    
    //double read(int index) {
        //int t = index + counter;
        //if (t > this->size) {
        //    t -= this->size;
        //}
        //return this->buffer[t];
    //}
};

/*
FIR filter class. Stores latest input values in a buffer object
can be implemented as adaptive or non adaptive filter
*/
class FIRfilter {
public:
    FIRfilter(double* coefficients, int len) :_buffer{len} {
        this->_coefficients = coefficients;
    };
    double dofilter(double v) {
        this->_buffer.append(v);
        double filtered=0;
        for (int i = 0; i < _buffer.size; i++) {
            filtered += this->_buffer.buffer[i] * this->_coefficients[i];
        }
        return filtered;
    };
    double doFilterAdaptive(double v, double noise, double rate) {
        double canceller = this->dofilter(noise);
        double output = v - canceller;
        for (int i = 0; i < this->_buffer.size; i++) {
            this->_coefficients[i] = this->_coefficients[i] + output * rate * this->_buffer.buffer[i];
        }
        return output;
    }
private:
    double* _coefficients; // FIR coefficients
    Buffer _buffer;
};


/*
Uses filtered data with distinct r-peaks to determine heartrate
*/
class HeartbeatDetector {

public:
    float fs; // sampling rate
    double threshold; // threshold of an r-peak
    int maxbpm; // maximum bpm used to avoid unreliable data and r peaks of several samples wide
    float bpm;
    HeartbeatDetector(float fs, double threshold = 3*pow(10, -11), int maxbpm = 250) {
        this->fs = fs;
        this->threshold = threshold;
        this->maxbpm = maxbpm;

        this->_counter = 0; 
        this->_cooldown = 1; 
        this->bpm = 0;
        this->_ref_beat = false; 
    };
    /*Detect function works as following:
    1. Checks if recent detection has not been too recent.
    2. If there were no too recent detections, check for a heartbeat.
    3. If a heartbeat is detected check if there is a reference heartbeat (i.e. latest previous beat).
    4. If there is a reference heartbeat calculate the time between them to find bpm.
    */
    float detect(double signal) {
        this->_counter++;
        if (this->_cooldown > 0) {
            this->_cooldown--;
        }
        else  {
            if (signal > this->threshold) {
                if (this->_ref_beat) {
                    this->bpm = (this->_counter / this->fs) * 60.0;
                }
                else {
                    this->_ref_beat = true;
                }
                this->_counter = 0;
                this->_cooldown = int((60.0 / this->maxbpm) * this->fs);
            }
        }
        return this->bpm;
    }
private:
    int _counter; // used to determine bpm
    int _cooldown; // checks if the bpm is not over maximum
    bool _ref_beat; // used to check if there was a reference beat already (false - no reference)
};



int main()
{
    //Loading Example Data
    ifstream file("ecg.dat");
    float num;
    vector<float> data;
    while (file >> num) {
        data.push_back(num);
    }
    
    //Initialising Highpass Filter
    FIRfilter highpass_filter(highpassDesign(FS, 0.7), FS*SAMPLES_RATIO);
    
    // Initialising LMS filter with empty untrained coefficients
    double empty_arr[FS * SAMPLES_RATIO] = { };
    FIRfilter adaptive_filter(empty_arr, FS*SAMPLES_RATIO);

    ofstream file_filtered;
    file_filtered.open("filtered.dat", ios::out);
    ofstream file_rpeaks;
    file_rpeaks.open("rpeaks.dat", ios::out);
    ofstream file_bpms;
    file_bpms.open("bpms.dat", ios::out);

    //Mexican hat wavelet matched filter construction for r-peak detection
    double x_mhat[50];
    double wavelet_mhat[50];
    for (int i = 0; i < 50; i++) {
        x_mhat[i] = (i - 50/2) / 4.5;
        wavelet_mhat[49-i] = 0.0013 * (1 - pow(x_mhat[i], 2)) * exp(-(pow(x_mhat[i], 2)) / 2);
    }
    
    FIRfilter mhat_filter(wavelet_mhat, 50);

    HeartbeatDetector detector(FS);
    
    // Signal Sampling Simulation
    double no_DC, noise, lms, rp;
    float bpm;
    for (int i = 0; i < size(data); i++) {
        no_DC = highpass_filter.dofilter(data[i]);
        noise = sin(2 * M_PI * 50/FS * i);
        lms = adaptive_filter.doFilterAdaptive(no_DC, noise, LEARNING_RATE);
        file_filtered << lms << endl;
        rp = mhat_filter.dofilter(lms);
        rp = pow(rp, 2);
        file_rpeaks << rp << endl;
        // Heartbeat detection
        bpm = detector.detect(rp);
        file_bpms << bpm << endl;
    }
    
    file_filtered.close();
    file_rpeaks.close();
    file_bpms.close();

    return 0;
}
