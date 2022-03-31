# fir-heartbeat-detector
This code is implementation of FIR, LMS and Matched filters for ECG data filtering and heartbeat detection. FIR coefficients were found using inverse fast fourier transform (IFFT) with the help of [fftw](https://www.fftw.org/) library.
## Important files:
- [main.cpp](fir-heartbeat-detector/main.cpp) - contains all classes and functions
- [ecg.dat](fir-heartbeat-detector/ecg.dat) - raw ecg sample data
- [filtered.dat](fir-heartbeat-detector/filtered.dat) - filtered ecg data
- [rpeaks.dat](fir-heartbeat-detector/rpeaks.dat) - filtered data with amplified R-peaks
- [bpms.dat](fir-heartbeat-detector/bpms.dat) - detectected heartrates
# FIR filter
Here FIR filter is implemented to remove DC offset or in other words low frequency noise.
![image](https://user-images.githubusercontent.com/78025384/161136296-30fd070e-52a8-4982-b0a8-8cc128077baa.png)
# LMS filter
LMS filter is an adaptive filter, which adjusts coeffiecients itself. 
It uses 50 HZ noise as a reference which needs to be removed from the signal.
![image](https://user-images.githubusercontent.com/78025384/161136718-b2ba3f60-14fe-42ec-bcd7-aedadbdcc93f.png)
![image](https://user-images.githubusercontent.com/78025384/161136763-3f9edfeb-f127-4de0-870a-fd4e2aec3167.png)
# Matched filter
For heartbeat detection a matched filter is used.
A wavelet is used as coefficients for FIR filter which greatly amplifies R-peaks (heartbeats).
To determine best wavelet shape 3 different wavelets were tested.
![image](https://user-images.githubusercontent.com/78025384/161137295-43fb772f-b6b8-44a2-a493-f4a4435130eb.png)

Mexican hat wavelet showed best filtering performance.
![image](https://user-images.githubusercontent.com/78025384/161137659-4fba5456-16f1-4cc7-8335-657098b66df7.png)
# Heartbeat detection
For heartbeat detection a threshold was applied on final filtered data.
It has been accounted for the fact that a single R-peak contains several samples above the threshold.
Here is the final detected heartrate.
![image](https://user-images.githubusercontent.com/78025384/161139176-ca93fc70-13b9-4fac-ab2f-0221ce6f7f2d.png)
