/* tone_detector.c - C code for Project */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#define PI 3.141592653589793

// Helper function to perform FFT (Cooley-Tukey Radix-2 DIT, power-of-2 only)
typedef struct {
    double real;
    double imag;
} Complex;

void fft(Complex* x, int N) {
    if (N <= 1) return;

    // Divide
    Complex* even = (Complex*)malloc(N / 2 * sizeof(Complex));
    Complex* odd = (Complex*)malloc(N / 2 * sizeof(Complex));
    for (int i = 0; i < N / 2; i++) {
        even[i] = x[i*2];
        odd[i] = x[i*2 + 1];
    }

    // Conquer
    fft(even, N / 2);
    fft(odd, N / 2);

    // Combine
    for (int k = 0; k < N / 2; k++) {
        double t = -2 * PI * k / N;
        Complex twiddle = { cos(t), sin(t) };
        Complex odd_term = {
            twiddle.real * odd[k].real - twiddle.imag * odd[k].imag,
            twiddle.real * odd[k].imag + twiddle.imag * odd[k].real
        };
        x[k].real = even[k].real + odd_term.real;
        x[k].imag = even[k].imag + odd_term.imag;
        x[k + N / 2].real = even[k].real - odd_term.real;
        x[k + N / 2].imag = even[k].imag - odd_term.imag;
    }
    free(even);
    free(odd);
}

void cfar(double sample_rate, double* magnitudes, int N, double thresholdOffset_dB, double* peak_freq) {
    int numGuard = 2;   // Guard cells around CUT
    int numTrain = 10;  // Training cells on each side
    int CUT;
    double noiseSum = 0;
    int numNoiseCells = 0;
    double noiseAvg;
    double threshold;

    double max_peak_value = -INFINITY;
    int peak_bin = -1;

    double freq_low = 500.0;
    double freq_high = 4000.0;

    for (int k = 1; k < N / 2; ++k) {
        CUT = k;
        double freq = (CUT * sample_rate) / N;
        if (freq < freq_low || freq > freq_high)
            continue;  // Skip bins outside the desired frequency range

        noiseSum = 0;
        numNoiseCells = 0;

        for (int i = CUT - numTrain - numGuard; i < CUT - numGuard; ++i) {
            if (i >= 0) {
                noiseSum += magnitudes[i];
                ++numNoiseCells;
            }
        }

        for (int i = CUT + numGuard + 1; i <= CUT + numGuard + numTrain && i < N / 2; ++i) {
            noiseSum += magnitudes[i];
            ++numNoiseCells;
        }

        if (numNoiseCells == 0) continue;

        noiseAvg = noiseSum / numNoiseCells;
        threshold = noiseAvg + thresholdOffset_dB;

        if (magnitudes[CUT] > threshold && magnitudes[CUT] > max_peak_value) {
            max_peak_value = magnitudes[CUT];
            peak_bin = CUT;
        }
    }

    // Output: frequency of peak (or -1 if none found)
    if (peak_bin != -1) {
        *peak_freq = (peak_bin * sample_rate) / N;
    } else {
        *peak_freq = -1.0;  // No peak found
    }
}

	


__declspec(dllexport)
void ProcessAudio(double* input, int length, double sample_rate, double* mag_out, double* peaks, int* tone_detected,double* time_xaxis) {
    // Convert input to complex format
    Complex* x = (Complex*)malloc(length * sizeof(Complex));
	//double* peak_array	= (double*)malloc(length * sizeof(double));
	
    for (int i = 0; i < length; i++) {
        x[i].real = input[i];
        x[i].imag = 0.0;
    }
	for (int i = 0; i < length; i++) {
        time_xaxis[i] = i/sample_rate;
    }

    // Perform FFT
    fft(x, length);

    // Compute magnitude and phase
    for (int i = 0; i < length / 2; i++) {
        mag_out[i] = 20*log10(sqrt(x[i].real * x[i].real + x[i].imag * x[i].imag));
    }
	
	cfar(sample_rate, mag_out, length, 23.0, peaks);
	//peaks = peak_array;
	*tone_detected = (*peaks == -1.0) ? 0 : 1;
	
    free(x);
	//free(peak_array);
}

__declspec(dllexport)
void EstimateDOA(double* mic1, double* mic2, int length, double sample_rate, double mic_distance, double* doa_angle_deg, int* tone_detected, double* tone_frequency,double* out_tone,double* doaplot) {
    Complex* x1 = (Complex*)malloc(length * sizeof(Complex));
    Complex* x2 = (Complex*)malloc(length * sizeof(Complex));
	
    for (int i = 0; i < length; i++) {
        x1[i].real = mic1[i]; x1[i].imag = 0.0;
        x2[i].real = mic2[i]; x2[i].imag = 0.0;
    }

    fft(x1, length);
    fft(x2, length);


    int peak_bin;
    // cfar(sample_rate, mag1, length, 23.0, tone_frequency, &peak_bin);
    // *tone_detected = (*tone_frequency == -1.0) ? 0 : 1;
	peak_bin = (int) ((length / sample_rate) * (*tone_frequency));

    if (*tone_detected && peak_bin > 0) {
        Complex X1 = x1[peak_bin];
        Complex X2 = x2[peak_bin];

        double delta_phase = atan2(X1.imag * X2.real-X1.real * X2.imag ,
                                   X1.real * X2.real + X1.imag * X2.imag);
		
        double lambda = 340.0 / (*tone_frequency);
        double cos_theta = (delta_phase * lambda) / (2 * PI * mic_distance);
		*out_tone = delta_phase;
        if (cos_theta < -1.0) cos_theta = -1.0;
        if (cos_theta > 1.0) cos_theta = 1.0;
        *doa_angle_deg = acos(cos_theta) * (180.0 / PI);
    } else {
        *doa_angle_deg = -1.0; // invalid angle
    }
	
	for (int i = 0; i <  10; i++) {
        doaplot[i] = -i;
    }
    
    free(x1);
    free(x2);
}
