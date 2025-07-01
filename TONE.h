#ifndef TONE_H
#define TONE_H
 __declspec(dllexport) void ProcessAudio(double* input, int length, double sample_rate, double* mag_out, double* peaks, int* tone_detected,double* time_xaxis);
 __declspec(dllexport) void EstimateDOA(double* mic1, double* mic2, int length, double sample_rate, double mic_distance, double* doa_angle_deg, int* tone_detected, double* tone_frequency,double* out_tone, double* doaplot);
#endif 