#include <stdio.h>

// Function to apply the Butterworth filter to input data
void butterworth_filter(double *input, double *output, int n, double cutoff, double sample_rate) {
    double RC = 1.0 / (cutoff * 6.283185307179586);
    double dt = 1.0 / sample_rate;
    double alpha = dt / (RC + dt);
    
    output[0] = input[0];
    for (int i = 1; i < n; i++) {
        output[i] = output[i - 1] + (alpha * (input[i] - output[i - 1]));
    }
}

int main() {
    int n = 10; // Number of samples
    double input_data[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double output_data[n];
    
    double cutoff_frequency = 1; // Cutoff frequency in Hz
    double sample_rate = 10; // Sampling frequency in Hz
    
    butterworth_filter(input_data, output_data, n, cutoff_frequency, sample_rate);
    
    printf("Filtered data: \n");
    for (int i = 0; i < n; i++) {
        printf("%lf\n", output_data[i]);
    }
    
    return 0;
}
