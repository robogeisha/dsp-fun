#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to compute the Butterworth filter coefficients
void butter_coefficients(int N, double Wc, double *a, double *b, int filter_section) {
    double k1 = tan(M_PI * Wc);
    double angle = M_PI * (2 * filter_section + 1) / (2 * N);
    double k2 = 1.0 + 2.0 * k1 * cos(angle) + k1 * k1;

    b[0] = k1 * k1 / k2;
    a[0] = 1.0;
    a[1] = 2.0 * (1.0 - k1 * k1) / k2;
    a[2] = (1.0 - 2.0 * k1 * cos(angle) + k1 * k1) / k2;
}

// Function to apply the IIR filter to an input signal
void iir_filter(double *input, double *output, int n_samples, double *a, double *b, int N) {
    double y;
    for (int i = 0; i < n_samples; ++i) {
        y = b[0] * input[i];
        for (int j = 1; j <= N; ++j) {
            if (i - j >= 0) {
                y += b[j] * input[i - j] - a[j] * output[i - j];
            }
        }
        output[i] = y;
    }
}

int main() {
    int order = 4;
    double cutoff_frequency = 0.1; // Normalized to Nyquist frequency

    // Filter coefficients for each 2nd-order section
    double a1[3], b1[3], a2[3], b2[3];
    butter_coefficients(order, cutoff_frequency, a1, b1, 0);
    butter_coefficients(order, cutoff_frequency, a2, b2, 1);

    // Test input data
    int n_samples = 10;
    double input[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double intermediate_output[10];
    double final_output[10];

    // Apply the two 2nd-order Butterworth filter sections in cascade
    iir_filter(input, intermediate_output, n_samples, a1, b1, 2);
    iir_filter(intermediate_output, final_output, n_samples, a2, b2, 2);

    printf("Filtered output:\n");
    for (int i = 0; i < n_samples; ++i) {
        printf("%.2f\n", final_output[i]);
    }

    return 0;
}
