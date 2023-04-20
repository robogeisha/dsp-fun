#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define N 100 // Number of input data samples
#define M 10  // Length of the adaptive filter

// Function prototypes
double complex *generate_random_complex_data(int n);
void rls_equalizer(double complex *input, double complex *desired, int n, int m, double lambda, double complex *output);
void print_complex_data(const char *label, double complex *data, int n);

int main() {
    // Generate random input and desired data samples
    double complex *input = generate_random_complex_data(N);
    double complex *desired = generate_random_complex_data(N);

    // RLS equalization parameters
    double lambda = 0.99; // Forgetting factor

    // Perform RLS equalization
    double complex *output = malloc(N * sizeof(double complex));
    rls_equalizer(input, desired, N, M, lambda, output);

    // Print input, desired, and output data
    print_complex_data("Input", input, N);
    print_complex_data("Desired", desired, N);
    print_complex_data("Output", output, N);

    // Free allocated memory
    free(input);
    free(desired);
    free(output);

    return 0;
}

double complex *generate_random_complex_data(int n) {
    double complex *data = malloc(n * sizeof(double complex));
    srand(time(NULL));

    for (int i = 0; i < n; i++) {
        double real = (double)rand() / RAND_MAX;
        double imag = (double)rand() / RAND_MAX;
        data[i] = real + I * imag;
    }

    return data;
}

void rls_equalizer(double complex *input, double complex *desired, int n, int m, double lambda, double complex *output) {
    double complex *weights = calloc(m, sizeof(double complex));
    double complex *inversed_R = calloc(m * m, sizeof(double complex));

    // Initialize the inversed_R matrix as identity matrix divided by delta
    double delta = 1e-6;
    for (int i = 0; i < m; i++) {
        inversed_R[i * m + i] = 1.0 / delta;
    }

    for (int k = 0; k < n; k++) {
        // Construct the input vector u
        double complex u[m];
        for (int i = 0; i < m; i++) {
            if (k - i >= 0) {
                u[i] = input[k - i];
            } else {
                u[i] = 0;
            }
        }

        // Calculate k-th output
        double complex y = 0;
        for (int i = 0; i < m; i++) {
            y += weights[i] * u[i];
        }
        output[k] = y;

        // Calculate error
        double complex e = desired[k] - y;

        // Update the inversed_R matrix
        double complex inversed_R_times_u[m];
        for (int i = 0; i < m; i++) {
            inversed_R_times_u[i] = 0;
            for (int j = 0; j < m; j++) {
                inversed_R_times_u[i] += inversed_R[i * m + j] * u[j];
            }
        }

        double complex u_times_inversed_R_times_u = 0;
               for (int i = 0; i < m; i++) {
            u_times_inversed_R_times_u += conj(u[i]) * inversed_R_times_u[i];
        }

        double complex factor = 1.0 / (lambda + u_times_inversed_R_times_u);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                inversed_R[i * m + j] = (inversed_R[i * m + j] - factor * inversed_R_times_u[i] * conj(inversed_R_times_u[j])) / lambda;
            }
        }

        // Update the weights
        for (int i = 0; i < m; i++) {
            double complex delta_w = 0;
            for (int j = 0; j < m; j++) {
                delta_w += inversed_R[i * m + j] * conj(u[j]);
            }
            weights[i] += e * delta_w;
        }
    }

    // Free allocated memory
    free(weights);
    free(inversed_R);
}

void print_complex_data(const char *label, double complex *data, int n) {
    printf("%s:\n", label);
    for (int i = 0; i < n; i++) {
        printf("%f + %fi\n", creal(data[i]), cimag(data[i]));
    }
}


