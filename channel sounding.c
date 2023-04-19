#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <fftw3.h>
#include <string.h>

// Function prototypes
void generate_linear_chirp(double f0, double f1, double T, int N, double complex *chirp_signal);
// Generates a linear chirp signal with start frequency f0, end frequency f1, duration T, and N samples.

void channel_sounding(double complex *tx_signal, int N, double complex *channel_impulse_response, int M, double complex *rx_signal);
// Performs channel sounding by convolving the transmitted signal tx_signal with the channel impulse response channel_impulse_response, and storing the result in rx_signal.

void cross_correlation_fft(double complex *x, double complex *y, int N, double complex *result);
// Calculates the cross-correlation of the complex signals x and y using FFTs, and stores the result in result.

double complex *generate_gaussian_noise(int N, double sigma);

// Generates complex Gaussian noise with standard deviation sigma and N samples.

void add_noise(double complex *signal, int N, double complex *noise);
// Adds complex noise to the complex signal signal with N samples.

double calculate_nmse(double complex *actual, double complex *estimated, int N);

// Calculates the Normalized Mean Square Error (NMSE) between the actual and estimated signals, both with N samples.



int main() {
    // Chirp signal parameters
    double f0 = 24.25e9; // Start frequency: 24.25 GHz (5G frequency range)
    double f1 = 27.5e9;  // End frequency: 27.5 GHz (5G frequency range)
    double T = 1e-6;     // Chirp duration: 1 microsecond
    int N = 2048;        // Number of samples

    // Channel impulse response parameters
    int M = 9; // Number of channel taps
    double channel_tap_delays[] = {0.0, 10e-9, 30e-9, 50e-9, 70e-9, 90e-9, 110e-9, 130e-9, 150e-9}; // Tap delays in seconds
    double channel_tap_gains[] = {1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05}; // Tap gains

    // Noise parameters
    double SNR_dB = 20.0; // Signal-to-Noise Ratio in dB

    // Generate the chirp signal
    double complex *chirp_signal = malloc(N * sizeof(double complex));
    generate_linear_chirp(f0, f1, T, N, chirp_signal);

// Generate the channel impulse response
double complex *channel_impulse_response = calloc(N, sizeof(double complex));
for (int i = 0; i < M; i++) {
    int tap_index = (int)(channel_tap_delays[i] * N / T);
    channel_impulse_response[tap_index] = channel_tap_gains[i];
}

// Print the actual impulse response
printf("Actual Impulse Response:\n");
for (int i = 0; i < N; i++) {
    printf("%d: %.6f + j%.6f\n", i, creal(channel_impulse_response[i]), cimag(channel_impulse_response[i]));
}



// Channel sounding: Convolve the chirp signal with the channel impulse response to generate the received signal
double complex *rx_signal = calloc(N + M - 1, sizeof(double complex));
channel_sounding(chirp_signal, N, channel_impulse_response, M, rx_signal);


    // Generate noise and add it to the received signal
    double signal_power = 0.0;
    for (int i = 0; i < N; i++) {
        signal_power += pow(creal(chirp_signal[i]), 2) + pow(cimag(chirp_signal[i]), 2);
    }

    double noise_power = signal_power / pow(10, SNR_dB / 10);
    double noise_sigma = sqrt(noise_power / N);
    double complex *noise = generate_gaussian_noise(N + M - 1, noise_sigma);
    add_noise(rx_signal, N + M - 1, noise);




// Estimate the channel impulse response using cross-correlation between the transmitted and received signals
double complex *estimated_impulse_response = malloc((2 * N - 1) * sizeof(double complex));
cross_correlation_fft(chirp_signal, rx_signal, N, estimated_impulse_response);

// Normalize the estimated impulse response
double normalization_factor = cabs(chirp_signal[0]);
for (int i = 0; i < 2 * N - 1; i++) {
    estimated_impulse_response[i] /= normalization_factor;
}

// Calculate the Normalized Mean Square Error (NMSE) between the actual and estimated impulse responses
double nmse = calculate_nmse(channel_impulse_response, estimated_impulse_response + N - 1, N);
printf("Normalized Mean Square Error (NMSE): %.6f\n", nmse);

// Print the estimated impulse response
printf("Estimated Impulse Response:\n");
for (int i = 0; i < N; i++) {
    printf("%d: %.6f + j%.6f\n", i, creal(estimated_impulse_response[N - 1 + i]), cimag(estimated_impulse_response[N - 1 + i]));
}


// Free allocated memory
free(chirp_signal);
free(channel_impulse_response);
free(rx_signal);
free(noise);
free(estimated_impulse_response);

return 0;
}

void generate_linear_chirp(double f0, double f1, double T, int N, double complex *chirp_signal) {
double k = (f1 - f0) / T;
double delta_t = T / N;

for (int i = 0; i < N; i++) {
    double t = i * delta_t;
    double omega_t = 2.0 * M_PI * (f0 * t + 0.5 * k * t * t);
    chirp_signal[i] = cexp(I * omega_t);
}
}

void channel_sounding(double complex *tx_signal, int N, double complex *channel_impulse_response, int M, double complex *rx_signal) {
for (int i = 0; i < N; i++) {
for (int j = 0; j < M; j++) {
if (i + j < N + M - 1) {
rx_signal[i + j] += tx_signal[i] * channel_impulse_response[j];
}
}
}
}// generates the received signal by convolving the transmitted signal with the channel impulse response.

void cross_correlation_fft(double complex *x, double complex *y, int N, double complex *result) {
int N_fft = 2 * N - 1;
int N_padded = (int)pow(2, ceil(log2(N_fft)));

double complex *x_padded = calloc(N_padded, sizeof(double complex));
double complex *y_padded = calloc(N_padded, sizeof(double complex));

memcpy(x_padded, x, N * sizeof(double complex));
memcpy(y_padded, y, N * sizeof(double complex));

fftw_complex *x_fft = (fftw_complex *)fftw_malloc(N_padded * sizeof(fftw_complex));
fftw_complex *y_fft = (fftw_complex *)fftw_malloc(N_padded * sizeof(fftw_complex));
fftw_complex *result_fft = (fftw_complex *)fftw_malloc(N_padded * sizeof(fftw_complex));

fftw_plan x_plan = fftw_plan_dft_1d(N_padded, x_padded, x_fft, FFTW_FORWARD, FFTW_ESTIMATE);
fftw_plan y_plan = fftw_plan_dft_1d(N_padded, y_padded, y_fft, FFTW_FORWARD, FFTW_ESTIMATE);
fftw_plan result_plan = fftw_plan_dft_1d(N_padded, result_fft, result, FFTW_BACKWARD, FFTW_ESTIMATE);


fftw_execute(x_plan);
fftw_execute(y_plan);

for (int i = 0; i < N_padded; i++) {
    result_fft[i] = x_fft[i] * conj(y_fft[i]);
}

fftw_execute(result_plan);

double scale_factor = 1.0 / N_padded;
for (int i = 0; i < N_padded; i++) {
    result[i] *= scale_factor;
}

fftw_destroy_plan(x_plan);
fftw_destroy_plan(y_plan);
fftw_destroy_plan(result_plan);

fftw_free(x_fft);
fftw_free(y_fft);
fftw_free(result_fft);

free(x_padded);
free(y_padded);
}

double complex *generate_gaussian_noise(int N, double sigma) {
double complex *noise = malloc(N * sizeof(double complex));
srand(time(NULL));

for (int i = 0; i < N; i++) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    double z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
    noise[i] = sigma * (z0 + I * z1);
}

return noise;
}  //calculates the cross-correlation between two signals using the FFT algorithm.

void add_noise(double complex *signal, int N, double complex *noise) {
for (int i = 0; i < N; i++) {
signal[i] += noise[i];
}
} // generates additive white Gaussian noise with the variance.

double calculate_nmse(double complex *actual, double complex *estimated, int N) {
double mse = 0.0;
double signal_power = 0.0;

for (int i = 0; i < N; i++) {
    mse += pow(creal(actual[i] - estimated[i]), 2) + pow(cimag(actual[i] - estimated[i]), 2);
    signal_power += pow(creal(actual[i]), 2) + pow(cimag(actual[i]), 2);
}
//calculates the Normalized Mean Square Error between two signals.

return mse / signal_power;


}
