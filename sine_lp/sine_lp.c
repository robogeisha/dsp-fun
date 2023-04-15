#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846
#define SAMPLE_RATE 48000
#define DURATION 1
#define NUM_SAMPLES (SAMPLE_RATE * DURATION)
#define FREQ 440
#define LOW_PASS_FREQ 2000
#define ALPHA 0.1

void generate_sine_wave(float *buffer, int num_samples, float freq, float sample_rate);
void low_pass_filter(float *buffer, int num_samples, float alpha);

int main() {
    float buffer[NUM_SAMPLES];
    generate_sine_wave(buffer, NUM_SAMPLES, FREQ, SAMPLE_RATE);
    low_pass_filter(buffer, NUM_SAMPLES, ALPHA);

    // Output the filtered samples for visualization or further processing
    for (int i = 0; i < NUM_SAMPLES; i++) {
        printf("%f\n", buffer[i]);
    }

    return 0;
}

void generate_sine_wave(float *buffer, int num_samples, float freq, float sample_rate) {
    float angular_freq = 2 * PI * freq;
    for (int i = 0; i < num_samples; i++) {
        buffer[i] = sinf((float)i * angular_freq / sample_rate);
    }
}

void low_pass_filter(float *buffer, int num_samples, float alpha) {
    float prev_sample = buffer[0];
    for (int i = 1; i < num_samples; i++) {
        buffer[i] = alpha * buffer[i] + (1 - alpha) * prev_sample;
        prev_sample = buffer[i];
    }
}
