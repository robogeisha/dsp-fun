#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sndfile.h>

#define NUM_COMB_FILTERS 8

typedef struct {
    float *buffer;
    int bufferSize;
    int writeIndex;
    int readIndex;
    float feedback;
    float lowpassFeedback;
    float lowpassOutput;
} CombFilter;

float processCombFilter(CombFilter *combFilter, float input) {
    float output = combFilter->buffer[combFilter->readIndex];
    float newVal = input + (output * combFilter->feedback);

    combFilter->lowpassOutput = (output * combFilter->lowpassFeedback) + (combFilter->lowpassOutput * (1.0f - combFilter->lowpassFeedback));

    combFilter->buffer[combFilter->writeIndex] = combFilter->lowpassOutput;

    combFilter->writeIndex = (combFilter->writeIndex + 1) % combFilter->bufferSize;
    combFilter->readIndex = (combFilter->readIndex + 1) % combFilter->bufferSize;

    printf("Input: %f, Buffer: %f, Output: %f\n", input, newVal, output);

    return output;
}

float processReverb(CombFilter *combFilters, float input, float gain) {
    float outputSum = 0.0f;

    for (int i = 0; i < NUM_COMB_FILTERS; i++) {
        CombFilter *filter = &combFilters[i];
        
        // Read the value from the buffer
        float output = input + gain * filter->buffer[filter->readIndex];

        // Write the new value to the buffer
        filter->buffer[filter->writeIndex] = output;

        // Update the read and write indices
        filter->readIndex = (filter->readIndex + 1) % filter->bufferSize;
        filter->writeIndex = (filter->writeIndex + 1) % filter->bufferSize;

        // Accumulate the output values from each comb filter
        outputSum += output;
    }

    return outputSum / NUM_COMB_FILTERS;
}


int main() {
    const char *inputFilename = "input.wav";
    const char *outputFilename = "output_reverb.wav";
    float gain = 0.5;


    CombFilter combFilters[NUM_COMB_FILTERS] = {0};

    for (int i = 0; i < NUM_COMB_FILTERS; i++) {
        int delayInSamples = (int)(0.1 * 44100 * (i + 1) / NUM_COMB_FILTERS);
        combFilters[i].bufferSize = delayInSamples;
        combFilters[i].buffer = (float *)calloc(delayInSamples, sizeof(float));
        combFilters[i].writeIndex = 0;
       combFilters[i].readIndex = (combFilters[i].writeIndex + delayInSamples - 1) % combFilters[i].bufferSize;

        combFilters[i].feedback = 0.7;
        combFilters[i].lowpassFeedback = 0.7;
        combFilters[i].lowpassOutput = 0.0;
    }

    SF_INFO fileInfo = {0};
    SNDFILE *inputFile = sf_open(inputFilename, SFM_READ, &fileInfo);

    if (!inputFile) {
        printf("Error: could not open input file '%s'\n", inputFilename);
        return 1;
    }

    SNDFILE *outputFile = sf_open(outputFilename, SFM_WRITE, &fileInfo);

    if (!outputFile) {
        printf("Error: could not open output file '%s'\n", outputFilename);
        sf_close(inputFile);
        return 1;
    }

    int bufferSize = 4096;
    float *inputBuffer = (float *)malloc(bufferSize * sizeof(float));
    float *outputBuffer = (float *)malloc(bufferSize * sizeof(float));

    sf_count_t numFrames;
    int frameCount = 0;
    while ((numFrames = sf_read_float(inputFile, inputBuffer, bufferSize)) > 0) {
        for (int i =
        0; i < numFrames; i++) {
            outputBuffer[i] = processReverb(combFilters, inputBuffer[i], gain);
        }
        sf_write_float(outputFile, outputBuffer, numFrames);
        frameCount += numFrames;
    }
    printf("Processed %d frames\n", frameCount);

    for (int i = 0; i < NUM_COMB_FILTERS; i++) {
        free(combFilters[i].buffer);
    }

    free(inputBuffer);
    free(outputBuffer);

    sf_close(inputFile);
    sf_close(outputFile);

    return 0;
}
