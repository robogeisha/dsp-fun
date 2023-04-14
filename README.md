# dsp-fun

comb filter reverb 


This C code implements a simple reverb effect on an input audio file.

First, the code defines a struct called CombFilter which contains information about each individual filter used to create the reverb effect. Each CombFilter struct contains a buffer for storing delayed audio samples, as well as several parameters that control the behavior of the filter, such as feedback and lowpassFeedback.

The processCombFilter function takes an input sample and a CombFilter struct and applies the filter to the input sample, returning the resulting output. The function updates the writeIndex and readIndex of the filter's buffer so that the filter can be applied to the next input sample.

The processReverb function takes an array of CombFilter structs, an input sample, and a gain parameter. It applies each CombFilter to the input sample using the processCombFilter function and then sums the resulting outputs. The function returns the average of the output values across all the CombFilter structs.

In the main function, the code creates an array of CombFilter structs and initializes each one with a different delay time and set of parameters. It then opens an input audio file and an output audio file and reads the audio data from the input file into a buffer. The code then processes each block of audio data in the buffer using the processReverb function and writes the resulting output to the output file.

Finally, the code cleans up by freeing memory allocated for the CombFilter buffers and closing the input and output files.
