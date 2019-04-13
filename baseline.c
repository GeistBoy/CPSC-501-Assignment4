#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE 44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE 16

/*  Standard sample size in bytes  */
#define BYTES_PER_SAMPLE (BITS_PER_SAMPLE / 8)

/*  Number of channels  */
#define MONOPHONIC 1
#define STEREOPHONIC 2

// =====================================================================
// Structure
struct header_file
{
    char chunk_id[4];
    int chunk_size;
    char format[4];
    char subchunk1_id[4];
    int subchunk1_size;
    short int audio_format;
    short int num_channels;
    int sample_rate;
    int byte_rate;
    short int block_align;
    short int bits_per_sample;
    char subchunk2_id[4];
    int subchunk2_size;
};

int count;

// =====================================================================
// Functions
double *readFile(char *name);
void convolve(double *x, int N, double *h, int M, double *y, int P);
void writeFile(char *name, double *outputSignal, int size);
size_t fwriteShortLSB(short int data, FILE *stream);
size_t fwriteIntLSB(int data, FILE *stream);
void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile);

// =====================================================================
// Function implements
double *readFile(char *name)
{
    short int *buffer;
    unsigned long fileLen;
    struct header_file meta;
    //	Open file
    FILE *file = fopen(name, "rb");
    if (!file)
    {
        fprintf(stderr, "Unable to open file %s", name);
        return NULL;
    }
    //	Read header
    fread(&meta, sizeof(meta), 1, file);

    //	Get file length
    fseek(file, 0, SEEK_END);
    fileLen = ftell(file);
    count = fileLen / 2;

    fseek(file, 0, SEEK_SET);

    //	Allocate memory
    buffer = (short int *)malloc(fileLen + 1);
    if (!buffer)
    {
        fprintf(stderr, "Memory error!");
        fclose(file);
        return NULL;
    }

    //	Read file contents into buffer
    fread(buffer, fileLen, 1, file);
    double *input_signal = malloc(count * sizeof(double) + 1);
    for (int i = 0; i <= count; i++)
    {
        int value = ((short int *)buffer)[i];
        input_signal[i] = (double)value / 32768;
    }

    fclose(file);
    free(buffer);
    return input_signal;
}

void writeFile(char *name, double *outputSignal, int size)
{
    FILE *outputStrem = fopen(name, "wb");
    if (outputStrem == NULL)
    {
        fprintf(stderr, "File %s cannot be opened for writing\n", name);
        return;
    }

    printf("Now writing the output file \"%s\" ...\n", name);
    writeWaveFileHeader(MONOPHONIC, size, SAMPLE_RATE, outputStrem);

    // Normalization
    double min = outputSignal[0];
    double max = outputSignal[0];

    // find max
    for (int i = 0; i < size; i++)
    {
        if (outputSignal[i] < min)
            min = outputSignal[i];
    }

    // find min
    for (int i = 0; i < size; i++)
    {
        if (outputSignal[i] > max)
            max = outputSignal[i];
    }

    for (int i = 0; i < size; i++)
        outputSignal[i] = (1 - (-1)) * (outputSignal[i] - min) / (max - min) + (-1); // !!!!!!! Code tunning - Jamming

    short *output = (short *)malloc(size * sizeof(short));
    for (int i = 0; i < size; i++)
        output[i] = (short)(outputSignal[i] * 32767);

    for (int i = 0; i <= size; i++)
    {
        fwriteShortLSB(output[i], outputStrem);
    }

    // clean up
    free(output);
    fclose(outputStrem);
}

// =====================================================================
// Main function

int main(int argc, char *argv[])
{
    // Check the command line args
    if(argc < 4){
        printf("Usage: convolve inputfile IRfile outputfile\n");
        return -1;
    }

    //Double checking input informations
    printf("The inputs are:\n");
    printf("inputfile: %s\n", argv[1]);
    printf("IRfile: %s\n", argv[2]);
    printf("outputfile: %s\n", argv[3]);

    // ================================================================
    // Read Inputs
    printf("======================\n");
    printf("Now starting reading files ...\n");

    double *inputSignal = readFile(argv[1]); 
    int inputSize = count;

    double *irSignal = readFile(argv[2]); 
    int irSize = count;

    printf("End file reading.\n");
    printf("======================\n");

    // ================================================================
    // Convolve
    printf("Now start convolving...\n");

    double *outputSignal = (double *)malloc((inputSize + irSize) * sizeof(double) + 1);

    int outputSize = inputSize + irSize - 1;
    convolve(inputSignal, inputSize, irSignal, irSize, outputSignal, outputSize);

    printf("Finished convolving.\n");

    // ================================================================
    // Write to the output
    printf("======================\n");
    char *outputFilename = argv[3];
    writeFile(outputFilename, outputSignal, outputSize);
    printf("Finished writing.\n");
    printf("======================\n");

    free(inputSignal);
    free(irSignal);
    free(outputSignal);
    printf("Done.\n");
    return 0;
}

void convolve(double *x, int N, double *h, int M, double *y, int P)
{
    int n, m;

    /*  Make sure the output buffer is the right size: P = N + M - 1  */
    if (P != (N + M - 1))
    {
        printf("Output signal vector is the wrong size\n");
        printf("It is %-d, but should be %-d\n", P, (N + M - 1));
        printf("Aborting convolution\n");
        return;
    }

    /*  Clear the output buffer y[] to all zero values  */
    for (n = 0; n < P; n++)
        y[n] = 0.0;

    /*  Do the convolution  */
    /*  Outer loop:  process each input value x[n] in turn  */
    for (n = 0; n < N; n++)
    {
        /*  Inner loop:  process x[n] with each sample of h[]  */
        for (m = 0; m < M; m++)
        {
            y[n + m] += x[n] * h[m];
        }
    }
}

size_t fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}

size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}

void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = channels * numberSamples * BYTES_PER_SAMPLE;

    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;

    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * BYTES_PER_SAMPLE;

    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);

    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);

    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);

    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((short)channels, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}