#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SIZE 8
#define PI 3.141592653589793
#define TWO_PI (2.0 * PI)
#define SWAP(a, b) \
    tempr = (a);   \
    (a) = (b);     \
    (b) = tempr

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
const double LOG2 = 0.69314718;

// =====================================================================
// Functions
double *readFile(char *name);
void writeFile(char *name, double *outputSignal, int size);
size_t fwriteShortLSB(short int data, FILE *stream);
size_t fwriteIntLSB(int data, FILE *stream);
void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile);

void four1(double data[], int nn, int isign);

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

void writeFile(char *name, double *outputSignal, int size){
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

    // !!!!!!! Code tunning - Jamming - done
    // !!!!!!! Code tunning - minimize work inside loop - done
    // //!!!!!!! Code tunning - eliminate common subexpression - done
    // find max and min
    for(int i=0; i<size; i++){
        if(outputSignal[i] < min)
            min = outputSignal[i];
        if(outputSignal[i] > max)
            max = outputSignal[i];
    }
    short *output = (short*)malloc(size * sizeof(short));
    double denom = max - min;
    for(int i=0; i<size; i++){
        outputSignal[i] = 2 * (outputSignal[i] - min) / denom + (-1);
        output[i] = (short)(outputSignal[i] * 32767);
    }

    for (int i = 0; i <= size; i++){
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
    clock_t begin = clock();
    //Check the command line args
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

    // Get the nn
    int outputSize = inputSize + irSize - 1;

    clock_t start = clock();
    int nn = pow(2, ceil(log(outputSize) / LOG2)); //!!!!!!! Code tunning - initialize at compile time (change of base rule) - done
                                                //!!!!!!! Code tunning - Use proper data type - done
    

    // build arrays of size of nn * 2
    double *X = (double *)malloc(2 * nn * sizeof(double)); //!!!!!!! Code tunning - use caching
    double *H = (double *)malloc(2 * nn * sizeof(double));
    memset(X, 0.0, 2*nn);
    memset(H, 0.0, 2*nn);

    // set values for arrays
    for (int i = 0, j = 0; j < inputSize; i += 2, j++)
        X[i] = inputSignal[j];

    for (int i = 0, j = 0; j < irSize; i += 2, j++)
        H[i] = irSignal[j];

    // Do the FFT
    four1(X-1, nn, 1);
    four1(H-1, nn, 1);

    // Claculate the Y
    double *Y = (double *)malloc(2 * nn * sizeof(double));

    for(int i=0; i<2*nn; i+=2){
        Y[i] = (X[i]*H[i] - X[i+1]*H[i+1]);
        Y[i+1] = X[i+1]*H[i] + X[i]*H[i+1];
    }

    four1(Y-1, nn, -1);

    clock_t finish = clock();
    double time = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("Cost: %fs\n", time);

    for(int i=0, j=0; i<outputSize; i++, j+=2){
        outputSignal[i] = Y[j];
    }

    printf("Finished convolving.\n");

    // ================================================================
    // Write to the output
    printf("======================\n");
    char *outputFilename = argv[3]; // arg[3]
    writeFile(outputFilename, outputSignal, outputSize);
    printf("Finished writing.\n");
    printf("======================\n");

    free(inputSignal);
    free(irSignal);
    free(outputSignal);
    free(X);
    free(H);
    free(Y);
    printf("Done.\n");
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time: %fs\n", time_spent);
    return 0;
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
    fwriteShortLSB((short) channels, outputFile);

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

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    
    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2)
    {
        if (j > i)
        {
            SWAP(data[j], data[i]);
            SWAP(data[j + 1], data[i + 1]);
        }
        m = nn;
        while (m >= 2 && j > m)
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    mmax = 2;
    while (n > mmax)
    {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2)
        {
            for (i = m; i <= n; i += istep)
            {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j + 1];
                tempi = wr * data[j + 1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j + 1] = data[i + 1] - tempi;
                data[i] += tempr;
                data[i + 1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}
