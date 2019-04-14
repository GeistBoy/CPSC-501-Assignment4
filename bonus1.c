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
double *readWavFile(char *name, int channel);
void writeFile(char *name, double *outputSignalL, double *outputSignalR, int size);
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

double *readWavFile(char *name, int channel)
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
    count = fileLen / 4;

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
    for (int i = channel; i <= count; i+=2)
    {
        int value = ((short int *)buffer)[i];
        input_signal[i] = (double)value / 32768;
    }

    fclose(file);
    free(buffer);
    return input_signal;
}

void writeFile(char *name, double *outputSignalL, double *outputSignalR, int size)
{
    FILE *outputStrem = fopen(name, "wb");
    if (outputStrem == NULL)
    {
        fprintf(stderr, "File %s cannot be opened for writing\n", name);
        return;
    }

    printf("Now writing the output file \"%s\" ...\n", name);
    writeWaveFileHeader(STEREOPHONIC, size, SAMPLE_RATE, outputStrem);

    // Normalization
    double minL = outputSignalL[0];
    double maxL = outputSignalL[0];

    double minR = outputSignalR[0];
    double maxR = outputSignalR[0];

    // find max and min
    for(int i=0; i<size; i++){
        if (outputSignalL[i] < minL)
            minL = outputSignalL[i];
        if (outputSignalL[i] > maxL)
            maxL = outputSignalL[i];

        if (outputSignalR[i] < minR)
            minR = outputSignalR[i];
        if (outputSignalR[i] > maxR)
            maxR = outputSignalR[i];
    }

    short *outputL = (short*)malloc(size * sizeof(short));
    short *outputR = (short *)malloc(size * sizeof(short));
    double denomL = maxL - minL;
    double denomR = maxR - minR;
    for(int i=0; i<size; i++){
        outputSignalL[i] = 2 * (outputSignalL[i] - minL) / denomL + (-1);
        outputL[i] = (short)(outputSignalL[i] * 32767);

        outputSignalR[i] = 2 * (outputSignalR[i] - minR) / denomR + (-1);
        outputR[i] = (short)(outputSignalR[i] * 32767);
    }

    for (int i = 0; i <= size; i++){
        fwriteShortLSB(outputL[i], outputStrem);
        fwriteShortLSB(outputR[i], outputStrem);
    }

    // clean up
    free(outputL);
    free(outputR);
    fclose(outputStrem);
}

// =====================================================================
// Main function

int main(int argc, char *argv[])
{
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

    // Left channel == 0
    // Right channel == 1
    double *irSignalL = readWavFile(argv[2], 0);
    double *irSignalR = readWavFile(argv[2], 1);
    int irSize = count;

    printf("End file reading.\n");
    printf("======================\n");

    // ================================================================
    // Convolve
    printf("Now start convolving...\n");

    double *outputSignalL = (double *)malloc((inputSize + irSize) * sizeof(double) + 1);
    double *outputSignalR = (double *)malloc((inputSize + irSize) * sizeof(double) + 1);

    // Get the nn
    int outputSize = inputSize + irSize - 1;

    int nn = pow(2, ceil(log(outputSize) / LOG2));
    
    int arraySize = 2*nn;
    // build arrays of size of nn * 2
    double *X = (double *)malloc(arraySize * sizeof(double));
    double *HL = (double *)malloc(arraySize * sizeof(double));
    double *HR = (double *)malloc(arraySize * sizeof(double));
    memset(X, 0.0, arraySize);
    memset(HL, 0.0, arraySize);
    memset(HR, 0.0, arraySize);

    // set values for arrays
    for (int i = 0, j = 0; j < inputSize; i += 2, j++)
        X[i] = inputSignal[j];

    for (int i = 0, j = 0; j < irSize; i += 2, j++){
        HL[i] = irSignalL[j];
        HR[i] = irSignalR[j];
    }
        

    // Do the FFT
    four1(X-1, nn, 1);
    four1(HL-1, nn, 1);
    four1(HR-1, nn, 1);

    // Claculate the Y
    double *YL = (double *)malloc(arraySize * sizeof(double));
    double *YR = (double *)malloc(arraySize * sizeof(double));

    for(int i=0; i<arraySize; i+=2){
        YL[i] = (X[i]*HL[i] - X[i+1]*HL[i+1]);
        YL[i+1] = X[i+1]*HL[i] + X[i]*HL[i+1];

        YR[i] = (X[i] * HR[i] - X[i + 1] * HR[i + 1]);
        YR[i + 1] = X[i + 1] * HR[i] + X[i] * HR[i + 1];
    }

    four1(YL-1, nn, -1);
    four1(YR - 1, nn, -1);

    for(int i=0, j=0; i<outputSize; i++, j+=2){
        outputSignalL[i] = YL[j];
        outputSignalR[i] = YR[j];
    }

    printf("Finished convolving.\n");

    // ================================================================
    // Write to the output
    printf("======================\n");
    char *outputFilename = argv[3]; // arg[3]
    writeFile(outputFilename, outputSignalL, outputSignalR, outputSize);
    printf("Finished writing.\n");
    printf("======================\n");

    free(inputSignal);
    free(irSignalL);
    free(irSignalR);
    free(outputSignalL);
    free(outputSignalR);
    free(X);
    free(HL);
    free(HR);
    free(YL);
    free(YR);
    printf("Done!!!!!!!\n");
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
