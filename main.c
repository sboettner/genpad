#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <popt.h>
#include <fftw3.h>
#include <sndfile.h>

int arg_numchannels=1;
int arg_numharmonics=64;
int arg_numperiods=16;
float arg_decay=1.0f;
float arg_bandwidth=0.01f;

struct poptOption argtable[]={
    { "mono",       'M', POPT_ARG_VAL,      &arg_numchannels,   1, "output mono sample (one channel)", NULL },
    { "stereo",     'S', POPT_ARG_VAL,      &arg_numchannels,   2, "output stereo sample (two channels)", NULL },
    { "harmonics",  'h', POPT_ARG_INT,      &arg_numharmonics,  0, "number of harmonics to produce", "num" },
    { "periods",    'p', POPT_ARG_INT,      &arg_numperiods,    0, "number of quasi-periods to produce", "num" },
    { "decay",      'd', POPT_ARG_FLOAT,    &arg_decay,         0, "exponent at which the harmonics decay", "exp" },
    { "bandwidth",  'b', POPT_ARG_FLOAT,    &arg_bandwidth,     0, "harmonic bandwidth" },
    POPT_AUTOHELP
    POPT_TABLEEND
};

int numsamples;
double** samples;

int main(int argc, const char* argv[])
{
    poptContext pctx=poptGetContext(NULL, argc, argv, argtable, 0);
    poptSetOtherOptionHelp(pctx, "<output file name>");

    poptGetNextOpt(pctx);

    const char* filename=poptGetArg(pctx);
    if (!filename || poptPeekArg(pctx)) {
        poptPrintUsage(pctx, stderr, 0);
        return 0;
    }

    poptFreeContext(pctx);


    // allocate sample buffers
    numsamples=2*arg_numharmonics*arg_numperiods;
    samples=(double**) malloc(sizeof(double*)*arg_numchannels);
    for (int i=0;i<arg_numchannels;i++) {
        samples[i]=(double*) fftw_malloc(sizeof(double)*numsamples);

        for (int j=0;j<numsamples;j++)
            samples[i][j]=0.0;
    }


    for (int i=1;i<arg_numharmonics;i++) {
        double magnitude=pow((double) i, -arg_decay);
        double bandwidth=arg_bandwidth * i * arg_numperiods;

        int binwidth=(int) ceil(3.5*bandwidth);

        for (int j=0;j<arg_numchannels;j++) {
            for (int k=-binwidth;k<=binwidth;k++) {
                double phase=M_PI*ldexp(rand()&0xffffff, -23);
                double s=(double) k/bandwidth;
                double amp=magnitude * exp(-s*s) / bandwidth;
                int idx=i*arg_numperiods+k;
                if (idx<1 || idx>=numsamples/2) continue;

                samples[j][idx]+=cos(phase) * amp;
                samples[j][numsamples-idx]+=sin(phase) * amp;
            }
        }
    }

    // perform inverse FFT transform
    for (int i=0;i<arg_numchannels;i++) {
        fftw_plan ifft=fftw_plan_r2r_1d(numsamples, samples[i], samples[i], FFTW_HC2R, FFTW_ESTIMATE);
        fftw_execute(ifft);
        fftw_destroy_plan(ifft);
    }


    // normalize sample data
    double maxabs=0.0;
    for (int i=0;i<arg_numchannels;i++)
        for (int j=0;j<numsamples;j++) {
            double v=fabs(samples[i][j]);
            if (v>maxabs)
                maxabs=v;
        }


    // write output file
    SF_INFO sfinfo;
    sfinfo.frames=numsamples;
    sfinfo.samplerate=440 * arg_numharmonics;
    sfinfo.channels=arg_numchannels;
    sfinfo.format=SF_FORMAT_WAV | SF_FORMAT_PCM_32;
    sfinfo.sections=1;
    sfinfo.seekable=0;

    SNDFILE* file=sf_open(filename, SFM_WRITE, &sfinfo);
    if (!file) {
        fprintf(stderr, "Error opening output file '%s'\n", filename);
        return 1;
    }

    double* writebuf=(double*) malloc(numsamples*arg_numchannels*sizeof(double));
    for (int i=0;i<numsamples;i++)
        for (int j=0;j<arg_numchannels;j++)
            writebuf[i*arg_numchannels+j]=samples[j][i] / maxabs;
    
    sf_writef_double(file, writebuf, numsamples);
    free(writebuf);

    sf_close(file);

    return 0;
}
