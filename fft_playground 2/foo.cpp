#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "kiss_fftr.h"

std::vector<float> file_to_samples(const std::string & filename)
{
  std::vector<float> res;
  std::ifstream f {filename};
  if(!f.good()) {
    throw std::runtime_error{"Failed to open " + filename};
  }

  f.seekg(0, f.end);
  res.resize(f.tellg() / 4, 0);
  f.seekg(0, f.beg);
  f.read((char *)&(res[0]), res.size() * 4);
  return res;
}

const int N = 2048, skip = 128;

int main() {
  auto samples = file_to_samples("in.pcm");
  std::vector<float> res;
  res.resize(samples.size());
  for(size_t i = 0; i < res.size(); i++) {
    res[i] = 0.0;
  }

  kiss_fftr_cfg kcfg = kiss_fftr_alloc(N, 0, NULL, NULL);
  kiss_fftr_cfg kicfg = kiss_fftr_alloc(N, 1, NULL, NULL);
  kiss_fft_scalar * time = new kiss_fft_scalar[N];
  kiss_fft_cpx * freq = new kiss_fft_cpx[N/ 2 + 1];
  kiss_fft_cpx * phi = new kiss_fft_cpx[N/ 2 + 1];

  std::vector<float> blackman;
  blackman.resize(N);
  {
    double a = 0.16, a0 = (1.0 - a) / 2.0, a1 = 0.5, a2 = a / 2.0;
    for(size_t i = 0; i < N; i++) {
      blackman[i] = a0
        - a1 * std::cos(2.0 * M_PI * (double)i / ((double)N - 1.0))
        + a2 * std::cos(4.0 * M_PI * (double)i / ((double)N - 1.0));
    }
  }

  size_t len = samples.size(), pos = 0;
  while(pos + N < len) {
    for(size_t i = 0; i < N; i++) {
      time[i] = samples[pos + i] * blackman[i];
    }

    kiss_fftr(kcfg, time, freq);

    for(size_t i = 0; i < N / 2 + 1; i++) {
      auto & rr = freq[i].r;
      auto & ii = freq[i].i;
      phi[i].r = std::sqrt(rr * rr + ii * ii);
      phi[i].i = std::atan2(ii, rr);
    }




    // do cool stuff in each block here


    //k/n*samplerate/2 hetrz
    //r is the scale of the sign wave
    //i is phase offset


    // //adds noise
    // kiss_fft_cpx * tmpBin = new kiss_fft_cpx[N/ 2 + 1];
    // for(size_t i = 0; i < N / 4 + 1; i++) {

    //   if (i<3)
    //   {
    //     tmpBin[i] = phi[i];
    //   }
    //   else {
    //     tmpBin[i] = phi[i*2];
    //     tmpBin[i*2] = phi[i];
    //   }
    //   phi[i].r += ((double)rand() / RAND_MAX) * 0.4 - 0.2;
    //   phi[i].i += ((double)rand() / RAND_MAX) * 5.0 - 2.5;
    //   phi[i] = phi[i*2];
    // }
    // phi = tmpBin;



    // //adds sin boops
    // makes it monstery
    // for(size_t i = 0; i < N / 4 + 1; i++) {


    //   if (i<2) {
    //     continue;
    //   }
    //   phi[i].r -= ((double)rand() / RAND_MAX) * 0.4 - 0.2;
    //   phi[i].i -= ((double)rand() / RAND_MAX) * 5.0 - 2.5;
    //   phi[i] = phi[i*2];
    // }


    // for(size_t i = N/4; i < N / 2 + 1; i++) {
    //   if (i<0) {
    //     continue;
    //   }
    //   //phi[i].r *= phi[i].r;
    //   phi[i].i = (phi[i*2].i);
    //   // phi[i] = phi[i*2];
    // }


    // shift pitch
    int shift = 20;
    kiss_fft_cpx * tmpBin = new kiss_fft_cpx[N/ 2 + 1];
    for(size_t i = 0; i < N / 2 + 1; i++) {

      int j = (i-shift)%N;
      if (i < 3 || i > N/2 + 1 - 3)
      {
        j = i;
      }
      tmpBin[i] = phi[j];
    }
    phi = tmpBin;


    // decrease phase to distort more
    for(size_t i = 0; i < N / 2 + 1; i++) {
      if (i < 3)
      {
        continue;
      }
       //phi[i].i = 1/((double)(i+1));
      phi[i].i *= 1/((double)(i)/5);
      //phi[i].r *= phi[i].r;
    }









    for(size_t i = 0; i < N / 2 + 1; i++) {
      freq[i].r = phi[i].r * std::cos(phi[i].i);
      freq[i].i = phi[i].r * std::sin(phi[i].i);
    }

    kiss_fftri(kicfg, freq, time);

    for(size_t i = 0; i < N; i++) {
      res[pos + i] += time[i];
    }

    pos += skip;
  }

  float vnorm = 0.00001;
  for(size_t i = 0; i < len; i++) {
    vnorm = std::max(vnorm, std::abs(res[i]));
  }

  for(size_t i = 0; i < len; i++) {
    res[i] /= vnorm;
  }

  std::cout.write((char*)(&res[0]), res.size() * 4);
}
