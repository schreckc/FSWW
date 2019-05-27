#ifndef INPUTPOINT_HPP
#define INPUTPOINT_HPP

#include <list>
#include <vector>
#include "definitions.hpp"

class InputPoint {
private:
  VEC2 pos;
  
  uint window_size;
  FLOAT dt;
  FLOAT sample_rate;
  uint nb_frequencies;
  FLOAT frequency_step;
  
  std::list<FLOAT> samples;
  FLOAT *spectrum_re;
  FLOAT *spectrum_im;

  std::list<std::vector<FLOAT> > spectrogram;

  uint t;

  std::string name;
  
public:
  InputPoint();
  InputPoint(uint nb_samples, FLOAT time_step);
  ~InputPoint();
  
  void setPos(VEC2 p);
  void setPos(FLOAT x, FLOAT y);

  VEC2 getPos() const;
  
  void update(FLOAT next_sample);
  void computeSpectrum();

  COMPLEX amplitude(FLOAT frequency) const;

  void plotSpectrum() const;
  void plotSpectrogram() const;
  void plotSamples() const;

  void move(VEC2 trans);
};

#endif
