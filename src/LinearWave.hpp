#ifndef LINEARWAVE_HPP
#define LINEARWAVE_HPP

#include "Wave.hpp"

class LinearWave : public Wave {
protected: 
  FLOAT wave_length;
  FLOAT wave_number;
  FLOAT angular_freq;
  VEC2 dir;
  COMPLEX ampli;

public:
  LinearWave();
  LinearWave(FLOAT wl);

  void reset();
  
  void update(FLOAT dt);
  FLOAT height(FLOAT x, FLOAT y, FLOAT time) const;
  FLOAT height(VEC2 p, FLOAT time) const;
  COMPLEX heightc(FLOAT x, FLOAT y, FLOAT time) const;
  COMPLEX heightc(VEC2 p, FLOAT time) const;
  VEC2C gradHeightc(VEC2 p, FLOAT time) const;
  VEC2C gradHeightc(FLOAT x, FLOAT y, FLOAT time) const;
  
  void setDir(VEC2 d);
  void setDir(FLOAT x, FLOAT y);
  void setAmplitude(COMPLEX a);

  FLOAT getDistance(FLOAT x, FLOAT y) const;
  COMPLEX getAmpli(int i) const;
  COMPLEX getAmpli() const;
  uint getIndex() const;
  VEC2 getPos() const;
};


#endif
