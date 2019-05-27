#ifndef WAVYOBSTACLE_HPP
#define WAVYOBSTACLE_HPP

#include "Obstacle.hpp"

class WavyObstacle: public Obstacle {
  uint pattern; // 0: cosine, 1:triangle
  FLOAT freq;
  FLOAT ampli;
  
  void setBoundaries(uint w);
  void setEquivalentSources(uint w);
  
public:
  WavyObstacle(uint p, FLOAT f, FLOAT a);
  ~WavyObstacle();

  FLOAT getGrid(int i, int j) const;

};

#endif
