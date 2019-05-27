#ifndef CIRCULAROBSTACLE_HPP
#define CIRCULAROBSTACLE_HPP

#include "Obstacle.hpp"

class CircularObstacle: public Obstacle {

protected:
  FLOAT radius;

  void setBoundaries(uint w);
  void setEquivalentSources(uint w);

public:
  CircularObstacle(FLOAT r);
  CircularObstacle(VEC2 p, FLOAT r);
  ~CircularObstacle();

  FLOAT getGrid(int i, int j) const;
};

#endif
