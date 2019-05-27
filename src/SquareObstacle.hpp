#ifndef SQUAREOBSTACLE_HPP
#define SQUAREOBSTACLE_HPP

#include "Obstacle.hpp"

class SquareObstacle: public Obstacle {

protected:
  FLOAT size;

  void setBoundaries(uint w);
  void setEquivalentSources(uint w);

public:
  SquareObstacle(FLOAT s);
  SquareObstacle(VEC2 p, FLOAT s);
  ~SquareObstacle();

  FLOAT getGrid(int i, int j) const;
};

#endif
