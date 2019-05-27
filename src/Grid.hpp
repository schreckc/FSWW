#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include "definitions.hpp"
#include <SDL2/SDL_image.h>

class Grid {

private:
  uint n_rows, n_cols, n_nodes;
  FLOAT cell_size;
  std::vector<FLOAT> nodes;

  inline uint index(uint i, uint j) const;
  inline uint row(uint ind) const ;
  inline uint col(uint ind) const;
  
public:
  Grid();
  Grid(uint n_rows, uint n_cols, FLOAT cs);
  ~Grid();
  
  bool isEmpty() const;
  
  uint getNbRows() const;
  uint getNbCols() const;
  FLOAT getCellSize() const;
  void setCellSize(FLOAT cs);
  uint getHeight() const;
  uint getWidth() const;
  VEC2 toWorld(uint i, uint j) const;

  FLOAT operator()(int i, int j) const;
  FLOAT &operator()(int i, int j);
  FLOAT value(FLOAT x, FLOAT y) const;
  FLOAT interpolatedValue(FLOAT x, FLOAT y) const;

  void reset(FLOAT val);
  void setValues(const SDL_Surface *texture);

  void getPixels(unsigned char *pixels) const;
  std::ostream& exportObj(std::ostream& os) const;
  
  Grid& operator=(const Grid& g);
  Grid& operator+=(const Grid& g);
  friend std::ostream& operator<<(std::ostream& os, const Grid& g);
  friend std::istream& operator>>(std::istream& is, Grid& F);
};

#endif
