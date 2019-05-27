#ifndef PROJECTED_GRID_HPP
#define PROJECTED_GRID_HPP

#include "settings.hpp"
#include <vector>

class ProjectedGrid {
private:
  uint n_rows, n_cols, n_nodes;
  std::vector<VEC3> viewer_pos;
  std::vector<VEC2> displacement;
  std::vector<FLOAT> sizes;
  FLOAT min_size;

  inline uint index(uint i, uint j) const;
  inline uint row(uint ind) const ;
  inline uint col(uint ind) const;

public:
  ProjectedGrid();
  ProjectedGrid(uint nr, uint nc);
  ~ProjectedGrid();

  VEC3 operator()(uint i, uint j) const;
  //VEC3 &operator()(uint i, uint j);
  VEC3 operator()(uint i) const;
  //VEC3 &operator()(uint i);
  VEC2 getPosWorld(int i, int j) const;
  VEC2 getPosWorld(int i) const;
  
  void setPosOnThePlane(uint i, uint j, FLOAT x, FLOAT y);
  void setHeight(uint i, uint j, FLOAT z);
  void setDisplacement(uint i, uint j, VEC2 d);
  void setDisplacement(uint i, uint j, FLOAT x, FLOAT y);
  
  void setSizes();
  FLOAT minSize();
  FLOAT size(uint i);

  std::ostream& exportObj(std::ostream& os) const;
  
  friend std::ostream& operator<<(std::ostream& os, const ProjectedGrid& g);
  //  friend std::istream& operator>>(std::istream& is, Grid& F);
};

#endif
