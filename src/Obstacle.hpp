#ifndef OBSTACLE_HPP
#define OBSTACLE_HPP

#include <SDL2/SDL_image.h>
#include "InputPoint.hpp"
#include "EquivalentSource.hpp"
#include "MovingEquivalentSource.hpp"
#include "Grid.hpp"
#include "Wave.hpp"
#include <Eigen/SVD>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <fstream>

class WaterSurface;
class ComposedObstacle;

class Obstacle {

protected:
  std::vector<FLOAT> wave_lenghts;
  std::vector<std::vector<EquivalentSource*> >sources_l;
  std::vector<std::vector<InputPoint*> >boundaries_l;
  std::vector<std::vector<VEC2C> >normals_l;
  std::vector<Eigen::BDCSVD<MatrixXcf> >svd_l;

  VEC2 pos;
  uint nb_wl;

  std::vector<int> nb_old_sources;
  bool movable_;

  virtual void setBoundaries(uint w) = 0;
  virtual void setEquivalentSources(uint i) = 0;
  
public:

  //interactive
  void setAmpliSources(const std::list<Wave*> &waves, uint index, int t);
  void setTransfeMatrix(uint index, int t);

  // frequency
  void setAmpliSources(const std::list<Wave*> &waves, uint index);
  void setTransfeMatrix(uint index);

  Obstacle();
  Obstacle(VEC2 p);
  virtual ~Obstacle();

  void reset(const std::vector<std::list<Wave*> > &waves_l,const std::vector<FLOAT> &wl);
  
  void setPos(VEC2 p);
  void setPos(FLOAT x, FLOAT y);
  VEC2 getPos();
  
  void getBoundariesPos(std::list<VEC2> &boundaries) const;

  void update(const WaterSurface *surface);

  FLOAT height(FLOAT x, FLOAT y, FLOAT time) const;
  COMPLEX heightc(FLOAT x, FLOAT y, FLOAT time) const;
  COMPLEX heightc_wl(FLOAT x, FLOAT y, uint wl) const;
  COMPLEX heightc_wl(FLOAT x, FLOAT y, uint w, FLOAT time) const;

  void getSources(std::list<Wave*> &waves, uint w);
  void getSourcesM(std::list<MovingEquivalentSource*> &waves, uint w);
 
  std::vector<InputPoint*> &getBoundaries(uint w);
  std::vector<VEC2C> &getNormales(uint w);
  std::vector<EquivalentSource*> &getSources(uint w);
  std::vector<EquivalentSource*> &getSourcesApprox(uint w);
  int getBoundariesSize(uint w) const;
  int getSourcesSize(uint w) const;
  int getSourcesApproxSize(uint w) const;

  
  virtual FLOAT getGrid(int i, int j) const = 0;
  virtual std::ofstream& exportMitsuba(std::ofstream &file) const;

  void move(VEC2 dir, FLOAT speed, int time);

  void setMovable(bool mv);
  bool isMovable() const;
  

};

#endif
