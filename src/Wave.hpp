#ifndef WAVE_HPP
#define WAVE_HPP

#include "definitions.hpp"

class Wave {
protected:
  mutable bool is_active;
  //  mutable bool changed_active;
  COMPLEX phase_exp;
  FLOAT phase;

public:
  Wave();
  
  virtual void update(FLOAT dt) = 0;
  virtual FLOAT height(FLOAT x, FLOAT y, FLOAT time) const = 0;
  virtual FLOAT height(VEC2 p, FLOAT time) const = 0;
  virtual COMPLEX heightc(FLOAT x, FLOAT y, FLOAT time) const = 0;
  virtual COMPLEX heightc(VEC2 p, FLOAT time) const = 0;
  virtual VEC2C gradHeightc(FLOAT x, FLOAT y, FLOAT time) const = 0;
  virtual VEC2C gradHeightc(VEC2 p, FLOAT time) const = 0;

  virtual FLOAT getDistance(FLOAT x, FLOAT y) const = 0;
  virtual COMPLEX getAmpli(int i) const = 0;
  virtual COMPLEX getAmpli() const = 0;
  virtual uint getIndex() const = 0;
  virtual VEC2 getPos() const = 0;

  bool isActive() const;
  virtual void setActive(int t) const;

  void setPhase(FLOAT p);
  COMPLEX getPhaseExp() const;
  FLOAT getPhase() const;
};


#endif
