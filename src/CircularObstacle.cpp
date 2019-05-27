#include "CircularObstacle.hpp"
#include "settings.hpp"
#include <iostream>

using namespace settings;

void CircularObstacle::setBoundaries(uint w) {
  FLOAT wl = wave_lenghts[w];
  FLOAT sample_rate_boundary = step_sampling_*wl/2.0;
  FLOAT angle_step = sample_rate_boundary/radius;
  uint nb_sources = 2*M_PI/angle_step +1;
  if (nb_sources < 4) {
    nb_sources = 4;
  }

  angle_step = 2*M_PI/(FLOAT)nb_sources;
  VEC2 dir(1, 0);
  MAT2 rotation;
  rotation <<
    cos(angle_step), -sin(angle_step),
    sin(angle_step), cos(angle_step);
  
  for (uint i = 0; i < nb_sources; ++i) {
    InputPoint* ip = new InputPoint(128, dt_);
    ip->setPos(pos + radius*dir);
    dir = rotation*dir;
    dir.normalize();
    VEC2C n(COMPLEX(dir(0), 0), COMPLEX(dir(1), 0));
    normals_l[w].push_back(n);
    boundaries_l[w].push_back(ip);
  }
}

void CircularObstacle::setEquivalentSources(uint w) {
  FLOAT wl = wave_lenghts[w];
  FLOAT r = radius - offset_;//*wl;
  if (r < 0.1*radius) {
    r = radius*0.1;
  }
  FLOAT sample_rate_boundary = step_sampling_*wl;
  FLOAT angle_step = sample_rate_boundary/r;
  uint nb_sources = 2*M_PI/angle_step +1;
  if (nb_sources < 4) {
    nb_sources = 4;
  }

  angle_step = 2*M_PI/(FLOAT)nb_sources;
  VEC2 dir(1, 0);
  MAT2 rotation;
  rotation <<
    cos(angle_step), -sin(angle_step),
    sin(angle_step), cos(angle_step);

  EquivalentSource* es;

  for (uint i = 0; i < nb_sources; ++i) {
    es = new EquivalentSource(wl, ampli_steps[w]);
    es->setPos(pos + r*dir);
    sources_l[w].push_back(es);
      
    dir = rotation*dir;
  }
}

CircularObstacle::CircularObstacle(FLOAT r): radius(r), Obstacle() {}

CircularObstacle::CircularObstacle(VEC2 p, FLOAT r) : radius(r), Obstacle(p) {}

CircularObstacle::~CircularObstacle() {}

FLOAT CircularObstacle::getGrid(int i, int j) const {
  FLOAT out = 0;
  VEC2 v(i*cell_size_obs - pos(0), j*cell_size_obs - pos(1));
  if (v.norm() < radius) {
    out = 1;
  }
  return out;
}
