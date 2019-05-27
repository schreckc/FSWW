#include "ui_parameters.hpp"

namespace ui_parameters {
  bool show_in_field = true;
  bool show_scattered_field = true;

  bool circle_ = false;
  bool square_ = false;
  bool island_ = false;
  bool harbour_ = false;
  bool line_ = false;
  bool test_ = false;
  bool wavy_ = false;

  bool point_source1_ = false;
  bool point_source2_ = false;
  bool point_source3_ = false;
  bool linear_wave1_ = false;
  bool linear_wave2_ = false;
  bool linear_wave3_ = false;
  bool user_def_source_ = false;

  VEC2 user_def_pos_(0, 0);
  
  int dipole_ = 0;

  bool running_ = false;

  FLOAT size_drop_ = 10;
}
