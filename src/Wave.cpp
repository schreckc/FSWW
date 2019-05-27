#include "Wave.hpp"
#include "definitions.hpp"

Wave::Wave() {
  phase_exp = 1;
  phase = 0;
}

bool Wave::isActive() const {
  return is_active;
}

void Wave::setActive(int t) const {
  is_active = true;
}

void Wave::setPhase(FLOAT p) {
  phase_exp = exp(definitions::i_*p);
  phase = p;
}
COMPLEX Wave::getPhaseExp() const {
  return phase_exp;
}

FLOAT Wave::getPhase() const {
  return phase;
}
