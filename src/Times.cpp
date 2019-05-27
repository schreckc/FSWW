#include "Times.hpp"
#include <time.h>

using namespace boost::posix_time;

Times *Times::TIMES(new Times());
Times *Times::TIMES_UP(new Times());
Times::Times() {
  init();
}

void Times::init() {
  for (unsigned int i = 0; i < nTimes; ++i) {
    init_time[i] = not_a_date_time;
    loop_time[i] = 0;
    time_sum[i] = 0;
  }
  frame = 0;
  step = 0;
}

void Times::tick(unsigned int i) {
  init_time[i] = microsec_clock::local_time();
}

void Times::tock(unsigned int i) {
  assert(init_time[i] != not_a_date_time);
  ptime t_end = microsec_clock::local_time();
  loop_time[i] += (t_end - init_time[i]).total_microseconds()*1e-6;
  init_time[i] = not_a_date_time;
}

double Times::getTime(unsigned int i) {
  return loop_time[i];
}

double Times::getAverageTime(unsigned int i) {
  return time_sum[i]/((double)frame);
}

double Times::getAverageTimeByStep(unsigned int i) {
  return time_sum[i]/((double)step);
}

double Times::getAverageNbStepsByFrame() {
  return (double) step / ((double) frame);
}

void Times::next_loop() {
  for (unsigned int i = 0; i < nTimes; ++i) {
    time_sum[i] += loop_time[i];
    init_time[i] = not_a_date_time;
    loop_time[i] = 0;
  }

  ++frame;
}

void Times::next_step() {
  ++step;
}
