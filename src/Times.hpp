#ifndef TIMES_HPP
#define TIMES_HPP

#include <boost/date_time/posix_time/posix_time.hpp>

class Times {
  
public :
  static Times *TIMES;
  static Times *TIMES_UP;

  enum timing_t {simu_time_ = 0, 
		 display_time_,
		 solve_time_,
		 sum_up_time_,
		 total_time_,
		 nTimes};

private :
  boost::posix_time::ptime init_time[nTimes];
  double loop_time[nTimes];
  double time_sum[nTimes];

  //count
  unsigned int frame;
  unsigned int step;

public :
  Times();
  
  void init();

  void tick(unsigned int i);
  void tock(unsigned int i);
  double getTime(unsigned int i);
  double getAverageTime(unsigned int i);
  double getAverageTimeByStep(unsigned int i);
  double getAverageNbStepsByFrame();
  
  void next_loop();
  void next_step();
};

#endif //TIMES_HPP
