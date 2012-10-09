#ifndef __TIMER_H_
#define __TIMER_H_

#include <string>
#include <sys/time.h>

class Timer
{
public:
	Timer(std::string format);
	~Timer();
private:
	static int depth;
	static std::string prefix_n, prefix_a;
	struct timeval time;
	std::string format;
};

#endif