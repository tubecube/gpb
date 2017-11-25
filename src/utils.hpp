#ifndef UTILS_H
#define UTILS_H

#include <stdarg.h>
#include <time.h>
#include <string>
#include <random>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>

using namespace std;

template<class URNG>
void set_seed(URNG& g)
{
	g.seed(random_device()());
}

class Logger
{
public:
	typedef enum {DEBUG=0, INFO, WARN, ERROR} Level;
	static Level _LEVEL;
	static const char *LevelMap[4];
	static void logging(Level l, const char* format, ...);
	static string current_time(const char* format="%Y-%m-%d %H:%M:%S");
	static int setup_log_dir(const string& dirname);
	static int setup_logfd(const string& filename);
	static FILE* _logfd;
};

#define DEBUG(format, ...) Logger::logging(Logger::DEBUG, format, ## __VA_ARGS__)
#define INFO(format, ...) Logger::logging(Logger::INFO, format, ## __VA_ARGS__)
#define WARN(format, ...) Logger::logging(Logger::WARN, format, ## __VA_ARGS__)
#define ERROR(format, ...) Logger::logging(Logger::ERROR, format, ## __VA_ARGS__)

#endif
