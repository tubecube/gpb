#include "utils.hpp"

Logger::Level Logger::_LEVEL = DEBUG;

const char * Logger::LevelMap[4] = {"DEBUG", "INFO", "WARN", "ERROR"};

string Logger::current_time(const char* format)
{
  time_t cur = time(NULL);
  char tmp[200];
  strftime(tmp, sizeof(tmp), format, localtime(&cur));
  return tmp;
}

int makedir(const char* pathname)
{
	return mkdir(pathname, S_IRWXU);
}


void Logger::logging(Level l, const char* format, ...)
{
	if (l >= _LEVEL)
	{
		va_list ap;
		va_start(ap, format);

		FILE *stream = stdout;
		string info = current_time()+" "+LevelMap[l]+": ";
		fprintf(stream, "%s", info.c_str());
		vfprintf(stream, format, ap);
		fflush(stream);

		if (l == ERROR)
		{
			fprintf(stderr, "%s", info.c_str());
			vfprintf(stderr, format, ap);
		}

		va_end(ap);
	}
}
