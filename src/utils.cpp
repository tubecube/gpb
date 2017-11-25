#include "utils.hpp"

Logger::Level Logger::_LEVEL = DEBUG;

const char * Logger::LevelMap[4] = {"DBG", "INF", "WRN", "ERR"};

FILE* Logger::_logfd = NULL;

string Logger::current_time(const char* format)
{
	time_t cur = time(NULL);
	char tmp[200];
	strftime(tmp, sizeof(tmp), format, localtime(&cur));
	return tmp;
}

void Logger::logging(Level l, const char* format, ...)
{
	if (l >= _LEVEL)
	{
		va_list ap;
		va_start(ap, format);

		string info = current_time()+" "+LevelMap[l]+": ";
		fprintf(stdout, "%s", info.c_str());
		vfprintf(stdout, format, ap);
		fflush(stdout);

		if (_logfd != NULL)
		{
			fprintf(_logfd, "%s", info.c_str());
			va_start(ap, format);
			vfprintf(_logfd, format, ap);
			fflush(_logfd);
		}

		va_end(ap);
	}
}

int Logger::setup_logfd(const string& filename)
{
	if (_logfd) {
		fclose(_logfd);
		_logfd = 0;
	}
	_logfd = fopen(filename.c_str(), "w");
	if (!_logfd)
		return -1;
	return 0;
}

int Logger::setup_log_dir(const string& dirname)
{
	struct stat dirstat;
	if (stat(dirname.c_str(), &dirstat) != 0) {
		if (errno == ENOENT) {
			mkdir(dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			if (stat(dirname.c_str(), &dirstat) != 0) {
				fprintf(stderr, "Warning: could not create dir %s\n", dirname.c_str());
				return -1;
			}
		} else {
			fprintf(stderr, "Warning: could not stat dir %s\n", dirname.c_str());
			return -1;
		}
	}
	return 0;
}
