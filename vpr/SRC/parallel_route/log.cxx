#include "log.h"

zlog_category_t *delta_log;

using namespace std;

void init_logging(int num_threads)
{
 	if (dzlog_init("log.conf", "default") == -1) {
		printf("failed to init zlog\n");
		return;
	}

	delta_log = zlog_get_category("delta");
}
