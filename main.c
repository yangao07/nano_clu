#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debwt_index.h"
#include "nano_clu.h"

#define VERSION "1.0.0"
char PROG[20] = "NanoClu";

static int usage(void)	//main usage
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: %s (Read clustering for Nanopore RNA-seq data)\n\n", PROG);
	fprintf(stderr, "Usage:   %s <command> [options]\n\n", PROG);
	fprintf(stderr, "Command: \n");
	fprintf(stderr, "         index      index gene region sequence\n");
	fprintf(stderr, "         clu        cluster RNA-seq based on each gene region\n");
	fprintf(stderr, "\n");
	return 1;
}

char index_suffix[20] = ".nanoclu";

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "index") == 0)      return debwt_index(argc-1, argv+1);
	else if (strcmp(argv[1], "clu") == 0)   return nano_clu(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
    return 0;
}
