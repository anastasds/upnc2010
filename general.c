#include "definitions.h"
#include "includes.h"

char * get_input_filename(int argc, char** argv)
{
  char * input_filename;
  if(argc < 2)
    {
      input_filename = (char *) malloc((strlen(DEFAULT_CONFIG_FILE) + 1) * sizeof(char));
      strcpy(input_filename, DEFAULT_CONFIG_FILE);
    }
  else
    input_filename = argv[1];
  return input_filename;
}
