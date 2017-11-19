#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define MAX_ELEMENTS 15000

typedef struct Chemical
{
    char* id;
    char* compound;
} chemical;

//Opens a tsv file and returns a list with the data.
void openFile(chemical *chemicalsList)
{
    FILE *fp;
    char buff[255], temp[50];
    fp = fopen("chemicals.tsv", "r");
    int i = 0;
    while(fgets(buff, sizeof(buff), fp) != NULL)
    {
        sscanf(buff, "%s\t%s\t%s\t%s", temp, chemicalsList[i].id, temp, chemicalsList[i].compound);
        printf("%s\n", chemicalsList[i].id);
        printf("%s\n", chemicalsList[i].compound);
        i++;
    }
    fclose(fp);
}

int main()
{
    chemical *chemicalsList = (chemical*)malloc(sizeof(chemical) * 12422);
    openFile(chemicalsList);
}
