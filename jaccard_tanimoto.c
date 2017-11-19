#include <stdio.h>
#include <stdlib.h>

#define MAX_ELEMENTS 16

typedef struct Chemical
{
    char id[12];
    char compound[150];
} chemical;

//Opens a tsv file and returns a list with the data.
void openFile(chemical *chemicalsList)
{
    FILE *fp;
    char buff[255], temp[10];
    fp = fopen("chemicals.tsv", "r");
    int i = 0;
    while(fgets(buff, sizeof(buff), fp) != NULL)
    {
        sscanf(buff, "%s\t%s\t%s\t%s", temp, chemicalsList[i].id, temp, chemicalsList[i].compound);
        printf("%s\n\n", chemicalsList[i].id);
        printf("%s\n\n", chemicalsList[i].compound);
        i++;
    }
    fclose(fp);
}

int main()
{
    chemical chemicalsList[MAX_ELEMENTS];
    openFile(chemicalsList);
}
