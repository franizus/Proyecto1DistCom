#include <stdio.h>

void openFile(struct ElementTuple *tuple);
{
    FILE *fptr;
    fptr = fopen("chemicals.tsv","r");
    fscanf(fptr,"%d", &num);
    fclose(fptr); 
}

int main()
{
    struct ElementTuple tuple[12000];
    return 0;
}
