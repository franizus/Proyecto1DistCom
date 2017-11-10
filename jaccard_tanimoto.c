#include <stdio.h>

struct ElementTuple
{
    char id[12];
    char chemical[150];
};

struct CharacterCount
{
    char character;
    int count;
};

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
