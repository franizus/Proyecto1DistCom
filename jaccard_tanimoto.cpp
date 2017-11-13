#include <vector>
#include <string>
#include <fstream>
#include <tuple>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>
#include <chrono>
#include <math.h>
#include <omp.h>

std::vector<std::tuple<std::string, std::string>> openFile()
{
    std::string chemicalString;
    std::ifstream infile("chemicals.tsv");
    std::vector<std::tuple<std::string, std::string>> list;

    while(std::getline(infile, chemicalString))
    {
        std::istringstream iss(chemicalString);
        std::vector<std::string> chemical;
        chemical.resize(4);
        (iss >> chemical[0] >> chemical[1] >> chemical[2] >> chemical[3]);
        list.push_back(std::make_tuple(chemical.at(1), chemical.at(3)));
    }

    infile.close();

    return list;
}

std::map<char, int> analyzeString(std::string chemicalCompound)
{
    std::map<char, int> analyzedString;
    for(char& c : chemicalCompound)
    {
        if (analyzedString.count(c) > 0)
        {
            analyzedString[c]++;
        }
        else
        {
            analyzedString[c] = 1;
        }
    }
    if (analyzedString.count('@') > 0)
    {
        analyzedString['@'] = 1;
    }
    return analyzedString;
}

int getNumberCommonElements(std::map<char, int> chemicalA, std::map<char, int> chemicalB)
{
    int numberElements = 0;
    for(auto& key : chemicalA)
    {
        if (chemicalB.count(key.first) > 0)
        {
            numberElements += std::min(key.second, chemicalB[key.first]);
        }
    }
    return numberElements;
}

int getNumberChemicalElements(std::map<char, int> chemicalCompound)
{
    int numberElements = 0;
    for(auto& key : chemicalCompound)
    {
        numberElements += key.second;
    }
    return numberElements;
}

float getJacTanCoefficient(int elementsA, int elementsB, int commonElements)
{
    return roundf(((float) commonElements / (elementsA + elementsB - commonElements)) * 100) / 100;
}

std::vector<int> getPivots(int chemicalsLength, int numberProcessors)
{
    std::vector<int> pivotsList;    
    pivotsList.push_back(0);
    for (int i = numberProcessors - 1; i >= 1; i--)
    {
        pivotsList.push_back(int(roundf(chemicalsLength - (sqrt((float) i / numberProcessors) * chemicalsLength))));
    }
    pivotsList.push_back(chemicalsLength - 1);
    return pivotsList;
}

void fillComparedList(std::vector<std::tuple<std::string, std::string>> chemicalsList, int pivotMin, 
    int pivotMax, std::vector<std::tuple<std::string, std::string, float>> &comparedChemicalsList)
{
    std::vector<std::tuple<std::string, std::string, float>> row;
    std::map<char, int> lettersA, lettersB;
    float coefficient;
    for (int i = pivotMin; i < pivotMax; i++)
    {
        for (int j = i + 1; j < chemicalsList.size(); j++)
        {
            lettersA = analyzeString(std::get<1>(chemicalsList.at(i)));
            lettersB = analyzeString(std::get<1>(chemicalsList.at(j)));
            coefficient = getJacTanCoefficient(getNumberChemicalElements(lettersA), 
                                               getNumberChemicalElements(lettersB), 
                                               getNumberCommonElements(lettersA, lettersB));
            comparedChemicalsList.push_back(std::make_tuple(std::get<0>(chemicalsList.at(i)),
                                            std::get<0>(chemicalsList.at(j)),
                                            coefficient));
        }
    }
}

void writeFile(std::vector<std::vector<std::tuple<std::string, std::string, float>>> comparedChemicalsList, 
    int numberThreads, std::chrono::duration<double> elapsedSeconds)
{
    std::ofstream ofs("chem_sim_total_C++.tsv", std::ofstream::out);
    ofs << "Chem_ID_1\tChem_ID_2\tTanimoto_similarity\n";
    for (int i = 0; i < numberThreads; i++)
    {
        for (int j = 0; j < comparedChemicalsList[i].size(); j++)
        {
             ofs << std::get<0>(comparedChemicalsList.at(i).at(j)) << "\t" << std::get<1>(comparedChemicalsList.at(i).at(j)) << "\t" << std::get<2>(comparedChemicalsList.at(i).at(j)) << "\n";
        }
    }
    ofs << "Total time = " << elapsedSeconds.count() << "\n";
}

int main()
{
    int numberThreads = omp_get_max_threads();
    std::vector<std::tuple<std::string, std::string>> chemicalsList = openFile();
    std::vector<int> pivots = getPivots(chemicalsList.size(), numberThreads);
    std::vector<std::vector<std::tuple<std::string, std::string, float>>> comparedChemicals;
    comparedChemicals.resize(8);
	auto start = std::chrono::system_clock::now();
    #pragma omp parallel num_threads(numberThreads)
    {
        int pid = omp_get_thread_num();
        int pivotMin = pivots[pid];
        int pivotMax = pivots[pid + 1];
        std::vector<std::tuple<std::string, std::string, float>> comparedChemicalsList;
        fillComparedList(chemicalsList, pivotMin, pivotMax, comparedChemicalsList);
        comparedChemicals[pid] = comparedChemicalsList;
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    
    writeFile(comparedChemicals, numberThreads, elapsedSeconds);

    return 0;
}
