#include <iostream>
#include <vector>
#include "solver.h"
using namespace std;

int main(){
    vector<vector<int> > sudoku = {
        {5,4,0,3,8,7,0,0,1},
        {0,0,0,0,0,0,7,0,0},
        {0,3,0,0,1,0,0,0,8},
        {0,0,0,7,2,3,8,6,9},
        {0,0,0,1,0,8,5,7,2},
        {7,0,8,9,5,0,0,1,4},
        {0,0,4,8,7,5,2,3,6},
        {8,0,2,6,3,0,0,9,0},
        {0,6,0,2,0,0,1,8,0}
    };

    solver* sudokuSolver = new solver(sudoku, 9, 3);

    
    return 0;
}