#ifndef SOLVER_H
#define SOLVER_H
#include <iostream>
#include <vector>
#include <set>
#include "omp.h"
#include <stack> 
#include <algorithm>
#include <chrono>
using namespace std;
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;
class solver{
    public:
        solver(vector<vector<int>> &sudokuBoard, int N, int miniGrid){
            high_resolution_clock::time_point begin = high_resolution_clock::now();
            vector<vector<int>> tmpSudoku = sudokuBoard; 
            this->sudokuBoard = sudokuBoard;
            this->N = N;
            this->miniGrid = miniGrid;
            possibleValuesForSudoku.resize(N);
            for(auto &i: possibleValuesForSudoku){
                i.resize(N);
            }
            updateAllPossibleValues();
            cout << "SUDOKU BOARD BEFORE:" << endl;
            printCurrentSudokuBoard();
            while(eliminationSearch() || loneRangerSearch());
            bruteForceDFS();
            auto time_span = duration_cast<duration<double>>(high_resolution_clock::now() - begin);
            if(solved)
            {
                cout << "FOUND SOLUTION!" << endl;
                cout << "DURATION WITH PARALLELISM: " << time_span.count() << endl; 
                high_resolution_clock::time_point begin = high_resolution_clock::now();
                noParallelism(tmpSudoku);
                auto time_span = duration_cast<duration<double>>(high_resolution_clock::now() - begin);
                cout << "DURATION WITH NO PARALLELISM: " << time_span.count() << endl; 
                printCurrentSudokuBoard();
            }
            else
            {
                cout << "NO SOLUTION!" << endl;
            }




        }
        bool noParallelism(vector<vector<int>> board, int row = 0, int column = 0)
        {
            if(row == N)
            {
                return true;
            }
            else if(column == N)
            {
                return solveRemainingBoardUsingBacktracking(board, row + 1, 0);
            }
            else if(board[row][column] != 0)
            {
                return solveRemainingBoardUsingBacktracking(board, row, column + 1);
            }
            else
            {
                for(int i = 1; i <= N; ++i)
                {
                    if(validColumn(board, row, column, i) && validRow(board, row, column, i) && validSubgrid(board, row, column, i))
                    {
                        board[row][column] = i;
                        if(solveRemainingBoardUsingBacktracking(board, row, column + 1))
                        {
                            return true;
                        }
                        board[row][column] = 0;
                    }
                }
                return false;
            }
        }   
    private:
        bool solved = false;
        vector<vector<set<int>>> possibleValuesForSudoku; //has all possible values for solver using set
        vector<vector<int>> sudokuBoard; //pointer to actual board being passed in
        int N; //dimensions of grid
        int miniGrid; //dimensions for subgrid, for 9x9 = 3
        void bruteForceDFS()
        {
            stack<vector<int>> combinations;
            pair<int,int> indexesForFirstSevenCells[7];
            int cnt = 0; 
            for(int i = 0; i < N; ++i){
                for(int j = 0; j < N; ++j){
                    if(*possibleValuesForSudoku[i][j].begin() == -1)
                    {
                        continue;
                    }
                    indexesForFirstSevenCells[cnt] = make_pair(i,j);
                    cnt++;
                    if(cnt == 7)
                    {
                        break;
                    }
                }
                if(cnt == 7)
                {
                    break;
                }
            }
            vector<vector<int>> generateAllCombinations;
            
            for(int i = 0; i < 7; ++i)
            {
                int row = indexesForFirstSevenCells[i].first;
                int column = indexesForFirstSevenCells[i].second;
                vector<int> convertSet( possibleValuesForSudoku[row][column].begin(), possibleValuesForSudoku[row][column].end() );
                generateAllCombinations.push_back(convertSet);
            }
            vector<vector<int>> resultingCombinations;
            vector<int> current;
            recursiveGenerateCombinations(generateAllCombinations, current, resultingCombinations, 0);
            for(int i  = 0; i < resultingCombinations.size(); ++i)
            {
                combinations.push(resultingCombinations[i]);
            }
            omp_lock_t stackLock, foundLock;
            omp_init_lock(&stackLock);
            omp_init_lock(&foundLock);

            #pragma omp parallel 
            {
                while(combinations.size() != 0)
                {
                    vector<vector<int>> resultingVector = sudokuBoard;
                    vector<int> popFromStack;
                    omp_set_lock(&stackLock);
                    if(!combinations.empty())
                    {
                        popFromStack = combinations.top();
                        combinations.pop();
                    }
                    else
                    {
                        omp_unset_lock(&stackLock);
                        break;
                    }
                    omp_unset_lock(&stackLock);

                    if(isSolution(resultingVector, popFromStack))
                    {
                        omp_set_lock(&foundLock);
                        solved = true;
                        sudokuBoard = resultingVector; 
                        omp_unset_lock(&foundLock);
                        break;
                    }
                }

            }
            omp_destroy_lock(&stackLock);
            omp_destroy_lock(&foundLock);
        }
        bool isSolution(vector<vector<int>> &resultingVector, vector<int> popFromStack)
        {
            pair<int, int> indexesForFirstSevenCells[7];
            int cnt = 0;
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    if (*possibleValuesForSudoku[i][j].begin() != -1) {
                        if (cnt == 7){
                            break;
                        }
                        indexesForFirstSevenCells[cnt] = make_pair(i, j);
                        cnt++;
                    }
                }
                if (cnt == 7)
                {
                    break;
                }

            }
            for (int i = 0; i < 7; ++i) {
                int row = indexesForFirstSevenCells[i].first;
                int column = indexesForFirstSevenCells[i].second;
                if(popFromStack[i] == -1)
                {
                    continue; 
                }
                resultingVector[row][column] = popFromStack[i];
                if(!validColumn(resultingVector,row, column ))
                {
                    return false;
                }
                if(!validRow(resultingVector, row, column))
                {
                    return false;
                }
                if(!validSubgrid(resultingVector, row, column))
                {
                    return false;
                }
            }
            if(solveRemainingBoardUsingBacktracking(resultingVector))
            {
                return true; 
            }
            else
            {
                return false; 
            }

        }
        bool solveRemainingBoardUsingBacktracking(vector<vector<int>> &board, int row = 0, int column = 0)
        {
            if(row == N)
            {
                return true;
            }
            else if(column == N)
            {
                return solveRemainingBoardUsingBacktracking(board, row + 1, 0);
            }
            else if(board[row][column] != 0)
            {
                return solveRemainingBoardUsingBacktracking(board, row, column + 1);
            }
            else
            {
                for(auto it = possibleValuesForSudoku[row][column].begin(); it != possibleValuesForSudoku[row][column].end(); ++it)
                {
                    if(validColumn(board, row, column, *it) && validRow(board, row, column, *it) && validSubgrid(board, row, column, *it))
                    {
                        board[row][column] = *it;
                        if(solveRemainingBoardUsingBacktracking(board, row, column + 1))
                        {
                            return true;
                        }
                        board[row][column] = 0;
                    }
                }
                return false;
            }
        }
        bool validColumn(vector<vector<int>> resultingVector, int row, int column, int k = 0)
        {
            if(k)
            {
                resultingVector[row][column] = k; 
            }
            set<int> seen;
            for (int i = 0; i < N; ++i) {
                if(resultingVector[i][column] == 0)
                {
                    continue;
                }
                if(seen.count(resultingVector[i][column]))
                 {
                    return false;
                }
                seen.insert(resultingVector[i][column]);
            }
            return true;
        }
        bool validRow(vector<vector<int>> resultingVector, int row, int column, int k = 0){
            if(k)
            {
                resultingVector[row][column] = k; 
            }
            set<int> seen;
            for (int j = 0; j < N; ++j) {
                if(resultingVector[row][j] == 0)
                {
                    continue;
                }
                if(seen.count(resultingVector[row][j]))
                {
                    return false;
                }
                seen.insert(resultingVector[row][j]);

            }

            return true;
        }
        bool validSubgrid(vector<vector<int>> resultingVector, int row, int column, int k = 0)
        {
            if(k)
            {
                resultingVector[row][column] = k; 
            }
            set<int> seen;
            int startRow = (row/miniGrid) * miniGrid;
            int startCol = (column/miniGrid) * miniGrid;
            for (int i = startRow; i < startRow + miniGrid; ++i) {
                for (int j = startCol; j < startCol + miniGrid; ++j) {
                    if(resultingVector[i][j] == 0)
                    {
                        continue;
                    }
                    if (seen.count(resultingVector[i][j])) 
                    {
                        return false;
                    }
                    seen.insert(resultingVector[i][j]);
                    
                }
            }

            return true;
        }
        void recursiveGenerateCombinations(vector<vector<int>>& generateAllCombinations, vector<int>& current, vector<vector<int>>& resultingCombinations, int index){
            if(index == generateAllCombinations.size()){
                resultingCombinations.push_back(current);
                return;
            }
            for(int i = 0; i < generateAllCombinations[index].size(); ++i)
            {
                current.push_back(generateAllCombinations[index][i]);
                recursiveGenerateCombinations(generateAllCombinations, current, resultingCombinations, index + 1);
                current.pop_back();
            }
            
        }
        void processRows(int j, set<int> &possibleValues){
            set<int> rowValues;
            for(int rowItr = 0; rowItr < N; ++rowItr){
                if(sudokuBoard[rowItr][j] != 0){
                    rowValues.insert(sudokuBoard[rowItr][j]);
                }
            }
            #pragma omp critical
            {
                for(int value: rowValues){
                    possibleValues.erase(value);
                }
            }
        }
        void processColumns(int i, set<int> &possibleValues)
        {
            set<int> colValues;
            for(int colItr = 0; colItr < N; ++colItr){
                if(sudokuBoard[i][colItr] != 0){
                    colValues.insert(sudokuBoard[i][colItr]);
                }
            }
            #pragma omp critical
            {
                for (int value : colValues) {
                    possibleValues.erase(value); 
                }
            }
        }
        void processSubgrid(int i, int j, set<int> &possibleValues)
        {
            set<int> subGridValues;
            //for index (4,4) 9x9 = (4/3) * 3 = 3
            int startRow = (i/miniGrid) * miniGrid;
            int startCol = (j/miniGrid) * miniGrid;
            for (int row = startRow; row < startRow + miniGrid; ++row) {
                for (int col = startCol; col < startCol + miniGrid; ++col) {
                    if (sudokuBoard[row][col] != 0) {
                        subGridValues.insert(sudokuBoard[row][col]);
                    }
                }
            }
            #pragma omp critical
            {
                for (int value : subGridValues) {
                    possibleValues.erase(value); 
                }
            }
        }
        void updateAllPossibleValues(){
            for(int i = 0; i < N; ++i){
                for(int j = 0; j < N; ++j){
                    if(sudokuBoard[i][j] != 0){
                        possibleValuesForSudoku[i][j] = {-1};
                        continue;
                    }
                    set<int> possibleValues;
                    for(int values = 1; values <= N; ++values){
                        possibleValues.insert(values);
                    }
                    #pragma omp parallel sections 
                    {
                        //process rows
                        #pragma omp section
                        {
                            processRows(j, possibleValues);
                        }
                        
                        //process columns
                        #pragma omp section
                        {
                            processColumns(i, possibleValues);
                        }
                    
                        // process subgrid
                        #pragma omp section
                        {
                            processSubgrid(i,j, possibleValues);
                        }
                    }

                    possibleValuesForSudoku[i][j] = possibleValues;

                }
            }
        }
        void printAllPossibleValuesFunction(){
            for(int i = 0; i < N; ++i){
                cout << "[";
                for(int j = 0; j < N; ++j){
                    cout << "{";
                    for(auto it = possibleValuesForSudoku[i][j].begin(); it != possibleValuesForSudoku[i][j].end(); ++it){
                        cout << *it; 
                    }
                    cout << "}, ";
                }
                cout << "]" << endl;
            }
            cout << endl;
        }
        void printCurrentSudokuBoard()
        {
            for(int i = 0; i < N; ++i){
                for(int j = 0; j < N; ++j){
                    cout << sudokuBoard[i][j] << " ";
                }
                cout << endl;

            }
        }
        bool eliminationSearch(){
            //first check using elimation, checking whole grid to see if there is one value in the set using parallelism
            bool hasModified = false;
            #pragma omp parallel 
            {
                set<int> findSoloValues;
                int splitThreads = N/omp_get_num_threads();
                int id = omp_get_thread_num();
                int startingRow = id * splitThreads;
                int endingRow = id == omp_get_num_threads() - 1? N-1: startingRow + splitThreads - 1; //has to be inclusive
                
                for(int i = startingRow; i <= endingRow; ++i){
                    for(int j = 0; j < N; ++j){
                        if(possibleValuesForSudoku[i][j].size() == 1 && *possibleValuesForSudoku[i][j].begin() != -1){
                            int possibleValues = *possibleValuesForSudoku[i][j].begin();
                            #pragma omp critical
                            {
                                hasModified = true;
                                sudokuBoard[i][j] = possibleValues;
                                possibleValuesForSudoku[i][j] = {-1};
                            }
                        }

                    }
                }
            }
            if(hasModified)
            {
                updateAllPossibleValues();
            }
            return hasModified; 
        }
        bool loneRangerSearch()
        {
            bool hasModified = false;
            for(int i = 0; i < N; ++i){
                for(int j = 0; j < N; ++j){
                    vector<int> convertSet(possibleValuesForSudoku[i][j].begin(), possibleValuesForSudoku[i][j].end());
                    for(int element = 0; element < convertSet.size(); ++element)
                    {
                        int valueToFind = convertSet[element]; 
                        bool rowIsLoneRanger = true;
                        bool colIsLoneRanger = true;
                        bool subgridIsLoneRanger = true; 

                        #pragma omp parallel sections
                        {
                            //checking for lone rangers in row
                            #pragma omp section
                            {
                                for(int rows = 0; rows < N; ++rows){
                                    if(rows == i)
                                    {
                                        continue;
                                    }
                                    if(possibleValuesForSudoku[rows][j].contains(valueToFind) || sudokuBoard[rows][j] == valueToFind){
                                        #pragma omp critical
                                        {
                                            rowIsLoneRanger = false;
                                        }
                                    }
                                }
                            }
                            //checking for lone rangers in column
                            #pragma omp section
                            {
                                for(int column = 0; column < N; ++column){
                                    if(column == j)
                                    {
                                        continue;
                                    }
                                    if(possibleValuesForSudoku[i][column].contains(valueToFind)|| sudokuBoard[i][column] == valueToFind){
                                        #pragma omp critical
                                        {
                                            colIsLoneRanger = false;
                                        }
                                    }
                                }
                            }
                            //checking for lone rangers in subgrid
                            #pragma omp section 
                            {
                                set<int> subGridValues;

                                //for index (4,4) 9x9 = (4/3) * 3 = 3
                                int startRow = (i/miniGrid) * miniGrid;
                                int startCol = (j/miniGrid) * miniGrid;
                                for(int row = startRow; row < startRow+miniGrid; ++row)
                                {
                                    for(int column = startCol; column < startCol + miniGrid; ++column )
                                    {
                                        if(row == i && column == j)
                                        {
                                            continue;
                                        }
                                        if(possibleValuesForSudoku[row][column].contains(valueToFind) || sudokuBoard[row][column] == valueToFind)
                                        {
                                            #pragma omp critical
                                            { 
                                                subgridIsLoneRanger = false;
                                            }

                                        }

                                    }
                                }
                            }

                        }
                    
                        if(rowIsLoneRanger || colIsLoneRanger || subgridIsLoneRanger)
                        {
                            sudokuBoard[i][j] = valueToFind;
                            possibleValuesForSudoku[i][j] = {-1};
                            hasModified = true;
                            break;    
                        }
                    }
                    
                }
                if(hasModified)
                {
                    break;
                }
            }
            if(hasModified)
            {
                updateAllPossibleValues();
            }
            return hasModified;
        }      
};
#endif  