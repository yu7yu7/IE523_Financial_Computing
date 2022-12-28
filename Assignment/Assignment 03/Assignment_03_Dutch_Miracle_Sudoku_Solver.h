#ifndef sudoku
#define sudoku

#include <iostream>
#include <vector>
#include <fstream>
#include <tuple>
#include <set>
#include <stdlib.h>


using std::vector;
using namespace std;
class Sudoku
{
    // Private
    int puzzle[9][9];
    int count = 0 ;


    bool sudoku_valid(int row, int col)
    {
        
        int x = row / 3; 
        int y = col / 3;
        for(int i = 0; i < 9; i++)
            {
                if(puzzle[row][i] == puzzle[row][col] && i != col)
                    return false;
                if(puzzle[i][col] == puzzle[row][col] && i != row)
                    return false;
                }
        for(int i = 3*x; i < 3*x+3; i++)
            for(int j = 3*y; j < 3*y+3; j++)
                    if(puzzle[i][j] == puzzle[row][col] && (i != row && j != col))
                        return false;
        return true;
    }
    
    std::tuple<bool, int> diagonals_valid()
    {
        bool res = true;
        int num;
        
        for(int i=1;i<9 ;i++){
            set<int>x;
            
            for (int row=i, col=0 ; row>=0 && col<9 ; row--, col++)
            {
                if(x.count(puzzle[row][col]) == 0 && puzzle[row][col] != 0)
                {
                    x.insert(puzzle[row][col]);
                    continue;
                }
                if(x.count(puzzle[row][col])!=0 && puzzle[row][col] != 0)
                {
                    res = false;
                    num = puzzle[row][col];
                    break;
                }
            }
            if(! res) break;
        }
        
        if(res)
        {
            for(int i=1;i<8 ;i++){
                set<int>x;
                for (int row=8, col=i ; row>=0 && col<9 ; row--, col++)
                {
                    if(x.count(puzzle[row][col])==0 && puzzle[row][col] != 0)
                    {
                        x.insert(puzzle[row][col]);
                        continue;
                    }
                    if(x.count(puzzle[row][col])!=0 && puzzle[row][col] != 0)
                    {
                        res = false;
                        num = puzzle[row][col];
                        break;
                    }
                }
                if(! res) break;
            }
        }
        return make_tuple(res,num);
    }
    
    std::tuple<bool, int, int> difference_valid()
    {
        bool res = true;;
        int r;
        int c;
        
        for(int i=1;i<9 ;i++){
            for (int row=i, col=0 ; row>=1 && col<8 ; row--, col++)
                if(puzzle[row][col]!=0 && puzzle[row-1][col+1]!=0)
                    if(abs(puzzle[row][col]-puzzle[row-1][col+1]) <4)
                    {
                        res = false;
                        r=row;
                        c=col;
                        break;
                    }
            if(! res) break;
        }
        if(res)
            for(int i=1;i<8 ;i++){
                set<int>x;
                for (int row=8, col=i ; row>=1 && col<8 ; row--, col++)
                    if(puzzle[row][col]!=0 && puzzle[row-1][col+1]!=0)
                        if(abs(puzzle[row][col]-puzzle[row-1][col+1]) <4)
                        {
                            res = false;
                            r=row;
                            c=col;
                            break;
                        }

                if(! res) break;
            }
                
 
        
        return make_tuple(res,r,c);
    }
    
    
public:
    void read_puzzle(int argc, char * const argv[])
    {
        ifstream input_file("/Users/yu-chingliao/Desktop/Assignment_03_Dutch_Miracle_Sudoku_Solver/input5");
        int value;
        vector <int> P;
        if (input_file.is_open())
        {
            while(input_file >> value)
            {
                P.push_back(value);
            }
            for(int i=0;i<9;i++)
                for(int j=0; j<9; j++)
                    puzzle[i][j]=P[j+9*i];
        }
        
        else
            cout << "Input file:" << argv[1] << " does not exist" << endl;
    }
    
    void print_puzzle()
    {
        std::cout << std::endl << "Board Position" << std::endl;
        for (int i = 0; i < 9; i++)
        {
            for (int j = 0; j < 9; j++)
            {
                if ((puzzle[i][j] >= 1) && (puzzle[i][j] <= 9))
                {
                    std::cout << puzzle[i][j] << " ";
                }
                else {
                    std::cout << "X ";
                }
            }
            std::cout << std::endl;
        }
    }
    
    bool Place_Valid(int row, int col){
        return sudoku_valid(row,col)  && get<0>(difference_valid()) &&  get<0>(diagonals_valid());
    }
    
    
    bool Solve(int row, int col)
    {
        
        if(row==0 && col==0)
        {
            if(! get<0>(diagonals_valid()) || ! get<0>(difference_valid())){
                print_puzzle();
                cout << endl;
                cout<<"Partially-Filled Sudoku Puzzle does not meet Dutch Miracle requirement... Quitting!"<<endl;
                
            }
                
            if(! get<0>(diagonals_valid()))
                cout << "The number " << get<1>(diagonals_valid()) << " appears multiple times along the positive-diagonal" << endl;
            if(! get<0>(difference_valid()))
            {
                cout << "puzzle[" << get<1>(difference_valid())  << "]["<< get<2>(difference_valid()) << "] = "<< puzzle[get<1>(difference_valid())][get<2>(difference_valid())] <<" and puzzle["<< get<1>(difference_valid())-1  << "]["<< get<2>(difference_valid())+1 << "] = "<< puzzle[get<1>(difference_valid())+1][get<2>(difference_valid())-1] << " has a diffrence smaller than 4"<< endl;
                return false;
            }
            if(! get<0>(diagonals_valid())) return false;
            
            cout << " " << endl;
            cout <<"Initial Sudoku Puzzle meets Dutch Miracle requirments" << endl;
            print_puzzle();
            cout << endl;
            cout << "Enumerating all solutions to the Dutch Mriacle Sudoku instance shown above" << endl;
            
        }
        
        if (row == 9)
        {
            count ++;
            cout << endl;
            cout <<"Solution #" << count;
            print_puzzle();
            return true;
        }
        bool res = false;
        
        if(0 != puzzle[row][col])
        {
              if (col ==8)
              {
                  res = Solve(row+1,0) || res;
              }
              else
                  res = Solve(row, col+1) || res;
              return res;
        }

        for (int num = 1; num < 10; num++)
        {
            puzzle[row][col]=num;
            if (! Place_Valid(row, col))
                continue;
            if (col ==8)
            {
                res = Solve(row+1,0) || res;
                
                
            }
            else
            {
                res = Solve(row, col+1) || res;
                
                
            }
            
        }
        puzzle[row][col]=0;
        return res;
    }
};

#endif

