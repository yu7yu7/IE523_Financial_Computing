//
// Created by LYC on 2022/09/07.
//

#ifndef N_queens
#define N_queens
#include <set>
using namespace std;

class Board{
    int size;
    int **chess_board;


    bool is_this_position_safe(int row,int col){
        int i, j;
      
        for (i = 0; i < col; i++)
            if (chess_board[row][i])
                return false;
      
        for (i = row, j = col; i >= 0 && j >= 0; i--, j--)
            if (chess_board[i][j])
                return false;
      
        for (i = row, j = col; j >= 0 && i < size; i++, j--)
            if (chess_board[i][j])
                return false;
      
        return true;
    }

    void initialize(int n){
        size = n;
        chess_board = new int*[n];
        for(int i = 0; i < size; i++)
            chess_board[i] = new int[size];
    }

    void print_board(){
        //print board by iterate the array here
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++)
                if(chess_board[i][j]) cout <<"Q"<<" ";
                else
                    cout <<"-"<<" ";
            cout<< endl;
        }
    }

    bool solve(int col){
        if (col>=size) return true;
        for (int i = 0; i < size; i++) {
            if (is_this_position_safe(i, col)) {
                chess_board[i][col] = 1;
                if (solve(col + 1))
                    return true;
                chess_board[i][col] = 0;
            }
        }
        return false;
    }

public:
    void nQueens(int n){
        initialize(n);
        if(solve(0)) {
            cout << n << "-Queens Problem Soultion" << endl;
            for(int i = 0; i < 2*n ; i++)
                cout << "-";
            cout << endl;
            print_board();
            for(int i = 0; i < 2*n ; i++)
                cout << "-";
                cout <<endl;
        }
          else cout<<"There is no solution to the " << n <<" Queens Problem"<< endl;
    }
};

#endif //N_queens
