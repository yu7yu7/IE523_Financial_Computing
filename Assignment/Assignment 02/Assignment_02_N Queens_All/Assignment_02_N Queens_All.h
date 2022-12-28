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
    int c;



    bool is_this_position_safe(int col, int row, int size){

        for (int i = 0; i < row; i++) {
            if (chess_board[i][col]) {
                return false;
            }
        }

        for (int i = row,j = col;i >= 0 && j >= 0; i--,j--) {
            if (chess_board[i][j]) {
                return false;
            }
        }

        for (int i = row, j = col; i >= 0 && j < size; j++, i--) { 
            if (chess_board[i][j]) {
                return false;
            }
        }
        return true;
    }

    void initialize(int n){
        size = n;
        c = 0;
        chess_board = new int*[n];
        for(int i = 0; i < size; i++)
            chess_board[i] = new int[size];
    }

    void print_board(int size){
        //print board by iterate the array here
        cout << size << "X"<< size << " Solution #: " << c << endl;
        for(int i = 0; i < 2*size ; i++)
            cout << "-";
            cout << endl;
        for (int i = 0;i <= size-1; i++) {
            for (int j = 0;j <= size-1; j++) {
                if (chess_board[i][j] == 0)
                    cout <<"-"<< " ";
                else
                    cout <<"Q"<<" ";
                
            }
            cout<<endl;
        }
        for(int i = 0; i < 2*size ; i++)
            cout << "-" ;
            cout << endl;
        cout<<endl;

    }

    bool solve (int size, int row) {
        if (size == row) {
            c ++;
            print_board(size);
            return true;
        }

        bool res = false;
        for (int i = 0;i <=size-1;i++) {
            if (is_this_position_safe(i, row, size)) {
                chess_board[row][i] = 1;
                res = solve(size, row+1) || res;
                chess_board[row][i] = 0;
//                print_board(size);
            }
            print_board(size);
        }
        return res;
    }

public:
    void nQueens(int n){
        initialize(n);
        if (n>3)
            cout << "The solutions of " << n << "-Queens Problem are:" << endl;
        
        bool res = solve(n, 0);
    
        if(res == false) {
            cout << "There is no possible solution." << endl;
        } else {
            cout<<"There are "<<c<<" different solutions to the "<<n<<"-Queens Problem"<<endl;
            cout << endl;
        }
    }
};

#endif //N_queens
