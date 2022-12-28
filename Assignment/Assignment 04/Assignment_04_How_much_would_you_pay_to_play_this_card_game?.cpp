//
//  Assignment_04_How_much_would_you_pay_to_play_this_card_game?.cpp
//  Assignment_04_How_much_would_you_pay_to_play_this_card_game?
//
//  Created by Yu-Ching Liao on 9/15/22.
//

#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;
double **cache;

void initialize(int n){
    cache = new double*[n];
    for(int i = 0; i < n; i++)
        cache[i] = new double[n];
}

double value(int red, int black){
    double r_f = red, b_f = black, s_f = red+black;
    double p_of_red = r_f/s_f;
    double p_of_black = b_f/s_f;
    
    if (cache[red][black] != 0)
        return cache[red][black];
    
    if (0 == red)
        return ((double)black);
    if (0 == black)
        return (0);
    vector<double> val = {p_of_red * value(red-1, black)+p_of_black * value(red, black-1), b_f - r_f};
    
    if (val[0] > val[1]){
        cache[red][black] = val[0];
        return val[0];
    }
    else{
        cache[red][black] = val[1];
        return val[1];
    }
    
    
}

int main(int argc,char* argv[])
{
    int size_of_card, red_cards_left, black_cards_left;
    
    cout << "Total Number of Cards = ";
    cin >> size_of_card;
    if(size_of_card % 2 != 0)
        cout << "Please try again." << endl;
    else{
        red_cards_left = size_of_card/2;
        black_cards_left = size_of_card/2;
        initialize(size_of_card);
        }
    cout << "Value of the game = " << value(red_cards_left, black_cards_left) << endl;;
}

