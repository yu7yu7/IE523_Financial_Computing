//
//  Assignment_07_Simulation_Part1.cpp
//  Assignment_07_Simulation
//
//  Created by Yu-Ching Liao on 10/19/22.
//
#include <cmath>
#include <iostream>
using namespace std;

class Coin_Toss_Simulator{
    int n; //Number of games
    double q; //Probability for Alice
    double p; //Probability for Bob
    
    long long NCR(int n, int r) {
        if(r > n - r) r = n - r; // because C(n, r) == C(n, n - r)
        long long ans = 1;
        int i;

        for(i = 1; i <= r; i++) {
            ans *= n - r + i;
            ans /= i;
        }

        return ans;
    }

public:
    double Prob_of_Alice_Win(int n, double q, double p){
        double f = 0;
        for (int r = 0; r <= n; r++){
            for (int s = r+1; s <= n; s ++){
                f += NCR(n, r)*pow(p, r)*pow(1-p, n-r)*NCR(n,s)*pow(q, s)*pow(1-q, n-s);
            }
        }
        return f;
    }
};



int main() {
    int n = 2;
    double q;
    double p;
    
    double x;
    double x_prev;
    double x_next;
    bool tester = false;
    
    cout << "Probability of head up for Alice = ";
    cin >> q;
    cout << "Probability of head up for Bob = ";
    cin >> p;
    
    Coin_Toss_Simulator X;
    
    while (!tester){
        x_prev = X.Prob_of_Alice_Win(n-1, q, p);
        x = X.Prob_of_Alice_Win(n, q, p);
        x_prev = X.Prob_of_Alice_Win(n+1, q, p);
        if ((x_prev <= x) && (x >= x_next) && ((x != x_prev)||(x != x_next))){
            cout << "The optimal number of coin toss is " << n << "." <<endl;
            tester = true;
        }
        n += 1;
    }
    

}
