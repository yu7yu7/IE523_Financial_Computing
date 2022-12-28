#include "pbPlots.h"
#include "supportLib.h"
#include <cmath>
#include <iostream>
#include <random>

using namespace std;
unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);

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
    

    


    double get_uniform()
    {
        // http://www.cplusplus.com/reference/random/exponential_distribution/
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        double number = distribution(generator);
        return (number);
    }
    
};


int main(){
    int n = 2;
    int no_of_trials;
    double q;
    double p;
    int countAlice;
    int countBob;
    int countAwin;
    

    double testerAlice;
    double testerBob;
    bool success;
    bool examT = true;
    bool examE = true;
    
    vector<double> xs;
    vector<double> ys;
    vector<double> storage;
    vector<double> alice;
    vector<double> bob;
    vector<double> xs2;
    vector<double> ys2;
    
    cout << "Number of games = ";
    cin >> n;
    cout << "Probability of head up for Alice = ";
    cin >> q;
    cout << "Probability of head up for Bob = ";
    cin >> p;
    cout << "Number of trials = ";
    cin >> no_of_trials;
    
    if (q>=p){
        cout << "Error: The probability of Alice getting head is larger than that of Bob's."<<endl;
        return 0;
    }
    
    if ((q>=1) || (p > 1)){
        cout << "Error: Probability should not be greater than 1. " << endl; 
        return 0;
    }
    
    Coin_Toss_Simulator X;
    
    for (int i = 1; i <= n; i++){
        vector<double> storage = {};
        
        ys.push_back(X.Prob_of_Alice_Win(i, q, p));
        xs.push_back(i);
    }


    for (int i = 1; i <= n; i++){
        countAwin = 0;
        for (int k = 0; k <= no_of_trials; k++){
            countAlice = 0;
            countBob = 0;
            for (int j = 1; j <= i; j++){
                testerAlice = X.get_uniform();
                testerBob = X.get_uniform();
                if (testerAlice < q) countAlice += 1; //calculate the number of Alice getting head
                if (testerBob < p) countBob += 1; //calculate the number of Bob getting head
            }
            if (countAlice > countBob) countAwin += 1;
        }
        ys2.push_back(double(countAwin)/double(no_of_trials));
        xs2.push_back(i);
        
    }
    for (int i = 0; i < n; i++){
        if ((ys[i] == *max_element(ys.begin(), ys.end())) && (examT)){
            cout << "Theoretically, " << xs[i] << " is the best number of games for Alice." << endl;
            cout << "The theoretical possibility of Alice winning the game is " << ys[i] << "." << endl;
            examT = false;
        }
        if ((ys2[i] == *max_element(ys2.begin(), ys2.end())) && (examE)) {
            cout << "Experimentally, " << xs2[i] << " is the best number of games for Alice." << endl;
            cout << "The experimental possibility of Alice winning the game is " << ys2[i] << "." << endl;
            examE = false;
        }
    }
    
    
    
    
    StringReference *errorMessage = new StringReference();
    RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();
    
    ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
    series->xs = &xs2;
    series->ys = &ys2;
    series->linearInterpolation = true;
    series->pointType = toVector(L"dots");
    series->color = CreateRGBColor(1, 0, 0);
    
    ScatterPlotSeries *series2 = GetDefaultScatterPlotSeriesSettings();
    series2->xs = &xs;
    series2->ys = &ys;
    series2->linearInterpolation = true;
    series2->lineType = toVector(L"solid");
    series2->lineThickness = 2;
    series2->color = CreateRGBColor(0, 0, 1);
    
    ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
    settings->width = 1080;
    settings->height = 720;
    settings->autoBoundaries = true;
    settings->autoPadding = true;
    settings->title = toVector(L"Blue line is Theoretical, and Red line is Experimental. ");
    settings->xLabel = toVector(L"Number of Coin Tosses in Each Game");
    settings->yLabel = toVector(L"Alice's Probability of Winning");
    settings->scatterPlotSeries->push_back(series);
    settings->scatterPlotSeries->push_back(series2);
    
    
    
    
    success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);
    
    if(success){
        vector<double> *pngdata = ConvertToPNG(imageReference->image);
        WriteToFile(pngdata, "/Users/yu-chingliao/Desktop/Best number of games for Alice.png");
        DeleteImage(imageReference->image);
    }else{
        cerr << "Error: ";
        for(wchar_t c : *errorMessage->string){
            wcerr << c;
        }
        cerr << endl;
    }
    
    return success ? 0 : 1;
}
