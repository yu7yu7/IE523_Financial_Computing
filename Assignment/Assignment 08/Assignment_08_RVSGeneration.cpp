#include "pbPlots.h"
#include "supportLib.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
double pi = 3.1415926535897932384626433;
#define CAUCHY_DENSITY(x) pow(pi*(1+x*x) ,-1)

using namespace std;
unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);

double get_uniform()
{
    // http://www.cplusplus.com/reference/random/exponential_distribution/
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

double get_cauchy(){
    return tan(pi*(get_uniform()-0.5));
}

int main (int argc, char* argv[])
{
    float y;
    int no_of_trials, count[100];
    vector<double> xs;
    vector<double> ys;
    vector<double> xs2;
    vector<double> ys2;
    sscanf (argv[1], "%d", &no_of_trials);
    ofstream pdf_comparison_file(argv[1]);
    
    bool success;
    
    for (int i = 0; i < 100; i++) {
        count[i] = 0;
    }
    
    for (int i = 0; i < no_of_trials; i++) {
        y = get_cauchy();
        for (int j = 0; j < 100; j++)
            if ( (y >= ((float) (j-51)/10)) && (y < ((float) (j-50)/10)) )
                count[j]++;
    }
    
    
    int sum = 0;
    for (int j = 0; j < 100; j++) {
        sum += count[j];
        pdf_comparison_file << ((float) (j-50)/10) << ", " <<
        ((float) count[j]/no_of_trials) << ", " <<
        (0.1*CAUCHY_DENSITY((float) (j-50)/10)) << endl;
        
        xs.push_back(((float) (j-50)/10));
        xs2.push_back(((float) (j-50)/10));
        ys.push_back(((float) count[j]/no_of_trials));
        ys2.push_back((0.1*CAUCHY_DENSITY((float) (j-50)/10)));
        
    }
    float average = accumulate( ys.begin(), ys.end(), 0.0)/ ys.size();
    cout << average << endl;
    
    StringReference *errorMessage = new StringReference();
    RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();
    
    ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
    series->xs = &xs;
    series->ys = &ys;
    series->linearInterpolation = true;
    series->lineType = toVector(L"solid");
    series->lineThickness = 2;
    series->color = CreateRGBColor(1, 0, 0);
    
    ScatterPlotSeries *series2 = GetDefaultScatterPlotSeriesSettings();
    series2->xs = &xs2;
    series2->ys = &ys2;
    series2->linearInterpolation = true;
    series2->lineType = toVector(L"solid");
    series2->lineThickness = 2;
    series2->color = CreateRGBColor(0, 0, 1);
    
    ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
    settings->width = 1080;
    settings->height = 720;
    settings->autoBoundaries = true;
    settings->autoPadding = true;
    settings->title = toVector(L"Blue line is Theoritical, and Red line is Experimental.");
    settings->xLabel = toVector(L"X");
    settings->yLabel = toVector(L"CDF");
    settings->scatterPlotSeries->push_back(series);
    settings->scatterPlotSeries->push_back(series2);
    
    success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);
    
    if(success){
        vector<double> *pngdata = ConvertToPNG(imageReference->image);
        WriteToFile(pngdata, "/Users/yu-chingliao/Desktop/Cauchy RVs Generation.png");
        DeleteImage(imageReference->image);
    }else{
        cerr << "Error: ";
        for(wchar_t c : *errorMessage->string){
            wcerr << c;
        }
        cerr << endl;
    }
    
}
