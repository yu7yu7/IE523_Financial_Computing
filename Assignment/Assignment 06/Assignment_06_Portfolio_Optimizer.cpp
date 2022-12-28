#include "gurobi_c++.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
using namespace std;

class Optimizer{
    const double min_error = 1e-3;

    vector<double> PV;
    vector<int> M;
    vector<vector<double>> cashFlows;
    vector<double> convexity;
    vector<double> duration;
    vector<double> ytm;

    double debt;
    double Duration;
    int num_of_cashflows;
    
    double F(double r, double m, vector<double> cf, double pv)
    {
        double discountCf = 0;
        for(int t = 1; t <= (int) m; t++)
        {
            double p = m - (double) t;
            discountCf += cf[t-1] * pow(1+r,p);
        }
        double res = pow((1+r), m)*pv-discountCf;
        return res;
    }
    
    double dF(double r, double m, vector<double> cf, double pv)
    {
        double Sum = 0;
        for(int t = 1;t <= (int)m- 1; t++)
        {
            Sum += cf[t - 1] * (m - (double) t) * pow(1+r, m-(double)t - 1);
        }
        double res=(m * pow(1 + r, m - 1)) * pv - Sum;
        return res;
    }
    
    double YTM(double r, double m, vector<double> cf, double pv)
    {
        while(double error = abs(F(r, m, cf, pv)) > min_error)
        {
            r = r-F(r,m,cf,pv)/dF(r,m,cf,pv);
        }
        return r;
    }
    
    double D(double r, double m, vector<double> cf, double pv)
    {
        double res = 0;
        for(int t = 1; t <= m; t++)
        {
            res += (cf[t - 1] / pow(1+r,t) * t) / pv;
        }
        return res;
    }
    
    double C(double r, double m, vector<double> cf, double pv)
    {
        double res = 0;
        for(int t = 1; t <= m; t++)
        {
            res += (t*(t+1) * cf[t-1])/pow(1+r,t+2);
        }
        return res / pv;
    }
    
    void readData()
    {
        int m;
        double PV_plan, cf;
        
        ifstream inputFile("/Users/yu-chingliao/Library/CloudStorage/GoogleDrive-josephliao0127@gmail.com/My Drive/Note/UIUC/Fall_2022/Financial_Computing/Assignment/Assignment_06_Portfolio_Optimization/Assignment_06_Portfolio_Optimization/input1");
        if(inputFile.is_open())
        {
            inputFile >> num_of_cashflows;
            
            for(int i = 0; i < num_of_cashflows; i++) {
                inputFile >> PV_plan;
                PV.push_back(PV_plan);
                inputFile >> m;
                M.push_back(m);
                vector<double> CF1;
                for(int j = 0; j < m; j++) {
                    inputFile >> cf;
                    CF1.push_back(cf);
                }
                cashFlows.push_back(CF1);
                double r = YTM(0, (double)m, CF1, PV_plan);
                duration.push_back(D(r, (double)m, CF1, PV_plan));
                convexity.push_back(C(r, (double)m, CF1, PV_plan));
                ytm.push_back(r);
            }
            inputFile >> debt;
            inputFile >> Duration;
        }
    }

    
    
public:
    void printResult()
    {
        readData();
        cout << "We owed " << debt << " in " << Duration << " years" << endl;
        cout << "Numnber of Cash flow: " << num_of_cashflows << endl;
        for(int n = 0; n < num_of_cashflows; n++)
        {
            cout << "-------------------------------" << endl;
            cout << "Cash Flow #" << n + 1 << endl;
            cout << "Price = " << PV[n] << endl;
            cout << "Maturity = " << M[n] << endl;
            cout << "Yield to Matruity = " << ytm[n] << endl;
            cout << "Duration = " << duration[n] << endl;
            cout << "Convexity = " << convexity[n] << endl;
        }
        
    }
    
    void GurobiOptimizer() {

        cout << "********************************" << endl;
        try {
            GRBEnv env = GRBEnv(true);
            env.set("LogFile", "PortfolioOptimization.log");
            env.start();

            GRBModel model = GRBModel(env);
            

            GRBVar* lambda = new GRBVar[num_of_cashflows];
            
            for (int i = 0; i < num_of_cashflows; i++) {
                lambda[i] = model.addVar(0.0, 1, 0.0, GRB_CONTINUOUS , "lambda_" + to_string(i));
            }

            GRBLinExpr le1 = 0;
            GRBLinExpr le2 = 0;
            GRBLinExpr bConvex = 0;
            for (int i = 0; i < num_of_cashflows; i++) {
                le1 += lambda[i];
                le2 += lambda[i] * duration[i];
                bConvex += lambda[i] * convexity[i];
            }
            model.setObjective(bConvex, GRB_MAXIMIZE);

            model.addConstr(le1 == 1);
            model.addConstr(le2 == Duration);

            model.optimize();

            int optimistatus = model.get(GRB_IntAttr_Status);

            cout << "*************************************************************************" << endl;

            if (optimistatus == GRB_OPTIMAL) {

                double objval = model.get(GRB_DoubleAttr_ObjVal);
                printf("Largest Convexity we can get is %.3f \n",objval);
                cout << "*************************************************************************" << endl;
                cout << "To immunize against small changes in 'r' for each $1 of PV, you should buy" << endl;
                GRBVar* vars = model.getVars();
                int i = 0;
                for (GRBVar* p = vars; i < model.get(GRB_IntAttr_NumVars); i++, p++)
                    if (p->get(GRB_DoubleAttr_X) > 0) printf("$%f of Cash-Flow %d \n", p->get(GRB_DoubleAttr_X), i+1);
                cout << "If you need to immunize for a larger PV-value, just buy at an approprtiate proportion." << endl;
                cout << "*************************************************************************" << endl;
                cout << "For example, if you want to immunize for $500 of PV, buy" << endl;
                int j = 0;
                for (GRBVar* k = vars; j < model.get(GRB_IntAttr_NumVars); j++, k++)
                    if (k->get(GRB_DoubleAttr_X) > 0) printf("$500 x %f of Cash-Flow %d \n", k->get(GRB_DoubleAttr_X), j+1);
                
            }
            else if (optimistatus == GRB_INFEASIBLE) {
                printf("There is no portfolio that meets the duration constraint of %.1f years", Duration);
                model.computeIIS();
                model.write("model.ilp");
            }
            else {
                cout << "Optimization was stopped with status = " << optimistatus << endl;
            }
            cout << "*************************************************************************" << endl;
        }
        catch (GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }
        catch (...) {
            cout << "Error during optimization" << endl;
        }
        
    }
};


int main()
{
    Optimizer X;
    X.printResult();
    X.GurobiOptimizer();
    return 0;
}


