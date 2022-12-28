//  Gurobi Read Model From File
//
//  Created by Ramavarapu Sreenivas on 8/2/22.
//  from
//  https://www.gurobi.com/documentation/9.5/examples/lp_cpp_cpp.html
//  https://www.gurobi.com/documentation/9.5/refman/cpp_attribute_examples.html
//  https://github.com/fuminori-ido/gurobi/blob/master/ext/gurobi/model.cpp
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gurobi_c++.h"
using namespace std;

int main(int argc, const char * argv[])
{
    if (argc < 2)
    {
        cout << "Usage: lp_c++ filename" << endl;
        return 1;
    }
    try
    {
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env, argv[1]);
        model.optimize();
        int optimstatus = model.get(GRB_IntAttr_Status);
        if (optimstatus == GRB_INF_OR_UNBD) {
          model.set(GRB_IntParam_Presolve, 0);
          model.optimize();
          optimstatus = model.get(GRB_IntAttr_Status);
        }
        if (optimstatus == GRB_OPTIMAL) {
            double objval = model.get(GRB_DoubleAttr_ObjVal);
            cout << "Optimal objective: " << objval << endl;
            cout << "Optimal Values:" << endl;
            GRBVar* vars = model.getVars();
            int i = 0;
            for(GRBVar *p = vars; i < model.get(GRB_IntAttr_NumVars); i++, p++)
                printf("%s = %f\n", p->get(GRB_StringAttr_VarName).c_str(),
                       p->get(GRB_DoubleAttr_X));
        } else if (optimstatus == GRB_INFEASIBLE) {
          cout << "Model is infeasible" << endl;
          // compute and write out IIS
          model.computeIIS();
          model.write("model.ilp");
        } else if (optimstatus == GRB_UNBOUNDED) {
          cout << "Model is unbounded" << endl;
        } else {
          cout << "Optimization was stopped with status = "
               << optimstatus << endl;
        }
      } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Error during optimization" << endl;
      }    
    return 0;
}
