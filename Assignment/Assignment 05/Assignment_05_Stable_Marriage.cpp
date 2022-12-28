#include <iostream>
#include <vector>
#include <fstream>
#include <map>
using namespace std;

class stable_marriage_instance
{
    int no_of_couples;
    vector <vector <int> > Preference_of_men;
    vector <vector <int> > Preference_of_women;
    vector <int> match_for_men;
    vector <int> match_for_women;

    int anybody_free(vector <bool> my_array)
    {
        for (int i = 0; i<my_array.size(); i++) if (my_array[i])
            return i;
        return -1;
    }

    bool rank_check (vector <int> my_array, int i1, int i2)
    {
        for (int i; i<my_array.size();i++){
            if (my_array[i] == i1)
                return true;
            if (my_array[i] == i2)
                return false;
        }
        return false;
    }

    void Gale_Shapley()
    {
        vector <bool> is_man_free;
        vector <bool> is_woman_free;
        vector <vector <bool> > has_this_man_proposed_to_this_woman;

        int man_index, woman_index;

        for (int i= 0; i < no_of_couples; i++){
            is_man_free.push_back(1);
            is_woman_free.push_back(1);
            vector<bool> currMan_proposed_to_this_woman;
            for (int j = 0; j < no_of_couples ; j++) {
                currMan_proposed_to_this_woman.push_back(0);
            }
            has_this_man_proposed_to_this_woman.push_back(currMan_proposed_to_this_woman);
        }

        // Gale-Shapley Algorithm
        while ( (man_index = anybody_free(is_man_free)) >= 0)
        {
            if (!is_man_free[man_index]){
                continue;
            }

            for (int currManPrefIndx = 0; currManPrefIndx < Preference_of_men[man_index].size();currManPrefIndx++) {
                woman_index = Preference_of_men[man_index][currManPrefIndx];
                if (has_this_man_proposed_to_this_woman[man_index][woman_index]){
                    continue;
                }
                else{
                    break;
                }
            }
            has_this_man_proposed_to_this_woman[man_index][woman_index] = 1;
            if (is_woman_free[woman_index]){
                is_woman_free[woman_index] = 0;
                is_man_free[man_index] = 0;
                match_for_men[man_index] = woman_index;
                match_for_women[woman_index] = man_index;
            }
            else {
                if (rank_check(Preference_of_women[woman_index],man_index,match_for_women[woman_index])){
                    is_man_free[match_for_women[woman_index]] = 1;
                    is_woman_free[woman_index] = 0;
                    is_man_free[man_index] = 0;
                    match_for_men[man_index] = woman_index;
                    match_for_women[woman_index] = man_index;
                }
            }
        }
    }


    void read_data(int argc, const char * argv[])
    {
        ifstream input_file("/Users/yu-chingliao/Library/CloudStorage/GoogleDrive-josephliao0127@gmail.com/My Drive/Note/UIUC/Fall_2022/Financial_Computing/Assignment/Assignment_05_Stable_Marriage_Problem/input3");
        int input_value;
        if (input_file.is_open()){
            input_file >> input_value;
            no_of_couples = input_value;
            for (int i = 0; i<no_of_couples; i++){
                match_for_men.push_back(-1);
                vector<int> currMenPref;
                for (int j = 0; j < no_of_couples; j++){
                    input_file >> input_value;
                    currMenPref.push_back(input_value);
                }
                Preference_of_men.push_back(currMenPref);
            }
            for (int i = 0; i<no_of_couples; i++){
                match_for_women.push_back(-1);
                vector<int> currWomenPref;
                for (int j = 0; j < no_of_couples; j++){
                    input_file >> input_value;
                    currWomenPref.push_back(input_value);
                }
                Preference_of_women.push_back(currWomenPref);
            }
        }
    }

    void print_soln()
    {
        cout << "Preferences of Men" << endl;
        cout << "------------------" << endl;
        cout << endl;

        for (int i = 0; i < Preference_of_men.size(); i++){
            cout << "(" << i << "): ";
            for (int j = 0; j<Preference_of_men[i].size(); j++){
                cout << Preference_of_men[i][j] <<" ";
            }
            cout << endl;
        }

        cout << endl;
        cout << "Preferences of Women"<< endl;
        cout << "--------------------" << endl;
        cout << endl;

        for (int i = 0; i<Preference_of_women.size(); i++){
            cout << "(" << i << "): ";
            for (int j = 0; j<Preference_of_women[i].size(); j++){
                cout << Preference_of_women[i][j] <<" ";
            }
            cout << endl;
        }
        cout << endl;
        cout << "Matching for Men" <<endl;
        cout << endl;

        for (int i = 0; i<match_for_men.size(); i++) cout << "Man: "<<i<< " -> Woman: "<< match_for_men[i] << endl;
        cout << endl;
        cout << "Matching for Women" << endl;
        cout << endl;
        for (int i = 0; i<match_for_women.size(); i++) cout << "Woman: "<<i<< " -> Man: "<< match_for_women[i] << endl;

    }

public:
    void solve_it(int argc, const char * argv[])
    {
        read_data(argc, argv);
        Gale_Shapley();
        print_soln();
    }
};

int main(int argc, const char * argv[])
{
    stable_marriage_instance x;
    x.solve_it(argc, argv);
}
