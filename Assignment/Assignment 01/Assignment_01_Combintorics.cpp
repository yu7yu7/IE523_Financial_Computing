//
//  Combinatorics_Assignment.cpp
//  Assignment 01
//
//  Created by Yu-Ching Liao on 2022/8/31.
//


#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

int main(int argc,char* argv[])
{
    int size_of_array, number_of_objects_to_be_picked;
    int *array;
    
    sscanf(argv[1],"%d",&size_of_array);
    sscanf(argv[2],"%d",&number_of_objects_to_be_picked);
    

    array = new int[size_of_array];
    for (int i = 0; i < size_of_array; i++)
        array[i] = i+1;



    vector <bool> did_i_pick_this(size_of_array);
    fill(did_i_pick_this.begin()+ number_of_objects_to_be_picked,
         did_i_pick_this.end(), true);
    

    do
    {
        cout << "{ ";
        for (int i = 0; i < size_of_array; ++i)
        {
            if (did_i_pick_this[i]==false)
                cout << array[i] << " ";
        }
        cout << "}" << endl;
    } while (next_permutation(did_i_pick_this.begin(), did_i_pick_this.end()));
        
    return 0;
}
