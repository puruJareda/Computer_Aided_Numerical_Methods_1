#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

vector<long double> newton_divided_difference(vector<long double>& x, vector<long double>& y){

    int n = x.size();

    //constructing divided difference table
    vector<vector<long double>> table(n, vector<long double>(n));

    //first column = y values
    for(int i=0;i<n;i++){
        table[i][0] = y[i];
    }

    //filling table using divided difference formula
    for(int j=1; j<n; j++){
        for(int i =0; i<n-j;i++){
            table [i][j] = (table[i+1][j-1] - table[i][j-1]) / (x[i+j] - x[i]);
        }
    }

    //extracting newton coefficients
    vector<long double> coeff(n);
    for(int j=0; j<n;j++){
        coeff[j] = table[0][j];
    }

    //converting to standard polyomial
    vector<long double> coeff_standard(n, 0);
    vector<long double> term = {1}; // represents current product term

    // add contribution of current term : a_i * T_i(x)
    for(int i=0;i<n;i++){
        for(int j=0;j<term.size() && j<n;j++){
            coeff_standard[j] += coeff[i] * term[j];
        }

        if(i == n-1){
            break;
        }

        //multiply current term by (x - x_i)
        vector<long double> new_term (term.size()+1, 0);//increasing the size of term to store the coefficient of x with higher power 

        for(int j=0;j<term.size(); j++){
            new_term[j] -= term[j] *x[i]; // -x_i * term
            new_term[j+1] += term[j]; // multiply by x
        }

        term = new_term; // update term
    }

    //printing standard polynomial
    cout << "\n divided difference Polynomial in standard form"<< endl;
    cout << "\n P(x) = ";
    
    long double EPS = 1e-12;
    bool first = true;

    for(int i = 0; i < n; i++) {
        if(fabsl(coeff_standard[i]) < EPS){
            coeff_standard[i] = 0;
            continue; 
        }

        if(!first){
            if(coeff_standard[i] > 0){
                cout << " + ";
            }
            else{
                cout << " - ";
            }
        } 
        else{
            if(coeff_standard[i] < 0){
                cout << "-";
            }
            first = false;
        }

        cout << fabsl(coeff_standard[i]);

        if(i >= 1) cout << "x";
        if(i >= 2) cout << "^" << i;
    }

    cout << endl;

    return coeff_standard;
}