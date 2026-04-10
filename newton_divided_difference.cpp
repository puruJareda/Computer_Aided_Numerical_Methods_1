#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

using namespace std;

int main(){

    //read data
    ifstream file("group2 canm dataset.csv");

    vector<long double> x,y;

    long double x_val, y_val;
    char comma;

    file.ignore(1000, '\n');

    while(file >> x_val >> comma >> y_val) {
        x.emplace_back(x_val);
        y.emplace_back(y_val);
    }
  
    file.close();

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

    //printing the divided difference table
    cout << "\n Divided Difference Table: \n";
    for(int i=0; i<n;i++){
        for(int j=0;j<n-i;j++){
            cout << table[i][j] << "\t";
        }
        cout << endl;
    }

    //printing the newton coefficients
    cout << "\n Coefficients:\n";
    for(int i=0;i<n;i++){
        cout << "a" << i << "=" << coeff[i] << endl;
    }

    //printing newton form of polynomial
    cout << "\n P(x) = ";
    cout << coeff[0];

    for(int i=1; i<n;i++){
        char sign1;
        if(-coeff[i] > 0){
            sign1 = '-';
        }
        else{
            sign1 = '+';
        }
        cout <<" " <<sign1 << " " <<fabsl(coeff[i]);

        for(int j=0;j<i;j++){
            char sign2;
            if(-x[j] > 0){
                sign2 = '-';
            }
            else{
                sign2 = '+';
            }
            cout << "(x" << sign2 << fabsl(x[j]) << ")";
        }
    }

    cout << endl;

    //converting to standard polyomial
    vector<long double> coeff_standard(n, 0); // final polynomial coefficients
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

    long double EPS = 1e-12; // removing floating point noise

    // export the standard coefficient file
    ofstream out_std("standard_polynomial_coeff.csv");
    out_std<< "power,coefficient\n";

    for(int i=0;i<n;i++){
        if(fabsl(coeff_standard[i]) < EPS){
            coeff_standard[i] = 0;
        }
        out_std<<i<<","<< coeff_standard[i] << endl;
    }
    out_std.close();

    //printing standard polynomial
    cout << "\n Polynomial in standard form"<< endl;
    cout << "\n P(x) = ";
    
    bool first = true;

    for(int i = 0; i < n; i++) {
        if(fabsl(coeff_standard[i]) < EPS){
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

    //export the values of the approximated polynomial at points inside the domain for plotting and visualization.
    ofstream out("divided.csv");

    long double step = (x[n-1] - x[0])/ 200.0;

    for(long double xi = x[0]; xi <= x[n-1]; xi += step){
        long double result = coeff[0];
        long double term = 1;

        for(int i=1;i<n;i++){
            term *= (xi - x[i-1]);
            result += coeff[i] * term;
        }

        out << xi << ","<<result << endl;
    }

    out.close();


    return 0;
}