#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

using namespace std;

//read data
void read_data(vector<long double>& x, vector<long double>& y){
    ifstream file("group2 canm dataset.csv");

    long double x_val, y_val;
    char comma;

    file.ignore(1000, '\n');

    while(file >> x_val >> comma >> y_val) {
        x.emplace_back(x_val);
        y.emplace_back(y_val);
    }
  
    file.close();
}

//constructing divided difference table
void table(vector<vector<long double>>& table_matrix,
           vector<long double>& x,
           vector<long double>& y){

    int n = x.size();

    //first column = y values
    for(int i=0;i<n;i++){
        table_matrix[i][0] = y[i];
    }

    //filling table using divided difference formula
    for(int j=1; j<n; j++){
        for(int i =0; i<n-j;i++){
            table_matrix [i][j] = (table_matrix[i+1][j-1] - table_matrix[i][j-1]) / (x[i+j] - x[i]);
        }
    }
}

//extracting newton coefficients
void get_coeff(vector<long double>& coeff,
               vector<vector<long double>>& table){

    int n = table.size();

    for(int j=0; j<n;j++){
        coeff[j] = table[0][j];
    }
}

//converting to standard polyomial
void convert_to_standard(vector<long double>& coeff_standard,
                         vector<long double>& coeff,
                         vector<long double>& x){

    int n = coeff.size();

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
}

// export the standard coefficient file
void export_standard(vector<long double>& coeff_standard){

    long double EPS = 1e-12; // removing floating point noise

    ofstream out_std("standard_polynomial_coeff.csv");
    out_std<< "power,coefficient\n";

    for(int i=0;i<coeff_standard.size();i++){
        if(fabsl(coeff_standard[i]) < EPS){
            coeff_standard[i] = 0;
        }
        out_std<<i<<","<< coeff_standard[i] << endl;
    }
    out_std.close();
}

//export the values of the approximated polynomial at points inside the domain for plotting and visualization.
void export_values(vector<long double>& x,
                   vector<long double>& coeff){

    int n = x.size();

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
}

int main(){

    vector<long double> x,y;

    read_data(x,y);

    int n = x.size();

    vector<vector<long double>> divided_difference_table(n, vector<long double>(n));

    table(divided_difference_table, x, y);

    vector<long double> coeff(n);
    get_coeff(coeff, divided_difference_table);

    //printing the divided difference table
    cout << "\n Divided Difference Table: \n";
    for(int i=0; i<n;i++){
        for(int j=0;j<n-i;j++){
            cout << divided_difference_table[i][j] << "\t";
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

     //converting to standard polynomial and export the coefficients of the standard polynomial 
    vector<long double> coeff_standard(n, 0);
    convert_to_standard(coeff_standard, coeff, x);
    export_standard(coeff_standard);

    //printing standard polynomial
    cout << "\n Polynomial in standard form"<< endl;
    cout << "\n P(x) = ";
    
    long double EPS = 1e-12;
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
    export_values(x, coeff);

    return 0;
}