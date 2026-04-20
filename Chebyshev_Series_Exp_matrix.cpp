#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <functional>
#include <string>
#include<fstream>

using namespace std;

void print_standard_polynomial(const vector<double>& c) {
    int d = c.size() - 1;

    // T[k] holds the standard coefficients of T_k(x)
    // T[k][i] = coefficient of x^i in T_k(x)
    vector<vector<double>> T(d + 1, vector<double>(d + 1, 0));

    T[0][0] = 1;
    if (d >= 1)
        T[1][1] = 1;

    for (int k=2; k<d+1; k++) {
        for (int i=1; i<=k; i++)
            T[k][i] += 2.0*T[k-1][i-1];
        for (int i=0; i<=k; i++)
            T[k][i] -= T[k-2][i];
    }

    vector<double> poly(d + 1, 0.0);
    for (int k = 0; k <= d; k++)
        for (int i = 0; i <= d; i++)
            poly[i] += c[k] * T[k][i];

    cout << "\nInterpolating Polynomial (Standard Form):\n";
    cout << "p(x) = ";

    int first_term = 1;
    for (int i = d; i >= 0; i--) {
        if (abs(poly[i]) < 1e-10) continue;
        double coeff = poly[i];
        if (first_term) {
            cout << fixed << setprecision(6) << coeff;
        } else {
            if (coeff >= 0)
                cout << " + " << fixed << setprecision(6) << coeff;
            else
                cout << " - " << fixed << setprecision(6) << abs(coeff);
        }
        if (i == 1)
            cout << "x";
        else if (i > 1)
            cout << "x^" << i;

        first_term = 0;
    }
    if (first_term) cout << "0";
    cout << "\n" << endl;
}

void read_data(vector<double>& x, vector<double>& y){
    ifstream file("group2.csv");

    double x_val, y_val;
    char comma;

    file.ignore(1000, '\n');

    while(file >> x_val >> comma >> y_val) {
        x.emplace_back(x_val);
        y.emplace_back(y_val);
    }
  
    file.close();
}

double T_x_k(int k,double x){
        if(k==0){
            return 1;
        }
        else if(k==1){
            return x;
        }
        double T_n,T_n_1,T_n_2;
        T_n_2=1;
        T_n_1=x;
        for(int i=2;i<=k;i++){
            T_n = 2*x*T_n_1 - T_n_2;
            T_n_2 = T_n_1;
            T_n_1 = T_n;
        }
        return T_n;
    }


void Chebyshev_Polynomial_Approximation(vector<double>& x,vector<double>& y,vector<double>& c,int d){
    int N = x.size();
    vector<vector<double>> A(N, vector<double>(d+1, 0));
    for(int i=0;i<N;i++){
        for(int k=0;k<d+1;k++){
            A[i][k]=T_x_k(k,x[i]);
        }
    }

    vector<vector<double>>AtA(d+1, vector<double>(d+1, 0));
    for(int k=0;k<d+1;k++){
        for(int j=0;j<d+1;j++){
            for(int i=0;i<N;i++){
                AtA[k][j] += A[i][j]*A[i][k];
            }
        }
    }

    vector<vector<double>>Aty(d+1,vector<double>(1,0));
    for(int i=0;i<d+1;i++){
        for(int j=0;j<N;j++){
            Aty[i][0] += A[j][i]*y[j];
        }
    }

    int n = d + 1;
    vector<vector<double>> M(n, vector<double>(n+1));
    for (int i =0; i<n; i++) {
        for (int j=0; j<n; j++)
            M[i][j] = AtA[i][j];
        M[i][n] = Aty[i][0];
    }

    for (int col=0; col<n; col++){
        int pivot = col;
        for (int row=col+1; row<n; row++)
            if (abs(M[row][col]) > abs(M[pivot][col]))
                pivot = row;
        swap(M[col], M[pivot]);

        for (int row=col+1; row<n; row++) {
            double factor = M[row][col]/M[col][col];
            for (int j=col; j<=n; j++)
                M[row][j] -= factor*M[col][j];
        }
        for (int i=n-1; i>=0; i--) {
            c[i] = M[i][n];
            for (int j=i+1; j<n; j++)
                c[i] -= M[i][j]*c[j];
                c[i] /= M[i][i];
        }
    }
}
int main(){
    
    int N=25;
    int d=9;
    vector<double>x,y;
    vector<double>c(d+1,0);

    read_data(x,y);

    Chebyshev_Polynomial_Approximation(x,y,c,d);

    for(int i=0;i<d+1;i++){
         cout << "c[" << i << "] = " << c[i] << endl;
    }

    print_standard_polynomial(c);

}