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

    
    void Gradient_Descent(vector<double>& x,vector<double>& y,vector<double>& c,int d,double alpha,int max_iter,double tol){
        int N = x.size();

        vector<vector<double>> A(N, vector<double>(d+1, 0));
        for(int i=0;i<N;i++){
            for(int k=0;k<d+1;k++){
                A[i][k]=T_x_k(k,x[i]);
            }
        }
        c.assign(d + 1, 0.0);

        cout << fixed << setprecision(8);
        cout << "\n--- Gradient Descent Log ---\n";
        cout << setw(8)  << "Iter"<< setw(18) << "Loss L(c)"<< setw(18) << "||gradient||" << "\n";
        cout << string(44, '-') << "\n";

        for (int iter = 0; iter <= max_iter; iter++) {

            vector<double> residual(N,0);

            for(int i=0; i<N; i++) {
                double s = 0.0;
                for(int j=0; j<d+1; j++){
                    s += c[j]*A[i][j];
                residual[i] = s - y[i];
                }
            }
            double loss = 0;
            for (int i= 0; i<N; i++){
                loss += residual[i]*residual[i];
            }
            loss /= (2*N);

            vector<double>grad(d+1,0);

            for(int i=0; i<d+1; i++){
                for(int j=0; j<N; j++){
                    grad[i] += A[j][i]*residual[j];
                }
            grad[i] /= N;
            }

            double grad_norm = 0.0;
            for (int k = 0; k <= d; k++){
                grad_norm += grad[k] * grad[k];
            }
            grad_norm = sqrt(grad_norm);

            if (iter % 100 == 0)
                cout << setw(8)  << iter << setw(18) << loss << setw(18) << grad_norm << "\n";

            if (grad_norm<tol) {
                cout <<"\nConverged at iteration "<< iter<<"  (||grad|| = "<<grad_norm<<" < tol = "<<tol<< ")\n";
                break;
            }

            for (int k=0;k<=d;k++){
                c[k] -= alpha*grad[k];
            }
        }

        double max_err=0, rms_err=0;
        for (int i=0; i<N; i++) {
            double s = 0.0;
            for (int j=0; j<d+1; j++){
                s += c[j] * A[i][j];
            }
            double e = fabs(s-y[i]);
            max_err  = max(max_err, e);
            rms_err += e * e;
        }
        rms_err = sqrt(rms_err / N);
        
        cout<<endl;
        cout << "Max error : " << max_err  << "\n";
        cout << "RMS error : " << rms_err  << "\n";

    }   
    int main(){
    vector<double> x, y;
    read_data(x, y);

    int d=9;       
    double alpha=0.01;
    int max_iter=10000; 
    double tol=1e-7;  

    vector<double> c;
    Gradient_Descent(x,y,c,d,alpha,max_iter,tol);

    cout<<"\nChebyshev Coefficients:\n";
    for (int i=0; i<=d; i++)
        cout<<"c[" << i << "] = "<<c[i]<<"\n";

    print_standard_polynomial(c);
    return 0;
    }


