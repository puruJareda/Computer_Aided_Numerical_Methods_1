#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;

// This function reads raw numbers line-by-line from your input files.
vector<long double> readRootsFromCSV(const string& filename) {
    vector<long double> roots;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "[ERROR] Could not open " << filename << " for reading.\n";
        return roots;
    }

    string line;
    // We start reading from the very first line
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        stringstream ss(line);
        long double val;
        // Extracts the first numeric value found on the line
        if (ss >> val) {
            roots.push_back(val);
        }
    }

    file.close();
    return roots;
}

// Successive error: |x_{n+1} - x_n|
vector<long double> successiveErrors(const vector<long double>& roots) {
    vector<long double> errors;
    for (int i = 0; i + 1 < (int)roots.size(); i++)
        errors.push_back(fabsl(roots[i + 1] - roots[i]));
    return errors;
}

// Convergence order: p_n = log(e_{n+1}) / log(e_n)
vector<long double> convergenceOrders(const vector<long double>& errors) {
    vector<long double> orders;
    for (int i = 0; i + 1 < (int)errors.size(); i++) {
        long double e_curr = errors[i];
        long double e_next = errors[i + 1];
        // Ensure values are above machine noise (1e-19 for long double)
        if (e_curr > 1e-19L && e_next > 1e-19L)
            orders.push_back(logl(e_next) / logl(e_curr));
        else
            orders.push_back(0.0L); // Limit of precision reached
    }
    return orders;
}

// Print table to console
void printTable(const string& label,
                const vector<long double>& roots,
                const vector<long double>& errors,
                const vector<long double>& orders) {
    cout << "\n" << "\n";
    cout << "  " << label << "\n";
    cout << "\n";
    cout << setw(6)  << "Iter"
         << setw(20) << "Root Estimate"
         << setw(20) << "Succ. Error"
         << setw(16) << "Conv. Order" << "\n\n";

    if (roots.empty()) return;

    // Row 1 (Initial guess/first iteration)
    cout << setw(6)  << 1
         << setw(20) << fixed << setprecision(12) << roots[0]
         << setw(20) << "---"
         << setw(16) << "---" << "\n";

    for (int i = 0; i < (int)errors.size(); i++) {
        cout << setw(6)  << i + 2
             << setw(20) << fixed << setprecision(12) << roots[i + 1]
             << setw(20) << scientific << setprecision(6) << errors[i];

        if (i < (int)orders.size() && orders[i] != 0.0L)
            cout << setw(16) << fixed << setprecision(6) << orders[i];
        else
            cout << setw(16) << "---";

        cout << "\n";
    }
    cout << "\n";
}

// Export Function: Individual File Export
void exportIndividualCSV(const string& filename, 
                         const string& methodName,
                         const vector<long double>& errors, 
                         const vector<long double>& orders) {
    
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "[ERROR] Could not write to " << filename << "\n";
        return;
    }

    // Header specific to this method
    file << "iteration,error,order\n";

    for (int i = 0; i < (int)errors.size(); i++) {
        file << i + 1 << "," 
             << scientific << setprecision(15) << errors[i] << ",";
        
        if (i < (int)orders.size())
            file << fixed << setprecision(6) << orders[i] << "\n";
        else
            file << "0.000000\n"; // Last row where order can't be computed
    }

    file.close();
}

int main() {
    
    string secantFile = "secant_iterations.csv";
    string newtonFile = "newton_iterations.csv";

    // 1. Data Import
    vector<long double> sec_roots = readRootsFromCSV(secantFile);
    vector<long double> nwt_roots = readRootsFromCSV(newtonFile);

    // 2. Safety Check
    if (sec_roots.empty() || nwt_roots.empty()) {
        cerr << "[!] Error: Ensure your root files exist and contain data.\n";
        return 1;
    }

    // 3. Analysis Calculation
    vector<long double> sec_err = successiveErrors(sec_roots);
    vector<long double> sec_ord = convergenceOrders(sec_err);
    
    vector<long double> nwt_err = successiveErrors(nwt_roots);
    vector<long double> nwt_ord = convergenceOrders(nwt_err);

    // 4. Console Presentation
    printTable("SECANT METHOD ANALYSIS", sec_roots, sec_err, sec_ord);
    printTable("NEWTON-RAPHSON ANALYSIS", nwt_roots, nwt_err, nwt_ord);

    // 5. Final Export
    exportIndividualCSV("secant_analysis.csv", "Secant", sec_err, sec_ord);
    exportIndividualCSV("newton_analysis.csv", "Newton-Raphson", nwt_err, nwt_ord);

    return 0;
}
