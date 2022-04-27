//
// Created by deafmist on 14.09.2021.
//

#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
using namespace std;

const vector<int> x = {1, 2, 3, 4, 5};
const vector<int> y = {-1, -2, -4, 0, -5};

double function1(double a, double b) {
    double myfunc = 0;
    for (int i = 0; i < 5; ++i) {
        myfunc += pow(a * x[i] + b - y[i], 2);
    }

    return myfunc;
}

double function2(double a, double b) {
    double myfunc = 0;
    for (int i = 0; i < 5; ++i) {
        myfunc += abs(a * x[i] + b - y[i]);
    }

    return myfunc;
}

vector<double> function(double gamma) {
    srand(time(NULL));

    double a0 = 1 + rand() % 10;
    double b0 = 1 + rand() % 10;
    double g1 = (double)rand() / (double) RAND_MAX;
    double g2 = 1 - pow(g1, 2);

    return {function1(a0 + gamma * g1, b0 + gamma * g2), function2(a0 + gamma * g1, b0 + gamma * g2)};
}

void methodDicho(double a0, double b0, double epsi) {
    double delta = epsi;
    double a = a0;
    double b = b0;
    double ek = (b - a) / 2;
    double c, d;
    vector<vector<double>> arr1;
    vector<vector<double>> arr2;
    double n = 0;

    while (ek > epsi) {
        c = (a + b - delta) / 2;
        d = (a + b + delta) / 2;
        arr1.push_back(vector<double>{n, a, b, ek, c, d, function(c)[0], function(d)[0]});

        if (function(c)[0] <= function(d)[0]) {
            b = d;
        }
        else {
            a = c;
        }
        ek = (b - a) / 2;
        n++;
    }

    a = a0;
    b = b0;
    ek = (b - a) / 2;
    n = 0;
    while (ek > epsi) {
        c = (a + b - delta) / 2;
        d = (a + b + delta) / 2;

        arr2.push_back(vector<double> {n, a, b, ek, c, d, function(c)[1], function(d)[1]});

        if (function(c)[1] <= function(d)[1]) {
            b = d;
        }
        else {
            a = c;
        }
        ek = (b - a) / 2;
        n++;
    }

    cout << "--------------------Dichotomy method, function #1--------------------" << endl;
    cout << "n" << "\t\t" << "an" << "\t\t" << "bn" << "\t\t" << "en" << "\t\t"
         << "c" << "\t\t" << "d" << "\t\t" << "phi(c)" << "\t\t" << "phi(d)" << endl;
    for (int i = 0; i < arr1.size(); ++i) {
        cout << double(int(arr1[i][0] * 10000)) / 10000 << "\t\t" << double(int(arr1[i][1] * 10000)) / 10000 << "\t\t"
             <<  double(int(arr1[i][2] * 10000)) / 10000 << "\t\t" << double(int(arr1[i][3] * 10000)) / 10000 << "\t\t"
             << double(int(arr1[i][4] * 10000)) / 10000 << "\t\t" << double(int(arr1[i][5] * 10000)) / 10000 << "\t\t"
             << double(int(arr1[i][6] * 10000)) / 10000 << "\t\t" << double(int(arr1[i][7] * 10000)) / 10000 << endl;
    }

    cout << endl;
    cout << "--------------------Dichotomy method, function #2--------------------" << endl;
    cout << "n" << "\t\t" << "an" << "\t\t" << "bn" << "\t\t" << "en" << "\t\t"
         << "c" << "\t\t" << "d" << "\t\t" << "phi(c)" << "\t\t" << "phi(d)" << endl;
    for (int i = 0; i < arr2.size(); ++i) {
        cout << double(int(arr2[i][0] * 10000)) / 10000 << "\t\t" << double(int(arr2[i][1] * 10000)) / 10000 << "\t\t"
             <<  double(int(arr2[i][2] * 10000)) / 10000 << "\t\t" << double(int(arr2[i][3] * 10000)) / 10000 << "\t\t"
             << double(int(arr2[i][4] * 10000)) / 10000 << "\t\t" << double(int(arr2[i][5] * 10000)) / 10000 << "\t\t"
             << double(int(arr2[i][6] * 10000)) / 10000 << "\t\t" << double(int(arr2[i][7] * 10000)) / 10000 << endl;
    }
}

void methodGold(double a0, double b0, double epsi) {
    double a = a0;
    double b = b0;
    double k = 0;
    double ek = pow((sqrt(5) - 1) / 2, k) * (b - a);
    double c = a + (3 - sqrt(5)) * (b - a) / 2;
    double d = a + (sqrt(5) - 1) * (b - a) / 2;
    vector<vector<double>> arr1;
    vector<vector<double>> arr2;

    while (ek > epsi) {
        arr1.push_back(vector<double> {k, ek, a, b, c, d, function(c)[0], function(d)[0]});

        if (function(c)[0] <= function(d)[0]) {
            b = d;
            d = c;
            c = a + (3 - sqrt(5)) * (b - a) / 2;
        }
        else {
            a = c;
            c = d;
            d = a + (sqrt(5) - 1) * (b - a) / 2;
        }
        k++;
        ek = pow((sqrt(5) - 1) / 2, k) * (b - a);
    }

    a = a0;
    b = b0;
    k = 0;
    ek = pow((sqrt(5) - 1) / 2, k) * (b - a);
    c = a + (3 - sqrt(5)) * (b - a) / 2;
    d = a + (sqrt(5) - 1) * (b - a) / 2;
    while (ek > epsi) {
        arr2.push_back(vector<double> {k, ek, a, b, c, d, function(c)[1], function(d)[1]});

        if (function(c)[1] <= function(d)[1]) {
            b = d;
            d = c;
            c = a + (3 - sqrt(5)) * (b - a) / 2;
        }
        else {
            a = c;
            c = d;
            d = a + (sqrt(5) - 1) * (b - a) / 2;
        }
        k++;
        ek = pow((sqrt(5) - 1) / 2, k) * (b - a);
    }

    cout << endl << endl;
    cout << "--------------------Golden ratio method, function #1--------------------" << endl;
    cout << "n" << "\t\t" << "en" << "\t\t" << "an" << "\t\t" << "bn" << "\t\t"
         << "c" << "\t\t" << "d" << "\t\t" << "phi(c)" << "\t\t" << "phi(d)" << endl;
    for (int i = 0; i < arr1.size(); ++i) {
        cout << double(int(arr1[i][0] * 10000)) / 10000 << "\t\t" << double(int(arr1[i][1] * 10000)) / 10000 << "\t\t"
             <<  double(int(arr1[i][2] * 10000)) / 10000 << "\t\t" << double(int(arr1[i][3] * 10000)) / 10000 << "\t\t"
             << double(int(arr1[i][4] * 10000)) / 10000 << "\t\t" << double(int(arr1[i][5] * 10000)) / 10000 << "\t\t"
             << double(int(arr1[i][6] * 10000)) / 10000 << "\t\t" << double(int(arr1[i][7] * 10000)) / 10000 << endl;
    }

    cout << endl;
    cout << "--------------------Golden ratio method, function #2--------------------" << endl;
    cout << "n" << "\t\t" << "en" << "\t\t" << "an" << "\t\t" << "bn" << "\t\t"
         << "c" << "\t\t" << "d" << "\t\t" << "phi(c)" << "\t\t" << "phi(d)" << endl;
    for (int i = 0; i < arr2.size(); ++i) {
        cout << double(int(arr2[i][0] * 10000)) / 10000 << "\t\t" << double(int(arr2[i][1] * 10000)) / 10000 << "\t\t"
             <<  double(int(arr2[i][2] * 10000)) / 10000 << "\t\t" << double(int(arr2[i][3] * 10000)) / 10000 << "\t\t"
             << double(int(arr2[i][4] * 10000)) / 10000 << "\t\t" << double(int(arr2[i][5] * 10000)) / 10000 << "\t\t"
             << double(int(arr2[i][6] * 10000)) / 10000 << "\t\t" << double(int(arr2[i][7] * 10000)) / 10000 << endl;
    }
}

int main() {
    methodDicho(-10, 10, 10e-6);
    methodGold(-10, 10, 10e-6);

    return 0;
}
