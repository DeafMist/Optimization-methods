#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
#include <map>
using namespace std;

const vector<int> x = {1, 2, 3, 4, 5};
const vector<int> y = {-1, -2, -4, 0, -5};

double a0;
double b0;
double g1;
double g2;


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

double function3(double a, double b) {
    double max = 0;
    double cur = 0;
    for (int i = 0; i < 5; ++i) {
        cur = abs(a * x[i] + b - y[i]);
        if (cur > max)
            max = cur;
    }

    return max;
}

vector<double> function(double gamma) {
    return {function1(a0 + gamma * g1, b0 + gamma * g2),
            function2(a0 + gamma * g1, b0 + gamma * g2),
            function3(a0 + gamma * g1, b0 + gamma * g2)};
}

double wrap1(double gamma) {
    return function(gamma)[0];
}

double wrap2(double gamma) {
    return function(gamma)[1];
}

double wrap3(double gamma) {
    return function(gamma)[2];
}

double diff(double gamma, double (*f)(double)) {
    double dx = 10e-6;
    return -abs((f(gamma + dx) - f(gamma)) / dx);
}

double wrapDiff1(double gamma) {
    return diff(gamma, wrap1);
}

double wrapDiff2(double gamma) {
    return diff(gamma, wrap2);
}

double wrapDiff3(double gamma) {
    return diff(gamma, wrap3);
}

double methodDicho(double (*f)(double), double a, double b, double epsi) {
    double delta = epsi;
    double ek = (b - a) / 2;
    double c, d;

    while (ek > epsi) {
        c = (a + b - delta) / 2;
        d = (a + b + delta) / 2;

        if (f(c) <= f(d)) {
            b = d;
        } else {
            a = c;
        }
        ek = (b - a) / 2;
    }

    return (a + b) / 2;
}

void polyLineMethod(double (*f)(double), double a, double b, double epsi, double L) {
    map <double, double> dots;
    map <double, double> :: iterator it;

    double xk = (f(a) - f(b) + L * (a + b)) / (2 * L);
    double pk = (f(a) + f(b) + L * (a - b)) / 2;
    dots[xk] = pk;
    double delta = (f(xk) - pk) / (2 * L);
    double x1, x2, lastPk, min, minIndex;
    int n = 2;

    cout << "n" << "\t\t" << "xn*" << "\t\t" << "pn*" << "\t\t" << "2*L*delta" << "\t\t" << "xn1" << "\t\t"
         << "xn2" << "\t\t" << "pn" << endl;

    while (2 * L * delta > epsi) {
        delta = (f(xk) - pk) / (2 * L);
        x1 = xk - delta;
        x2 = xk + delta;
        
        pk = (f(xk) + pk) / 2;
        dots[x1] = pk;
        dots[x2] = pk;
        lastPk = pk;
        
        it = dots.begin();
        min = it->second;
        minIndex = it->first;
        for (; it != dots.end(); ++it) {
            if (it->second < min) {
                min = it->second;
                minIndex = it->first;
            }
        }
        xk = minIndex;
        pk = min;

        if (n % 10 == 0) {
            cout << n << "\t\t" << double(int(xk * 10000)) / 10000 << "\t\t" << double(int(pk * 10000)) / 10000
                 << "\t\t" << double(int(2 * L * delta * 10000)) / 10000 << "\t\t"
                 << double(int(x1 * 10000)) / 10000 << "\t\t"
                 << double(int(x2 * 10000)) / 10000 << "\t\t" << double(int(lastPk * 10000)) / 10000 << endl;
        }

        n++;
        dots.erase(dots.find(xk));
    }
}

double differential(double gamma, double (*f)(double)) {
    double dx = 10e-6;
    return (f(gamma + dx) - f(gamma)) / dx;
}

double wrapDifferential1(double gamma) {
    return differential(gamma, wrap1);
}

void newtoneRaphson(double (*wrapF)(double), double a, double b, double epsi) {
    double df;
    double xk = a + rand() % int(b - a);
    int n = 0;
    bool flag = true;

    cout << "n" << "\t\t\t" << "xn" << "\t\t\t" << "f(xn)" << "\t\t\t" << "f'(xn)" << endl;

    while ((abs(wrapDifferential1(xk)) > epsi) && flag) {
        df = wrapDifferential1(xk);

        cout << n << "\t\t\t" << xk << "\t\t\t" << wrapF(xk) << "\t\t\t" << df << endl;

        if ((xk - df/wrapDifferential1(wrapDifferential1(xk)) <= -10) || (xk - df/wrapDifferential1(wrapDifferential1(xk)) >= 10))
            flag = false;
        else
            xk = xk - df/wrapDifferential1(wrapDifferential1(xk));
        n++;
    }
}

int main() {
    srand(time(NULL));
    a0 = 1 + rand() % 5;
    b0 = 1 + rand() % 5;
    g1 = (double)rand() / (double) RAND_MAX;
    g2 = 1 - pow(g1, 2);

    vector<double> L = {0, 0, 0};

    double Lx = methodDicho(wrapDiff1, -10, 10, 10e-6);
    L[0] = abs(wrapDiff1(Lx)) + 10e-5;

    Lx = methodDicho(wrapDiff2, -10, 10, 10e-6);
    L[1] = abs(wrapDiff2(Lx)) + 10e-5;

    Lx = methodDicho(wrapDiff3, -10, 10, 10e-6);
    L[2] = abs(wrapDiff3(Lx)) + 10e-5;

    cout << "--------------------Polyline method, function #1--------------------" << endl;
    polyLineMethod(wrap1, -10, 10, 10e-5, L[0]);
    cout << "\n\n";

    cout << "--------------------Newton-Raphson, function #1--------------------" << endl;
    newtoneRaphson(wrap1,-10, 10, 10e-6);
    cout << "\n\n";

    cout << "--------------------Polyline method, function #2--------------------" << endl;
    polyLineMethod(wrap2, -10, 10, 10e-5, L[1]);
    cout << "\n\n";

    cout << "--------------------Polyline method, function #3--------------------" << endl;
    polyLineMethod(wrap3, -10, 10, 10e-5, L[2]);
    cout << "\n\n";
}

