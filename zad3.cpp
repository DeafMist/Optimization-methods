#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>
using namespace std;

const vector<int> x = {1, 2, 3, 4, 5};
const vector<int> y = {-1, -2, -4, 0, -5};

double x_prev;
double y_prev;
double x_cur;
double y_cur;
vector<double> grad;

double function1(double a, double b) {
    double myfunc = 0;

    for (int i = 0; i < 5; ++i) {
        myfunc += pow(a * x[i] + b - y[i], 2);
    }

    return myfunc;
}

double funcAlfaX1(double alfa) {
    return function1(x_cur + alfa, y_cur);
}

double funcAlfaY1(double alfa) {
    return function1(x_cur, y_cur + alfa);
}

double diffX(double a, double b) {
    double delta = 10e-9;
    return (function1(a + delta, b) - function1(a, b)) / delta;
}

double diffY(double a, double b) {
    double delta = 10e-9;
    return (function1(a, b + delta) - function1(a, b)) / delta;
}

double function2(double a, double b) {
    double myfunc = 0;

    for (int i = 0; i < 5; ++i) {
        myfunc += abs(a * x[i] + b - y[i]);
    }

    return myfunc;
}

double funcAlfaX2(double alfa) {
    return function2(x_cur + alfa, y_cur);
}

double funcAlfaY2(double alfa) {
    return function2(x_cur, y_cur + alfa);
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

void coordDescentMethod(double (*f)(double, double), double (*fx)(double), double (*fy)(double), double epsi) {
    double alfa = 100;
    double gamma;
    int n = 1;

    cout << "n" << "\t\t" << "xn" << "\t\t" << "||xn - xn-1||" << "\t\t" << "f(xn)" << "\t\t" << "|f(xn) - f(xn-1)|" << endl;

    while (sqrt(pow(x_cur - x_prev, 2) + pow(y_cur - y_prev, 2)) > epsi && abs(f(x_cur, y_cur) - f(x_prev, y_prev)) > epsi) {
        x_prev = x_cur;
        y_prev = y_cur;
        gamma = methodDicho(fx, 0, alfa, epsi / 10);

        if (gamma < epsi)
            gamma = methodDicho(fx, -alfa, 0, epsi / 10);

        x_cur += gamma;

        gamma = methodDicho(fy, 0, alfa, epsi / 10);

        if (gamma < epsi)
            gamma = methodDicho(fy, -alfa, 0, epsi / 10);


        y_cur += gamma;

        cout << n << "\t\t" << '{' << double(int(x_cur * 10000)) / 10000 << ", " << double(int(y_cur * 10000)) / 10000
             << '}' << "\t\t" << double(int(sqrt(pow(x_cur - x_prev, 2) + pow(y_cur - y_prev, 2))) * 10000) / 10000
             << "\t\t" << double(int(f(x_cur, y_cur) * 10000)) / 10000
             << "\t\t" << double(int(abs(f(x_cur, y_cur) - f(x_prev, y_prev)) * 10000)) / 10000 << endl;

        n++;
    }
}


vector<double> gradient(double a, double b) {
    return {diffX(a, b), diffY(a, b)};
}

double funcAlfaXY(double alfa) {
    return function1(x_prev - alfa * grad[0], y_prev - alfa * grad[1]);
}

void gradSplitStep(double epsi) {
    double delta = (double)rand() / (double) RAND_MAX;
    double alfa = 1000;
    double a = -5000 + rand() % 10000;
    double b = -5000 + rand() % 10000;
    grad = gradient(a, b);
    int n = 0;

    cout << "n" << "\t\t" << "xn" << "\t\t\t" << "f(xn)" << "\t\t\t" << "||f'(xn)||" << endl;

    while (sqrt(pow(grad[0], 2) + pow(grad[1], 2)) > epsi) {
        grad = gradient(a, b);

        cout << n << "\t\t" << '{' << double(int(a * 10000)) / 10000 << ", " << double(int(b * 10000)) / 10000
             << '}' << "\t\t" << double(int(function1(a, b) * 10000)) / 10000
             << "\t\t" << double(int(sqrt(pow(grad[0], 2) + pow(grad[1], 2)) * 10000)) / 10000 << endl;

        while (function1(a - alfa * grad[0], b - alfa * grad[1]) - function1(a, b) >
                                -alfa * delta * pow(sqrt(pow(grad[0], 2) + pow(grad[1], 2)), 2))
            alfa /= 2;

        a -= alfa * grad[0];
        b -= alfa * grad[1];

        n++;
    }
}

void gradConstStep(double epsi) {
    double delta = (double)rand() / (double) RAND_MAX;
    double alfa = 2 * (1 - delta) / 63;
    double a = -5000 + rand() % 10000;
    double b = -5000 + rand() % 10000;
    grad = gradient(a, b);
    int n = 0;

    cout << "n" << "\t\t" << "xn" << "\t\t\t" << "f(xn)" << "\t\t\t" << "||f'(xn)||" << endl;

    while (sqrt(pow(grad[0], 2) + pow(grad[1], 2)) > epsi) {
        grad = gradient(a, b);

        cout << n << "\t\t" << '{' << double(int(a * 10000)) / 10000 << ", " << double(int(b * 10000)) / 10000
             << '}' << "\t\t" << double(int(function1(a, b) * 10000)) / 10000
             << "\t\t" << double(int(sqrt(pow(grad[0], 2) + pow(grad[1], 2)) * 10000)) / 10000 << endl;

        a -= alfa * grad[0];
        b -= alfa * grad[1];
        n++;
    }
}

void gradSteepestDescent(double epsi) {
    double alfa;
    grad = gradient(x_prev, y_prev);
    int n = 0;

    cout << "n" << "\t\t" << "xn" << "\t\t\t" << "f(xn)" << "\t\t\t" << "||f'(xn)||" << endl;

    while (sqrt(pow(grad[0], 2) + pow(grad[1], 2)) > epsi) {
        grad = gradient(x_prev, y_prev);

        cout << n << "\t\t" << '{' << double(int(x_prev * 10000)) / 10000 << ", " << double(int(y_prev * 10000)) / 10000
             << '}' << "\t\t" << double(int(function1(x_prev, y_prev) * 10000)) / 10000
             << "\t\t" << double(int(sqrt(pow(grad[0], 2) + pow(grad[1], 2)) * 10000)) / 10000 << endl;

        alfa = methodDicho(funcAlfaXY, 0, 100, epsi);

        x_prev -= alfa * grad[0];
        y_prev -= alfa * grad[1];
        n++;
    }
}

int main() {
    srand(time(NULL));
    x_prev = -5000 + rand() % 10000;
    y_prev = -5000 + rand() % 10000;
    x_cur = x_prev + 10;
    y_cur = y_prev + 10;

    cout << "--------------------Coord descent method, function #1--------------------" << endl;
    coordDescentMethod(function1, funcAlfaX1, funcAlfaY1, 10e-6);
    cout << "\n\n";

    cout << "--------------------Grad with split step, function #1--------------------" << endl;
    gradSplitStep(10e-6);
    cout << "\n\n";

    cout << "--------------------Grad with const step, function #1--------------------" << endl;
    gradConstStep(10e-6);
    cout << "\n\n";

    x_prev = -5000 + rand() % 10000;
    y_prev = -5000 + rand() % 10000;

    cout << "--------------------Grad steepest descent, function #1--------------------" << endl;
    gradSteepestDescent(10e-6);
    cout << "\n\n";

    x_prev = -5000 + rand() % 10000;
    y_prev = -5000 + rand() % 10000;
    x_cur = x_prev + 10;
    y_cur = y_prev + 10;

    cout << "--------------------Coord descent method, function #2--------------------" << endl;
    coordDescentMethod(function2, funcAlfaX2, funcAlfaY2, 10e-6);
    cout << "\n\n";
}