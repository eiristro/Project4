#include <iostream>
#include "armadillo"
#include <iomanip>
#include <cmath>

using namespace std;
using namespace arma;

void tridiag(double a, double b, double c, colvec y, colvec &u, int n);

int main()
{
    int n = 11, i;
    double deltax = 1.0/(n-1), deltat = 1.0/200;
    double alpha = double(deltat)/(deltax*deltax);
    double a, b, c;
    double maxT = 1.0, t;
    ofstream ofile;
    colvec u(n), v(n), v0(n), vv(n), us(n), r(n);

    // The steady-state solution. We solve v numerically and use u_s as initial condition
    // and to find u for plotting.
    for (i = 0; i<n; i++) {
        us(i) = 1 - i*deltax;
    }
    v0 = -us;


    // Forward Euler
    ofile.open("../Project4/fEuler.dat");

    v = v0;
    t = 0.0;
    while (t<maxT) {
        vv = v;
        for (i = 1; i<n-1; i++) {
            v(i) = alpha*vv(i+1) + (1 - 2*alpha)*vv(i) + alpha*vv(i-1);
        }
        v(0) = 0;
        v(n-1) = 0;

        // Finding u from v and us, then saving this for plotting
        u = v + us;
        for (int j = 0; j < n; j++) {
            ofile << u(j) << ' ';
        }
        ofile << t << ' ';
        ofile << endl;
        t += deltat;
    }
    ofile.close();

    // Backward Euler
    a = -alpha;
    b = 1 + 2*alpha;
    c = -alpha;

    v = v0;
    vv = v;

    ofile.open("../Project4/bEuler.dat");
    t = 0.0;
    while (t<maxT) {
        tridiag(a, b, c, vv, v, n-1);
        v(0) = 0;
        v(n-1) = 0;
        vv = v;

        // Saving
        u = v + us;
        for (int j = 0; j < n; j++) {
            ofile << u(j) << ' ';
        }
        ofile << t << ' ';
        ofile << endl;
        t += deltat;
    }
    ofile.close();

    // Crank-Nicolson
    ofile.open("../Project4/CN.dat");

    a = -alpha;
    b = 2 + 2*alpha;
    c = -alpha;

    v = v0;
    r.fill(0);

    t = 0.0;
    while (t<maxT) {
        // First find a new right hand side
        for (i = 1; i<n-1; i++) {
            r(i) = alpha*v(i+1) + (2 - 2*alpha)*v(i) + alpha*v(i-1);
        }
        r(0) = 0;
        r(n-1) = 0;

        // Then solve the tridiagonal problem
        tridiag(a, b, c, r, v, n-1);
        v(0) = 0;
        v(n-1) = 0;

        // Saving
        u = v + us;
        for (int j = 0; j < n; j++) {
            ofile << u(j) << ' ';
        }
        ofile << t << ' ';
        ofile << endl;

        t += deltat;
    }
    ofile.close();

    return 0;
}

void tridiag(double a, double b, double c, colvec uu, colvec &u, int n) {
    // A function implementing the tridiagonal solver from the first project.
    // It takes u by reference so it doesn't need to return anything

    // Forwards substitution,  starting at the second element
    double temp = double(a) / b;
    double bb = b - temp*c;
    for (int i = 1; i < n+1; i++)
    {
        uu(i) -= temp*uu(i-1);
    }

    //Backwards substitution
    for (int i = n-1; i>0; i--)
    {
        u(i) = (uu(i) - c*u(i+1))/bb;
    }
}
