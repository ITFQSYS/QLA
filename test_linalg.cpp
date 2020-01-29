#include "qla.hpp"

#include <iostream>
#include <cmath>
using namespace QLA;
using namespace std;

bool test_SOR_solver()
{
    Matrix
        a = {10., 3., 1., 2., 1.,
             1., 19., 2., -1., 5.,
             -1., 1., 30., 1., 10.,
             -2., 0., 1., 20., 5.,
             -3., 5., 1., -2., 25.};
    a.reshape(5, 5);
    //cout << a << endl;

    Matrix b = {-22, 27, 89, -73, 22};
    b.reshape(5, 1);
    Matrix x(5, 1);

    SOR_solver(a, x, b, 1.1, 1e-10);

    cout << "Error:" << sqrt((a.dot(x) - b).T().dot((a.dot(x) - b))(0, 0)) << endl;
    return true;
}

int main(int argc, char const *argv[])
{
    test_SOR_solver();
    return 0;
}
