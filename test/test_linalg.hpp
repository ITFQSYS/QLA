#include "../include/qla.hpp"
#include <gtest/gtest.h>
using namespace QLA;
using namespace std;

TEST(LinalgTest, SOR_test)
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

    solve_sor(a, x, b, 1.1, 1e-12);

    EXPECT_NEAR((b - a.dot(x)).l2norm(), 0, 1e-5);
}

TEST(LinalgTest, LU_decomposition_test)
{
    Matrix a = {0, 10, 20, 2, 1, 10, 10, 0, 5};
    a.reshape(3, 3);

    Matrix l(3,3), u(3,3), p(3,3);
    lu(a, l, u, p);

    EXPECT_NEAR((p.dot(a) -l.dot(u)).l2norm(),0,1e-8);
}

TEST(LinalgTest, solve_test)
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

    solve(a,x,b);
    
    EXPECT_NEAR((b - a.dot(x)).l2norm(), 0, 1e-8);
}
