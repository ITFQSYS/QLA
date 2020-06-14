#ifndef __QLA_LINALG_HPP__
#define __QLA_LINALG_HPP__ 1

#include "Types.hpp"
#include <cassert>
#include <iostream>
#include <limits>

namespace QLA
{
bool SOR_solver(Matrix const &a, Matrix &x, Matrix const &b, double omega, double eps = 1e-8)
{
    const int n = a.getCols();         //ベクトルbの求められるサイズ
    assert(1 < omega && omega < 2);    //SOR法の条件
    assert(a.getRows() == n && n > 0); //行列Aが有効な正方行列
    assert(x.getRows() == n);          //axが成立するか
    assert(x.getCols() == 1);          //xは行ベクトル
    assert(b.getCols() == 1);          //bは行ベクトル
    assert(b.getRows() == n);          //連立方程式が成立しているか

    Matrix x_new(n, 1, 2);
    double err = std::numeric_limits<double>::max();//doubleの最大値
    int count_iter = 0;

    while (err > eps && count_iter < 1e9)
    {

        for (int i = 0; i < n; i++)
        {
            double lx = 0;
            double udx = 0;

            for (int j = 0; j < i; j++)
            {
                lx += a(i, j) * x_new(j);
            }
            for (int j = i; j < n; j++)
            {
                udx += a(i, j) * x(j);
            }

            x_new(i) += omega / a(i, i) * (b(i) - lx - udx);
        }
        Matrix dx = x_new - x;
        err = dx.dot(dx.T())(0, 0);

        for (int i = 0; i < n; i++)
        {
            x(i) = x_new(i);
        }
        count_iter++;

        std::cout << "iteration\t" << count_iter << "\terr=" << err << std::endl;
    }

    if (count_iter >= 1e9) //収束しない場合
    {
        std::cerr << "SOR method failed!" << std::endl;
        return false;
    }
    return true;
}
} // namespace QLA

#endif