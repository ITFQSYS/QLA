#ifndef __QLA_LINALG_HPP__
#define __QLA_LINALG_HPP__ 1

#include "Types.hpp"
#include <cassert>
#include <iostream>
#include <limits>
#include <omp.h>

namespace QLA
{
    //SOR法を用いて解く omega=1ならばGauss-Seidel法になる
    bool solve_sor(Matrix const &a, Matrix &x, Matrix const &b, double omega, double eps = 1e-8, bool progress = false)
    {
        const int n = a.getCols();         //ベクトルbの求められるサイズ
        assert(1 < omega && omega < 2);    //SOR法の条件
        assert(a.getRows() == n && n > 0); //行列Aが有効な正方行列
        assert(x.getRows() == n);          //axが成立するか
        assert(x.getCols() == 1);          //xは行ベクトル
        assert(b.getCols() == 1);          //bは行ベクトル
        assert(b.getRows() == n);          //連立方程式が成立しているか

        Matrix x_new(n, 1, 2);
        double err = std::numeric_limits<double>::max(); //doubleの最大値
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

            if (progress)
            {
                std::cout << "iteration\t" << count_iter << "\terr=" << err << std::endl;
            }
        }

        if (count_iter >= 1e9) //収束しない場合
        {
            std::cerr << "SOR method failed!" << std::endl;
            return false;
        }
        return true;
    }

    //LU分解を行う ピボット探索の結果は置換行列Pに出力される
    //PA=LU の関係がある.
    // 参考: http://www.ritsumei.ac.jp/se/~hirai/edu/2019/algorithm/handout/linear_8up.pdf
    void lu(Matrix a, Matrix &l, Matrix &u, Matrix &p)
    {
        assert(a.getRows() == a.getCols()); //Aが正方行列
        int n = a.getCols();                //行列のサイズ

        //Lのサイズが正しいか?
        assert(l.getRows() == l.getCols());
        assert(l.getRows() == n);
        //Uのサイズが正しいか?
        assert(u.getRows() == u.getCols());
        assert(u.getRows() == n);
        //Pのサイズが正しいか?
        assert(p.getRows() == p.getCols());
        assert(p.getRows() == n);

        std::vector<int> index(n); //置換の記録

        //変数の初期化
        //#pragma omp parallel for
        for (int i = 0; i < n * n; i++)
        {
            l(i) = u(i) = p(i) = 0;
        }
        //#pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            l(i, i) = 1;
            index[i] = i;
        }

        for (int k = 0; k < n; k++)
        {
            int p = k;
            double a_max = abs(a(k, k));

            //ピボットを探す
            for (int j = k + 1; j < n; j++)
            {
                if (a_max < abs(a(k, j)))
                {
                    p = j;
                    a_max = abs(a(k, j));
                }
            }

            //Aの行の入れ替え
            for (int j = k; j < n; j++)
            {
                std::swap(a(k, j), a(p, j));
            }
            //Lの行の入れ替え
            for (int j = 0; j < k; j++)
            {
                std::swap(l(k, j), l(p, j));
            }

            std::swap(index[k], index[p]); //置換の記録

            for (int j = k; j < n; j++)
            {
                u(k, j) = a(k, j);
            }

            for (int i = k + 1; i < n; i++)
            {
                l(i, k) = a(i, k) / u(k, k);
            }

            for (int i = k + 1; i < n; i++)
            {
                for (int j = k + 1; j < n; j++)
                {
                    a(i, j) -= l(i, k) * u(k, j);
                }
            }
        }
        //Pの復元
        for (int k = 0; k < n; k++)
        {
            p(k, index[k]) = 1;
        }
    }

    //連立方程式Ax=bを直接解く
    bool solve(const Matrix &a, Matrix &x, const Matrix &b)
    {
        //連立方程式の条件チェック
        assert(a.getCols() == a.getRows());
        assert(a.getCols() == b.getRows());
        assert(b.getCols() == 1);

        int n = a.getCols(); //行列のサイズ

        if (x.getCols() != 1 || x.getRows() != n) //xのサイズ調整
        {
            x = Matrix(n, 1);
        }

        Matrix l(n, n), u(n, n), p(n, n);

        lu(a, l, u, p); //LU分解

        Matrix pb = p.dot(b);

        //Ax=bの両辺にpを掛けてpA=Pbを解く.
        //pA=LUだから LU=pbを解く

        Matrix y(n, 1, 0);

        //Ly=_bを解く
        for (int i = 0; i < n; i++)
        {
            double sum = 0;
#pragma omp parallel for
            for (int j = 0; j <= i - 1; j++)
            {
                sum += l(i, j) * y(j);
            }
            y(i) = pb(i) - sum;
        }

        //Ux=yを解く
        for (int i = n - 1; i >= 0; --i)
        {
            double sum = 0;

#pragma omp parallel for
            for (int j = i + 1; j < n; j++)
            {
                sum += u(i, j) * x(j);
            }
            x(i) = (y(i) - sum) / u(i, i);
        }
        return true;
    }

} // namespace QLA

#endif