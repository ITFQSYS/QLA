#ifndef __QLA_MATRIX_HPP__
#define __QLA_MATRIX_HPP__ 1

#include <iostream>
#include <vector>
#include <cassert>
#include <utility>
#include <cmath>
#include "Types.hpp"

namespace QLA
{

/*コンストラクタ*/
inline Matrix::Matrix()
{
    rows = cols = 0;
}

/*
コンストラクタ
rows:行の数
cols:列の数
*/
inline Matrix::Matrix(unsigned int rows, unsigned int cols)
{
    data = new double[rows * cols];
    this->rows = rows;
    this->cols = cols;
}

/*
コンストラクタ
rows:行の数
cols:列の数
v:初期化する値
*/
inline Matrix::Matrix(unsigned int rows, unsigned int cols, double v)
{
    data = new double[rows * cols];
    this->rows = rows;
    this->cols = cols;

#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < rows * cols; i++)
    {
        data[i] = v;
    }
}
/**/
inline Matrix::Matrix(std::initializer_list<double> const &init)
{
    cols = init.size();
    rows = 1;
    data = new double[rows * cols];

    int i = 0;
    for (auto itr = init.begin(); itr != init.end(); ++itr)
    {
        data[i] = *itr;
        i++;
    }
}

/*デストラクタ*/
inline Matrix::~Matrix()
{
    delete[] data;
}

/****************** 一般関数 ********************/
/*コピー*/
inline void Matrix::copy(Matrix const &m)
{
    assert(m.cols *m.rows != 0);
    if (m.cols * m.rows != cols * rows)
    {
        delete[] data;
        data = new double[rows * cols];
        cols = m.cols;
        rows = m.rows;
    }

    for (size_t i = 0; i < m.rows * m.cols; i++)
    {
        data[i] = m.data[i];
    }
}
/*行列の形*/
inline std::pair<unsigned int, unsigned int> Matrix::shape() const
{
    return std::make_pair(rows, cols);
}

/*行列の行数*/
inline unsigned int Matrix::getRows() const
{
    return rows;
}

/*行列の列数*/
inline unsigned int Matrix::getCols() const
{
    return cols;
}

/*行列の形の変更*/
inline void Matrix::reshape(unsigned int new_rows, unsigned int new_cols)
{
    assert(new_cols * new_rows == cols * rows);

    cols = new_cols;
    rows = new_rows;
}

/****************** 数学関数 ********************/

/*行列積を計算
b:*/
inline Matrix Matrix::dot(const Matrix &b) const
{
    assert(cols == b.rows);

    Matrix ret(rows, b.cols);

#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < b.cols; j++)
        {
            ret(i, j) = 0;

            for (int k = 0; k < b.rows; k++)
            {
                ret(i, j) += operator()(i, k) * b(k, j);
            }
        }
    }
    return ret;
}
/*転置*/
inline Matrix Matrix::T() const
{
    Matrix ret(cols, rows);

#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            ret(j, i) = operator()(i, j);
        }
    }
    return ret;
}
/*L2 normを計算*/
inline double Matrix::l2norm() const
{
    double ret = 0;

#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < cols * rows; i++)
    {
        ret += data[i] * data[i];
    }
    return sqrt(ret);
}
/****************** オペレーターのオーバーロード ********************/

/*要素の参照*/
inline double Matrix::operator()(unsigned int row, unsigned int col) const
{
    assert(row < rows);
    assert(col < cols);
    return data[row * cols + col];
}
/*要素の参照その2*/
inline double Matrix::operator()(unsigned int index) const
{
    assert(index < rows * cols);
    return data[index];
}

/*要素の代入*/
inline double &Matrix::operator()(unsigned int row, unsigned int col)
{
    assert(row < rows);
    assert(col < cols);
    return data[row * cols + col];
}
/*要素の代入その2*/
inline double &Matrix::operator()(unsigned int index)
{
    assert(index < rows * cols);
    return data[index];
}

/*行列の和*/
inline Matrix Matrix::operator+(Matrix const &m2) const
{
    Matrix ret(rows, cols);
    assert(rows == m2.rows && cols == m2.cols);

#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < rows * cols; i++)
    {
        ret.data[i] = data[i] + m2.data[i];
    }
    return ret;
}

/*行列の差*/
inline Matrix Matrix::operator-(Matrix const &m2) const
{
    return operator+(m2.operator*(-1.0));
}

/*行列の等号*/
inline bool Matrix::operator==(Matrix const &m2) const
{
    if (m2.rows != rows || m2.cols != cols)
    {
        return false;
    }

    bool isdiffrent = true;

#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < rows * cols; i++)
    {
        if (data[i] != m2.data[i])
        {
            isdiffrent = false;
        }
    }
    return isdiffrent;
}

/*行列のスカラー倍*/
inline Matrix Matrix::operator*(double a) const
{
    Matrix ret(rows, cols);

#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < rows * cols; i++)
    {
        ret.data[i] = data[i] * a;
    }
    return ret;
}
/*行列のスカラー倍その2*/
Matrix operator*(double a, Matrix const &b)
{
    return b * a;
}

/*行列のスカラーの逆数倍*/
inline Matrix Matrix::operator/(double a) const
{
    return operator*(1.0 / a);
}

/*iniitializerを用いた初期化(1次元)*/
inline void Matrix::operator=(std::initializer_list<double> const &init)
{
    cols = init.size();
    rows = 1;
    data = new double[rows * cols];

    int i = 0;
    for (auto itr = init.begin(); itr != init.end(); ++itr)
    {
        data[i] = *itr;
        i++;
    }
}


/*std::coutへの対応*/
std::ostream &operator<<(std::ostream &os, const Matrix &m)
{
    os << "Matrix([";
    for (size_t i = 0; i < m.getRows(); i++)
    {
        if (i > 0)
        {
            os << "        ";
        }

        os << "[";
        for (int j = 0; j < m.getCols(); j++)
        {
            os << m(i, j);
            if (j != m.getCols() - 1)
            {
                os << ",";
            }
        }
        os << "]";

        if (i != m.getRows() - 1)
        {
            os << "," << std::endl;
        }
    }

    os << "])";
    return os;
}
} // namespace QLA

#endif