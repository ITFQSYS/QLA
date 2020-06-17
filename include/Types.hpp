#ifndef __QLA_TYPES_HPP__
#define __QLA_TYPES_HPP__ 1
#include <initializer_list>
#include <utility>
#include <vector>

namespace QLA
{

class Matrix
{
protected: //データ
    unsigned int rows;
    unsigned int cols;
    double *data;

public: //メソッド
    Matrix();
    Matrix(unsigned int rows, unsigned int cols);
    Matrix(unsigned int rows, unsigned int cols, double v);
    Matrix(std::initializer_list<double> const &init);
    Matrix(const Matrix &rhs);
    
    ~Matrix(); //デストラクタ

    /*一般*/
    void copy(Matrix const &m);
    std::pair<unsigned int, unsigned int> shape() const;
    unsigned int getRows() const;
    unsigned int getCols() const;
    void reshape(unsigned int new_rows, unsigned int new_cols);
    //void resize(unsigned int new_rows, unsigned int new_cols);

    /*数学*/
    Matrix dot(const Matrix &B) const;
    Matrix T() const;
    double l2norm() const;
    //Matrix inverse() const; /*逆行列*/

    /*演算子オーバーロード*/
    Matrix& operator=(const Matrix &m);
    double operator()(unsigned int row, unsigned int col) const;
    double operator()(unsigned int index) const;
    double &operator()(unsigned int row, unsigned int col);
    double &operator()(unsigned int index);
    Matrix operator+(Matrix const &m2) const;
    Matrix operator-(Matrix const &m2) const;
    bool operator==(Matrix const &m2) const;
    Matrix operator*(double a) const;
    Matrix operator/(double a) const;
    void operator=(std::initializer_list<double> const &init);
};

// class Vector :public Matrix
// {
// private:
//     /* data */
// public:
//     Vector(/* args */);
//     Vector(unsigned int dims);
//     Vector(std::vector<double> const &v);
//     Vector(std::initializer_list<double> init);
//     ~Vector();

//     double operator()(unsigned int index) const;
//     //double &operator()(unsigned int index);
// };

} // namespace QLA
#endif