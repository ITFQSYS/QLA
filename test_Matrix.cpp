#include <iostream>
#include <vector>
#include <string>
#include "qla.hpp"

using namespace std;
using namespace QLA;

bool Matrix_eq(Matrix const &m, vector<vector<double>> &ref)
{
    int rows = m.getRows();
    int cols = m.getCols();

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            double m_ij = m(i, j);
            if (m_ij != ref[i][j])
            {
                return false;
            }
        }
    }
    return true;
}

string ok_ng(bool result)
{
    if (result)
    {
        return string("OK");
    }
    else
    {
        return string("NG");
    }
}

bool test_operator()
{
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    vector<vector<double>> b_ref = {{5, 6}, {7, 8}};

    Matrix a(2, 2);
    Matrix b(2, 2);

    /*代入のテスト*/
    bool result_assign = true;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
            b(i, j) = b_ref[i][j];
        }
    }
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            result_assign = a(i, j) == a_ref[i][j];
        }
    }
    cout << "Operator\tMatrix assign:\t";
    if (result_assign)
    {
        std::cout << "OK" << std::endl;
    }
    else
    {
        std::cout << "NG" << std::endl;
    }

    /*等号のテスト*/
    bool result_equal;

    result_equal = (a == a);
    std::cout << "Operator\tMatrix equal:\t";
    if (result_equal)
    {
        cout << "OK" << std::endl;
    }
    else
    {
        cout << "NG" << std::endl;
    }

    /*加算のテスト*/
    vector<vector<double>> apb_ref = {{1 + 5, 2 + 6}, {3 + 7, 4 + 8}};
    bool result_add = Matrix_eq(a + b, apb_ref);

    cout << "Operator\tMatrix add:\t";
    if (result_add)
    {
        std::cout << "OK" << std::endl;
    }
    else
    {
        cout << "NG" << std::endl;
    }

    /*減算のテスト*/
    vector<vector<double>> amb_ref = {{1 - 5, 2 - 6}, {3 - 7, 4 - 8}};
    bool result_sub = Matrix_eq(a - b, amb_ref);
    cout << "Operator\tMatrix sub:\t";
    if (result_sub)
    {
        std::cout << "OK" << std::endl;
    }
    else
    {
        cout << "NG" << std::endl;
    }

    /*スカラー倍のテスト*/
    double alpha = 2;

    vector<vector<double>> alpha_a_ref = {{1 * alpha, 2 * alpha}, {3 * alpha, 4 * alpha}};
    Matrix alpha_a = alpha * a;
    Matrix a_alpha = a * alpha;

    bool result_scalar_mul_ops = a_alpha == alpha_a;
    bool result_scalar_mul_calc = Matrix_eq(a_alpha, alpha_a_ref);
    bool result_scalar_mul = result_scalar_mul_calc & result_scalar_mul_ops;

    cout << "Operator\tScalar mul:\t";
    if (result_scalar_mul)
    {
        std::cout << "OK" << std::endl;
    }
    else
    {
        cout << "NG" << std::endl;
    }

    /*1/スカラー倍のテスト*/
    vector<vector<double>> b_alpha_ref = {{5.0 / alpha, 6.0 / alpha}, {7.0 / alpha, 8.0 / alpha}};
    Matrix b_alpha = b / alpha;
    bool result_scalar_div = Matrix_eq(b_alpha, b_alpha_ref);
    cout << "Operator\tScalar div:\t";
    if (result_scalar_div)
    {
        std::cout << "OK" << std::endl;
    }
    else
    {
        cout << "NG" << std::endl;
    }

    Matrix _3a_p_2b = 3 * a + 2 * b;
    Matrix b2_p_a3 = b / (1.0 / 2.0) + a / (1.0 / 3.0);
    vector<vector<double>> answer_integrated = {{13, 18}, {23, 28}};
    bool result_integrated = Matrix_eq(_3a_p_2b, answer_integrated);
    result_integrated = Matrix_eq(b2_p_a3, answer_integrated);
    cout << "Operator\tintegrated:\t";
    if (result_integrated)
    {
        std::cout << "OK" << std::endl;
    }
    else
    {
        cout << "NG" << std::endl;
    }

    /*初期化のテスト*/
    Matrix c;
    c = {1.0, 2.0, 3.0, 4.0};
    vector<vector<double>> c_ref = {{1, 2, 3, 4}};
    bool result_init = Matrix_eq(c, c_ref);
    cout << "Operator\tinintialize:\t";
    if (result_init)
    {
        std::cout << "OK" << std::endl;
    }
    else
    {
        cout << "NG" << std::endl;
    }

    /*表示のテスト*/
    cout << a<<endl;
    cout <<c<<endl;


    return result_assign && result_equal && result_add && result_sub && result_scalar_mul && result_scalar_div && result_integrated && result_init;
}

bool test_mathematics_dot()
{
    Matrix a(2, 2);
    Matrix b(2, 2);
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    vector<vector<double>> b_ref = {{5, 6}, {7, 8}};
    vector<vector<double>> answer = {{19, 22}, {43, 50}};

    bool result_assign = true;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
            b(i, j) = b_ref[i][j];
        }
    }
    Matrix ab = a.dot(b);
    return Matrix_eq(ab, answer);
}
bool test_mathematics_T()
{
    Matrix a(2, 2);
    Matrix b(2, 2);
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    vector<vector<double>> answer = {{1, 3}, {2, 4}};

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
        }
    }

    return Matrix_eq(a.T(), answer);
}

bool test_mathematics()
{
    bool ret = true;
    bool res;

    res = test_mathematics_dot();
    cout << "Mathematics\tdot\t" << ok_ng(res) << "\n";
    ret = ret && res;

    res = test_mathematics_T();
    cout << "Mathematics\tT\t" << ok_ng(res) << "\n";
    ret = ret && res;

    return res;
}
int main(int argc, char const *argv[])
{
    test_operator();
    test_mathematics();
    return 0;
}
