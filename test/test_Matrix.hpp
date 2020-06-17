#include "../include/qla.hpp"
#include <gtest/gtest.h>

#include <vector>

using namespace QLA;

/*代入のテスト*/
TEST(MatrixTest, operator_assign_Test)
{
    using namespace std;
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    Matrix a(2, 2);
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
        }
    }

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            EXPECT_TRUE(a(i, j) == a_ref[i][j]);
        }
    }
}

/*サイズ変更のテスト*/
TEST(MatrixTest, resize_test)
{
    Matrix a(1, 1);

    for (unsigned int i = 2; i < 1000; i++)
    {
        a.resize(i, i);

        EXPECT_EQ(a.getRows(),i);
        EXPECT_EQ(a.getCols(),i);
        
        for (size_t j = 0; j < i; j++)
        {
            for (size_t k = 0; k < i; k++)
            {
                a(j,k)=0;
            }
        }
    }
}

/*等号のテスト*/
TEST(MatrixTest, operator_equal_Test)
{
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    vector<vector<double>> b_ref = {{1, 1}, {3, 4}};

    Matrix a(2, 2);
    Matrix b(2, 2);

    bool result_assign = true;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
            b(i, j) = b_ref[i][j];
        }
    }
    EXPECT_TRUE(a == a);
    EXPECT_FALSE(a == b);
}

/*加算のテスト*/
TEST(MatrixTest, operator_add_Test)
{
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    vector<vector<double>> b_ref = {{5, 6}, {7, 8}};

    Matrix a(2, 2);
    Matrix b(2, 2);

    bool result_assign = true;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
            b(i, j) = b_ref[i][j];
        }
    }
    Matrix c = a + b;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            EXPECT_DOUBLE_EQ(c(i, j), a_ref[i][j] + b_ref[i][j]);
        }
    }
}

/*減算のテスト*/
TEST(MatrixTest, operator_sub_Test)
{
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    vector<vector<double>> b_ref = {{5, 6}, {7, 8}};

    Matrix a(2, 2);
    Matrix b(2, 2);

    bool result_assign = true;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
            b(i, j) = b_ref[i][j];
        }
    }
    Matrix c = a - b;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            EXPECT_DOUBLE_EQ(c(i, j), a_ref[i][j] - b_ref[i][j]);
        }
    }
}

/*スカラー倍のテスト*/
TEST(MatrixTest, scalar_prod_test)
{
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    vector<vector<double>> b_ref = {{5, 6}, {7, 8}};

    double alpha = 2;
    Matrix a(2, 2);
    Matrix b(2, 2);

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
            b(i, j) = b_ref[i][j];
        }
    }

    vector<vector<double>> alpha_a_ref = {{1 * alpha, 2 * alpha}, {3 * alpha, 4 * alpha}};
    Matrix alpha_a = alpha * a;
    Matrix a_alpha = a * alpha;

    EXPECT_TRUE(alpha_a == a_alpha);

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            EXPECT_DOUBLE_EQ(a_alpha(i, j), alpha_a_ref[i][j]);
        }
    }
}

/*1/スカラー倍のテスト*/
TEST(MatrixTest, scalar_div_test)
{
    Matrix b(2, 2);
    vector<vector<double>> b_ref = {{5, 6}, {7, 8}};
    double alpha = 2;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            b(i, j) = b_ref[i][j];
        }
    }
    vector<vector<double>> b_alpha_ref = {{5.0 / alpha, 6.0 / alpha}, {7.0 / alpha, 8.0 / alpha}};
    Matrix b_alpha = b / alpha;
}

/*オペレーターの結合テスト*/
TEST(MatrixTest, operator_integration_test)
{
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    vector<vector<double>> b_ref = {{5, 6}, {7, 8}};

    Matrix a(2, 2);
    Matrix b(2, 2);

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
            b(i, j) = b_ref[i][j];
        }
    }

    Matrix _3a_p_2b = 3 * a + 2 * b;
    vector<vector<double>> answer_integrated = {{13, 18}, {23, 28}};

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            EXPECT_DOUBLE_EQ(_3a_p_2b(i, j), answer_integrated[i][j]);
        }
    }
}

/*初期化のテスト*/
TEST(MatrixTest, initializer_test)
{
    Matrix m = {1.0, 2.0, 3.0, 4.0};
    vector<vector<double>> ref = {{1, 2, 3, 4}};
    int rows = m.getRows();
    int cols = m.getCols();

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            EXPECT_DOUBLE_EQ(m(i, j), ref[i][j]);
        }
    }
}

/*行列積のテスト*/
TEST(MatrixTest, mat_mul_test)
{
    Matrix a(2, 2);
    Matrix b(2, 2);
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    vector<vector<double>> b_ref = {{5, 6}, {7, 8}};
    vector<vector<double>> ref = {{19, 22}, {43, 50}};

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
            b(i, j) = b_ref[i][j];
        }
    }
    Matrix ab = a.dot(b);

    int rows = ab.getRows();
    int cols = ab.getCols();

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            EXPECT_DOUBLE_EQ(ab(i, j), ref[i][j]);
        }
    }
}

/*ベクトルとの積のテスト*/
TEST(MatrixTest, vec_mul_test)
{
    Matrix a(2, 2);
    Matrix b(2, 1);
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    vector<vector<double>> b_ref = {{5}, {7}};
    vector<vector<double>> ref = {{19}, {43}};

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
        }
        b(i) = b_ref[i][0];
    }

    Matrix ab = a.dot(b);

    int rows = ab.getRows();
    int cols = ab.getCols();

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            EXPECT_DOUBLE_EQ(ab(i, j), ref[i][j]);
        }
    }
}

/*ベクトルとの積のテスト*/
TEST(MatrixTest, transpose_test)
{
    Matrix a(2, 2);
    vector<vector<double>> a_ref = {{1, 2}, {3, 4}};
    vector<vector<double>> answer = {{1, 3}, {2, 4}};

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            a(i, j) = a_ref[i][j];
        }
    }
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            EXPECT_DOUBLE_EQ(a.T()(i, j), answer[i][j]);
        }
    }
}