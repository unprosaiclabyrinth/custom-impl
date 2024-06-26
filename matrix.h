/*
 * file: matrix.h
 * Implementing a matrix ADT to carry out linear algebraic operations
 * System: Mac using CLion
 * Author: Himanshu Dongre
 */

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <sstream>

using namespace std;

typedef vector<vector<Complex>> matrix;
typedef vector<Complex> Row;
typedef vector<Complex> Col;
typedef int rowno;
typedef int colno;
typedef pair<rowno, colno> Pos;

// Matrix is 00-indexed
class Matrix {
private:
    int nrows;
    int ncols;
    matrix M;

    //
    // max_lengths
    //
    // private member helper function
    // used in pretty printing the matrix
    // takes two boolean parameters:
    // b1 and b2 => consider the elements in cis form
    // b1 but not b2 => consider the elements in eulers form
    // not b1 and not b2 => consider the elements in cartesian form
    //
    vector<size_t> max_lengths(const bool &to_print = false, const bool &cis = false) const;

    struct iterator {
    private:
        Matrix *enclosing;
        rowno r;
        colno c;
    public:
        iterator(Matrix *enclosing, int r, int c) {
            this->enclosing = enclosing;
            this->r = r;
            this->c = c;
        }

        Complex operator*() {
            return enclosing->at(r, c);
        }

        bool operator!=(const iterator &rhs) {
            return this->r != rhs.r or this->c != rhs.c;
        }

        iterator operator++() {
            if (c < enclosing->ncols - 1) {
                ++c;
            } else if (r < enclosing->nrows - 1) {
                ++r;
                c = 0;
            } else {
                ++r;
                ++c;
            }
            return *this;
        }
    };

public:
    /// Constructors

    //
    // default constructor
    //
    Matrix();

    //
    // taking parameters for number of rows and number of columns
    // and initializing the matrix to a null matrix of that size
    //
    explicit Matrix(int nR, int nC);

    //
    // initializing a square matrix of given length to an identity matrix
    //
    explicit Matrix(int len);

    //
    // initializing the matrix from a vector of vectors
    // doesn't allow for (non-rectangular) jagged matrices by throwing a run-time error
    //
    template<typename T> Matrix(const vector<vector<T>> &M);
    // Matrix(const matrix &mat);

    /// Iterator methods

    //
    // begin()
    //
    iterator begin();

    //
    // end()
    //
    iterator end();

    /// Linear Algebraic & Logical operations

    //
    // operator+
    //
    // adding two matrices
    // throwing a runtime error if matrices are of unequal sizes
    //
    friend Matrix operator+(const Matrix &M1, const Matrix &M2);

    //
    // operator-
    //
    // subtracting one matrix from another
    // throwing a runtime error if matrices are of unequal sizes
    //
    friend Matrix operator-(const Matrix &M1, const Matrix &M2);

    //
    // operator*
    //
    // multiplying two matrices
    // throwing a runtime error if matrices don't allow multiplication
    //
    friend Matrix operator*(const Matrix &M1, const Matrix &M2);
    //
    // overloading operator* to be able to
    // multiply the matrix by a constant
    //
    friend Matrix operator*(const Complex &k, const Matrix &M1);

    //
    // operator/
    //
    // dividing the matrix by a constant
    //
    friend Matrix operator/(const Matrix &M1, const Complex &k);

    //
    // operator==
    //
    // checking whether two matrices are equal
    //
    friend bool operator==(const Matrix &M1, const Matrix &M2);
    //
    // operator!=
    //
    friend bool operator!=(const Matrix &M1, const Matrix &M2);

    /// Matrix operations

    //
    // operator()
    //
    // returning a reference to an element
    // with bounds checking
    // no difference between this and at()
    //
    Complex &operator()(const rowno &r, const colno &c);
    //
    // overloading operator() to work with a pair input
    //
    Complex &operator()(const Pos &pos);

    //
    // at
    //
    // method returning a reference to an element
    // with bounds checking
    //
    Complex &at(const rowno &r, const colno &c);
    //
    // overloading at to work with a pair input
    //
    Complex &at(const Pos &pos);

    //
    // set_row
    //
    // sets the row at row index r to the vector
    // passed in as a parameter
    //
    void set_row(const rowno &r, const Row &row);

    //
    // set_col
    //
    // sets the column at column index c to the vector
    // passed in as a parameter
    //
    void set_col(const colno &c, const Col &col);

    //
    // get_row
    //
    // method returning a reference to a row in the matrix
    //
    Row get_row(const rowno &r) const;

    //
    // get_col
    //
    // method returning a reference to a column in the matrix
    //
    Col get_col(const colno &c) const;

    //
    // posOf
    //
    // returns a pair of the row number and the column number of an element
    //
    Pos posOf(const Complex &elem) const;

    //
    // contains
    //
    // returns true if the element is present in the matrix
    //
    bool contains(const Complex &elem) const;

    //
    // order
    //
    // displaying the order of the matrix
    //
    void display_order(ostream &out) const;

    //
    // row_count
    //
    // returning the number of rows
    //
    int row_count() const;

    //
    // col_count
    //
    // returning the number of rows
    //
    int col_count() const;

    //
    // size
    //
    // returning the number of rows
    //
    int size() const;

    /// Advanced Matrix Operations

    //
    // transpose
    //
    // returning the transpose of the matrix
    //
    Matrix transpose();

    //
    // swap_rows
    //
    // swaps two rows in the matrix
    //
    void swap_rows(const rowno &r1, const rowno &r2);

    //
    // swap_cols
    //
    // swaps two columns in the matrix
    //
    void swap_cols(const colno &c1, const colno &c2);

    /* Square Matrix operations */

    //
    // trace
    //
    // returning the trace of the matrix
    //
    Complex trace();

    //
    // cofactor_matrix
    //
    // returning the cofactor matrix of a given element
    //
    Matrix cofactor_matrix(rowno r, colno c);

    //
    // cofactor
    //
    // returning the cofactor of a given element
    //
    Complex cofactor(rowno r, colno c);

    //
    // det
    //
    // returning the determinant of the matrix
    // recursion
    //
    Complex det();

    //
    // adj
    //
    // returning the adjoint of the matrix
    //
    Matrix adjoint();

    //
    // inverse
    //
    // returning the inverse of the matrix
    //
    Matrix inverse();

    //
    // conjugate
    //
    // returning the conjugate of the matrix
    //
    Matrix conjugate();

    //
    // hermitian_conjugate
    //
    // returning the hermitian conjugate of the matrix
    //
    Matrix hermitian_conjugate();

    /// Matrix I/O (Actually only O)

    //
    // print
    //
    // printing the matrix to an output stream
    // with all the elements in their polar form
    // prints the elements in their cis forms if the boolean
    // parameter, cis, is true
    //
    void print(ostream &out, const bool &cis = false);

    //
    // operator<<
    //
    // printing the matrix to an output stream
    // with all the elements in their cartesian form
    //
    friend ostream &operator<<(ostream &out, const Matrix &M1);
};


































































































//
// vvv IMPLEMENTATION ABSTRACTED AWAY vvv
//

































































































//
// FUNCTION DEFINITIONS
//

vector<size_t> Matrix::max_lengths(const bool &to_print, const bool &cis) const
{
    vector<size_t> max_lengths;
    for (colno c = 0; c < ncols; ++c) {
        size_t ml = 0;
        for (Complex e : get_col(c)) {
            stringstream elem;
            if (not to_print) elem << setprecision(4) << e;
            else if (not cis) e.print(elem);
            else e.print(elem, true);
            if (elem.str().size() > ml) {
                ml = elem.str().size();
            }
        }
        max_lengths.push_back(ml);
    }
    return max_lengths;
}

Matrix::Matrix() {
    nrows = 0;
    ncols = 0;
    M = {};
}

Matrix::Matrix(int nR, int nC)
{
    if (nR < 0 or nC < 0) {
        throw runtime_error("/*** DIMENSIONS OF MATRIX MUST BE NON-NEGATIVE ***/");
    } else if (nR == 0 or nC == 0) {
        nrows = 0;
        ncols = 0;
        M = {};
        return;
    }
    nrows = nR;
    ncols = nC;
    M = {};
    for (rowno r = 0; r < nrows; ++r) {
        Col column;
        for (colno c = 0; c < ncols; ++c) {
            column.emplace_back(0);
        }
        M.push_back(column);
    }
}

Matrix::Matrix(int len)
{
    if (len < 2) {
        throw runtime_error("/*** LENGTH OF IDENTITY MATRIX MUST BE AT LEAST 2 ***/");
    }
    nrows = len;
    ncols = len;
    for (rowno r = 0; r < nrows; ++r) {
        Col column;
        for (colno c = 0; c < ncols; ++c) {
            if (r == c) {
                column.emplace_back(1);
            } else {
                column.emplace_back(0);
            }
        }
        M.push_back(column);
    }
}

template<typename T>
Matrix::Matrix(const vector<vector<T>> &M)
// Matrix::Matrix(const matrix &mat)
{
    nrows = int(M.size());
    if (!M.empty()) {
        ncols = int(M[0].size());
        matrix Mv;
        for (rowno r = 0; r < nrows; ++r) {
            if (M[r].size() != ncols) {
                throw runtime_error("/*** MATRIX MUST BE RECTANGULAR ***/");
            }
            Row currv;
            for (colno c = 0; c < ncols; ++c) {
                currv.push_back(M[r][c]);
            }
            Mv.push_back(currv);
        }
        this->M = Mv;
    } else {
        ncols = 0;
    }
}

Matrix::iterator Matrix::begin() {
    return {this, 0, 0};
}

Matrix::iterator Matrix::end() {
    return {this, nrows, ncols};
}

Matrix operator+(const Matrix &M1, const Matrix &M2)
{
    if (M1.nrows != M2.nrows or M1.ncols != M2.ncols) {
        throw runtime_error("/*** sum: OPERANDS MUST HAVE EQUAL SIZES ***/");
    } else {
        matrix sumv;
        for (rowno r = 0; r < M1.nrows; ++r) {
            Row currsum;
            for (colno c = 0; c < M1.ncols; ++c) {
                currsum.push_back(M1.M[r][c] + M2.M[r][c]);
            }
            sumv.push_back(currsum);
        }
        Matrix sum(sumv);
        return sum;
    }
}

Matrix operator-(const Matrix &M1, const Matrix &M2) {
    if (M1.nrows != M2.nrows or M1.ncols != M2.ncols) {
        throw runtime_error("/*** difference: OPERANDS MUST HAVE EQUAL SIZES ***/");
    } else {
        matrix diffv;
        for (rowno r = 0; r < M1.nrows; ++r) {
            Row currdiff;
            for (colno c = 0; c < M1.ncols; ++c) {
                currdiff.push_back(M1.M[r][c] - M2.M[r][c]);
            }
            diffv.push_back(currdiff);
        }
        Matrix diff(diffv);
        return diff;
    }
}

Matrix operator*(const Matrix &M1, const Matrix &M2) {
    if (M1.ncols != M2.nrows) {
        throw runtime_error("/*** product: COLUMNS OF FIRST MUST BE EQUAL TO ROWS OF SECOND ***/");
    } else {
        matrix prodv;
        for (rowno r1 = 0; r1 < M1.nrows; ++r1) {
            Row currprod;
            for (colno c2 = 0; c2 < M2.ncols; ++c2) {
                Complex sum;
                for (rowno i = 0; i < M2.nrows; ++i) {
                    sum += M1.M[r1][i] * M2.M[i][c2];
                }
                currprod.push_back(sum);
            }
            prodv.push_back(currprod);
        }
        Matrix prod(prodv);
        return prod;
    }
}

Matrix operator*(const Complex &k, const Matrix &M1) {
    matrix resv;
    for (rowno r = 0; r < M1.nrows; ++r) {
        Row currres;
        for (colno c = 0; c < M1.ncols; ++c) {
            currres.push_back(k * M1.M[r][c]);
        }
        resv.push_back(currres);
    }
    Matrix res(resv);
    return res;
}

Matrix operator/(const Matrix &M1, const Complex &k) {
    if (k == "0") {
        throw runtime_error("/*** DIVISION BY 0 ***/");
    } else {
        matrix resv;
        for (rowno r = 0; r < M1.nrows; ++r) {
            Row currres;
            for (colno c = 0; c < M1.ncols; ++c) {
                currres.push_back(M1.M[r][c] / k);
            }
            resv.push_back(currres);
        }
        Matrix res(resv);
        return res;
    }
}

bool operator==(const Matrix &M1, const Matrix &M2)
{
    return M1.M == M2.M;
}

bool operator!=(const Matrix &M1, const Matrix &M2)
{
    return M1.M != M2.M;
}

Complex &Matrix::operator()(const rowno &r, const colno &c)
{
    if (r >= 0 and r < nrows and c >= 0 and c < ncols) {
        return M[r][c];
    } else {
        throw out_of_range("/*** at: OUT OF BOUNDS ***/");
    }
}

Complex &Matrix::operator()(const Pos &pos)
{
    int r = pos.first, c = pos.second;
    if (r >= 0 and r < nrows and c >= 0 and c < ncols) {
        return M[r][c];
    } else {
        throw out_of_range("/*** at: OUT OF BOUNDS ***/");
    }
}

Complex &Matrix::at(const rowno &r, const colno &c)
{
    if (r >= 0 and r < nrows and c >= 0 and c < ncols) {
        return M[r][c];
    } else {
        throw out_of_range("/*** at: OUT OF BOUNDS ***/");
    }
}

Complex &Matrix::at(const Pos &pos)
{
    int r = pos.first, c = pos.second;
    if (r >= 0 and r < nrows and c >= 0 and c < ncols) {
        return M[r][c];
    } else {
        throw out_of_range("/*** at: OUT OF BOUNDS ***/");
    }
}

void Matrix::set_row(const rowno &r, const Row &row)
{
    if (r < 0 or r > nrows - 1) {
        throw out_of_range("/*** set_row: r: OUT OF BOUNDS ***/");
    }
    if (row.size() != ncols) {
        throw out_of_range("/*** set_row: row: INCORRECT SIZE ***/");
    }
    M[r] = row;
}

void Matrix::set_col(const colno &c, const Col &col)
{
    if (c < 0 or c > ncols - 1) {
        throw out_of_range("/*** set_col: c: OUT OF BOUNDS ***/");
    }
    if (col.size() != nrows) {
        throw out_of_range("/*** set_col: col: INCORRECT SIZE ***/");
    }
    for (rowno r = 0; r < nrows; ++r) {
        M[r][c] = col[r];
    }
}

Row Matrix::get_row(const rowno &r) const
{
    if (r >= 0 and r < nrows) {
        return M[r];
    } else {
        throw out_of_range("/*** get_row: OUT OF BOUNDS ***/");
    }
}

Col Matrix::get_col(const colno &c) const
{
    if (c >= 0 and c < ncols) {
        Col column;
        for (rowno r = 0; r < nrows; ++r) {
            column.push_back(M[r][c]);
        }
        return column;
    } else {
        throw out_of_range("/*** get_col: OUT OF BOUNDS ***/");
    }
}

Pos Matrix::posOf(const Complex &elem) const
{
    for (rowno r = 0; r < nrows; ++r) {
        for (colno c = 0; c < ncols; ++c) {
            if (M[r][c] == elem) {
                return {r, c};
            }
        }
    }
    return {-1, -1};
}

bool Matrix::contains(const Complex &elem) const
{
    return posOf(elem).first > -1;
}

void Matrix::display_order(ostream &out) const
{
    out << nrows << " x " << ncols << endl;
}

int Matrix::row_count() const
{
    return nrows;
}

int Matrix::col_count() const
{
    return ncols;
}

int Matrix::size() const
{
    return nrows * ncols;
}


Matrix Matrix::transpose()
{
    matrix Mtv;
    for (rowno r = 0; r < ncols; ++r) {
        Row rowt;
        for (colno c = 0; c < nrows; ++c) {
            rowt.push_back(M[c][r]);
        }
        Mtv.push_back(rowt);
    }
    Matrix Mt(Mtv);
    return Mt;
}

void Matrix::swap_rows(const rowno &r1, const rowno &r2)
{
    if (r1 < 0 or r1 > nrows - 1) {
        throw out_of_range("/*** swap_rows: r1: OUT OF BOUNDS ***/");
    }
    if (r2 < 0 or r2 > nrows - 1) {
        throw out_of_range("/*** swap_rows: r2: OUT OF BOUNDS ***/");
    }
    Row temp = M[r1];
    M[r1] = M[r2];
    M[r2] = temp;
}

void Matrix::swap_cols(const colno &c1, const colno &c2)
{
    if (c1 < 0 or c1 > ncols - 1) {
        throw out_of_range("/*** swap_cols: c1: OUT OF BOUNDS ***/");
    }
    if (c2 < 0 or c2 > ncols - 1) {
        throw out_of_range("/*** swap_cols: c2: OUT OF BOUNDS ***/");
    }
    Col temp = get_col(c1);
    for (rowno r = 0; r < nrows; ++r) {
        M[r][c1] = M[r][c2];
        M[r][c2] = temp[r];
    }
}

Complex Matrix::trace()
{
    if (nrows != ncols) {
        throw runtime_error("*** TRACE UNDEFINED FOR A NON-SQUARE MATRIX ***/");
    }
    Complex tr;
    for (rowno r = 0; r < nrows; ++r) {
        for (colno c = 0; c < ncols; ++c) {
            if (r == c) {
                tr += M[r][c];
            }
        }
    }
    return tr;
}

Matrix Matrix::cofactor_matrix(rowno r, colno c)
{
    if (nrows != ncols) {
        throw runtime_error("/*** COFACTOR MATRIX UNDEFINED FOR A NON-SQUARE MATRIX ***/");
    }
    else if (r >= 0 and r < nrows and c >= 0 and c < ncols) {
        matrix cofv;
        for (rowno i = 0; i < nrows; ++i) {
            Row currcof;
            for (colno j = 0; j < ncols; ++j) {
                if (i != r and j != c) {
                    currcof.push_back(M[i][j]);
                }
            }
            if (!currcof.empty()) {
                cofv.push_back(currcof);
            }
        }
        Matrix cof(cofv);
        return cof;
    } else {
        throw out_of_range("/*** cofactor_matrix: OUT OF BOUNDS ***/");
    }
}

Complex Matrix::cofactor(rowno r, colno c)
{
    if (nrows != ncols) {
        throw runtime_error("/*** COFACTOR UNDEFINED FOR A NON-SQUARE MATRIX ***/");
    } else if (r >= 0 and r < nrows and c >= 0 and c < ncols) {
        return pow(-1, r + c) * cofactor_matrix(r, c).det();
    } else {
        throw out_of_range("/*** cofactor: OUT OF BOUNDS ***/");
    }
}

Complex Matrix::det()
{
    if (nrows != ncols) {
        throw runtime_error("/*** DETERMINANT UNDEFINED FOR A NON-SQUARE MATRIX ***/");
    }
    else if (size() == 1) {
        return M[0][0];
    } else {
        Complex D;
        for (colno c = 0; c < ncols; ++c) {
            Matrix cof = cofactor_matrix(0, c);
            D += M[0][c] * pow(-1, c) * cof.det();
        }
        return D;
    }
}

Matrix Matrix::adjoint()
{
    if (nrows != ncols) {
        throw runtime_error("/*** ADJOINT UNDEFINED FOR A NON_SQUARE MATRIX ***/");
    }
    matrix adjv;
    for (rowno r = 0; r < nrows; ++r) {
        Row curradj;
        for (colno c = 0; c < ncols; ++c) {
            curradj.push_back(cofactor(r, c));
        }
        adjv.push_back(curradj);
    }
    Matrix adj(adjv);
    adj = adj.transpose();
    return adj;
}

Matrix Matrix::inverse()
{
    if (nrows != ncols) {
        throw runtime_error("/*** INVERSE UNDEFINED FOR A NON-SQUARE MATRIX ***/");
    }
    Complex D = det();
    if (D == "0") {
        throw runtime_error("/*** INVERSE OF A SINGULAR MATRIX DOES NOT EXIST ***/");
    }
    return adjoint() / D;
}

Matrix Matrix::conjugate()
{
    matrix bar;
    for (rowno r = 0; r < nrows; ++r) {
        Row currbar;
        for(colno c = 0; c < ncols; ++c) {
            currbar.push_back(at(r, c).conjugate());
        }
        bar.push_back(currbar);
    }
    return bar;
}

Matrix Matrix::hermitian_conjugate() {
    return transpose().conjugate();
}

void Matrix::print(ostream &out, const bool &cis) {
    vector<size_t> ls = max_lengths(true, cis);
    for (rowno r = 0; r < nrows; ++r) {
        out << "[";
        for (colno c = 0; c < ncols; ++c) {
            stringstream elem;
            M[r][c].print(elem, cis);
            int space = c == 0 ? 1 : 2;
            out << setw(int(ls[c]) + space) << setprecision(4) << elem.str();
        }
        out << " ]" << endl;
    }
    out << endl;
}

ostream &operator<<(ostream &out, const Matrix &M1)
{
    if (M1.size() == 0) {
        out << "[]" << endl;
    } else {
        vector<size_t> ls = M1.max_lengths();
        for (rowno r = 0; r < M1.nrows; ++r) {
            out << "[";
            for (colno c = 0; c < M1.ncols; ++c) {
                stringstream elem;
                elem << setprecision(4) << M1.M[r][c];
                int space = c == 0 ? 1 : 2;
                out << setw(int(ls[c]) + space) << setprecision(4) << elem.str();
            }
            out << " ]" << endl;
        }
    }
    return out;
}

