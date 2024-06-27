/*
 * file: Complex.h
 * Implementing complex numbers
 * Author: Himanshu Dongre
 */

#pragma once

#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>

using namespace std;

#define pi M_PI
#define piby2 M_PI_2

// if the difference between a double and its nearest integer
// is greater than this, the double will be rounded to the int
const double SHORTEN_LIMIT = 0.0075;

// Implementation of complex numbers
class Complex {
private:
    double x;
    double y;
    double r;
    double A;

    //
    // shorten
    //
    // private member helper function
    // shorten a double very close to an integer to that integer
    //
    friend double shorten(const double &d);

    //
    // quadrant
    //
    // private member helper function
    // return an integer unique to the position of the
    // complex number in the plane
    //
    [[nodiscard]] int quadrant() const;

    //
    // magnitude
    //
    // private member helper function
    // calculate the magnitude of the complex number
    // from x and y
    //
    [[nodiscard]] double magnitude() const;

    //
    // theta
    //
    // private member helper function
    // calculate the argument of the complex number
    // from x and y
    //
    [[nodiscard]] double theta() const;

    //
    // extract_from_string
    //
    // private member helper function
    // extract a complex number from a string of one of the forms from:
    // x + yi, x+yi, x - yi, x-yi, re^Ai, rcisA
    //
    void extract_from_string(string num);

public:
    //
    // default constructor
    // convert a real number to a complex number by
    // initializing its complex part to 0
    //
    Complex(double num = 0);

    //
    // initialize a complex number from a string
    // initialize from the forms: x + yi, x+yi, x - yi, x-yi, re^Ai, rcisA
    //
    Complex(const char* num_s);

    //
    // initialize a complex number from either its rectangular coordinates
    // or from its polar coordinates depending upon a boolean parameter
    // if true then initialize to rect coordinates else to polar coordinates
    //
    Complex(double coord1, double coord2, const bool &rect_coords = true);

    //
    // initialize a complex number as cos(d) + isin(d)
    //
    friend Complex cis(double d);

    /* ACCESSORS: */

    //
    // re
    //
    // get the real part of the complex number
    //
    friend double re(const Complex &z);

    //
    // im
    //
    // get the imaginary part of the complex number
    //
    friend double im(const Complex &z);

    //
    // mod
    //
    // get the modulus/magnitude of the complex number
    //
    friend double mod(const Complex &z);

    //
    // arg
    //
    // get the argument of the complex number
    //
    friend double arg(const Complex &z);

    /* MUTATORS: */

    //
    // set_re
    //
    // set the real part of the complex number
    //
    void set_re(double x2);

    //
    // set_im
    //
    // set the imaginary part of the complex number
    //
    void set_im(double y2);

    //
    // set_mod
    //
    // set the modulus/magnitude of the complex number
    //
    void set_mod(double r2);

    //
    // set_arg
    //
    // set the argument of the complex number
    //
    void set_arg(double A2);

    //
    // set
    //
    // set both the rectangular coordinates or the polar coordinates
    // depending upon a boolean parameter
    // if true then set to rectangular coordinates else to polar coordinates
    //
    void set(double coord1, double coord2, const bool &rect_coords = true);

    /* OPERATIONS: */

    //
    // operator==
    //
    // check whether two Complex numbers are equal
    //
    friend bool operator==(const Complex &z1, const Complex &z2);
    //
    // operator!=
    //
    friend bool operator!=(const Complex &z1, const Complex &z2);

    //
    // operator+
    //
    // add two Complex numbers or a Complex number and a real number
    //
    friend Complex operator+(const Complex &z1, const Complex &z2);
    //
    // operator+=
    //
    friend void operator+=(Complex &z1, const Complex &z2);

    //
    // unary operator-
    //
    // negate a complex number
    //
    friend Complex operator-(const Complex &z);

    //
    // binary operator-
    //
    // subtract one complex number from another
    //
    friend Complex operator-(const Complex &z1, const Complex &z2);
    //
    // operator-=
    //
    friend void operator-=(Complex &z1, const Complex &z2);

    //
    // operator*
    //
    // multiply a Complex number to this
    //
    friend Complex operator*(const Complex &z1, const Complex &z2);
    //
    // operator*=
    //
    friend void operator*=(Complex &z1, const Complex &z2);

    //
    // operator/
    //
    // divide this by another Complex number
    //
    friend Complex operator/(const Complex &z1, const Complex &z2);
    //
    // operator/=
    //
    friend void operator/=(Complex &z1, const Complex &z2);

    //
    // dist
    //
    // return the distance between two complex numbers
    //
    friend double dist(const Complex &z1, const Complex &z2);

    //
    // sect
    //
    // return the complex number returned by the section formula
    //
    friend Complex sect(const Complex &z1, const Complex &z2, const double &k);

    //
    // pow
    //
    // raises the complex number to a given complex power
    //
    friend Complex pow(const Complex &z1, const Complex &z2);

    //
    // principal_arg
    //
    // return the principal argument i.e. in the range (-pi, pi]
    //
    friend double principal_arg(const Complex &z);

    //
    // conjugate
    //
    // return the conjugate of the Complex number
    //
    friend Complex conjugate(const Complex &z);

    //
    // is_real
    //
    // return whether the Complex number is real
    //
    friend bool is_real(const Complex &z);

    //
    // cross_ratio
    //
    // return the cross ratio of four Complex numbers
    //
    friend Complex cross_ratio(const Complex &z, const Complex &z1, const Complex &z2, const Complex &z3);

    //
    // rot
    //
    // Given three Complex numbers, z1, z2, z3 out of which z1 and z2 are known;
    // find z3 by rotating z1 about z2 since the magnitudes of the vectors
    // (z1 - z2) and (z3 - z2) are known and so is the angle between them
    //
    friend Complex rotate(const Complex &z, const Complex &about, const double &through,
                          const double &fin_len, const bool &anticlockwise);

    //
    // exp
    //
    // Complex exponential function defined as
    // exp(z = x + iy) = (e^x)(cosy + isiny) = (e^x)(e^iy)
    //
    friend Complex exp(const Complex &z);

    //
    // sin
    //
    // Complex sine function defined as
    // sin(z) = (1/2i)(exp(iz) - exp(-iz))
    //
    friend Complex sin(const Complex &z);

    //
    // cos
    //
    // Complex cosine function defined as
    // cos(z) = (1/2)(exp(iz) + exp(-iz))
    //
    friend Complex cos(const Complex &z);

    //
    // tan
    //
    // Complex tangent function defined as
    // tan(z) = sinz/cosz
    //
    friend Complex tan(const Complex &z);

    //
    // cosec
    //
    // Complex cosecant function defined as
    // cosec(z) = 1/sinz
    //
    friend Complex cosec(const Complex &z);

    //
    // sec
    //
    // Complex secant function defined as
    // sec(z) = 1/cosz
    //
    friend Complex sec(const Complex &z);

    //
    // cot
    //
    // Complex cotangent function defined as
    // cot(z) = cosz/sinz
    //
    friend Complex cot(const Complex &z);

    //
    // sinh
    //
    // Complex hyperbolic sine function defined as
    // sinh(z) = (1/2)(exp(z) - exp(-z))
    //
    friend Complex sinh(const Complex &z);

    //
    // cosh
    //
    // Complex hyperbolic cosine function defined as
    // cosh(z) = (1/2)(exp(z) + exp(-z))
    //
    friend Complex cosh(const Complex &z);

    //
    // tanh
    //
    // Complex hyperbolic tangent function defined as
    // tanh(z) = sinhz/coshz
    //
    friend Complex tanh(const Complex &z);

    //
    // cosech
    //
    // Complex hyperbolic cosecant function defined as
    // cosech(z) = 1/sinhz
    //
    friend Complex cosech(const Complex &z);

    //
    // sech
    //
    // Complex hyperbolic secant function defined as
    // sech(z) = 1/coshz
    //
    friend Complex sech(const Complex &z);

    //
    // coth
    //
    // Complex hyperbolic cotangent function defined as
    // coth(z) = coshz/sinhz
    //
    friend Complex coth(const Complex &z);

    //
    // print
    //
    // print the Complex number to an output stream
    // in either polar(cis) or Euler's form depending on
    // a boolean parameter (true and false respectively)
    //
    // prints in the forms: (r)e^(A), (r)cis(A)
    //
    void print(ostream &out, const bool &cis = false) const;

    //
    // operator<<
    //
    // print the Complex number to an output stream
    // in cartesian form
    //
    // prints in the forms: x + yi, x - yi
    //
    friend ostream &operator<<(ostream &out, const Complex &z1);

    //
    // operator>>
    //
    // extract a Complex number from an input stream
    // in cartesian, polar(cis), or Euler's form
    //
    // extracts correctly from the forms: x + yi, x+yi, x - yi, x-yi, re^Ai, rcisA
    //
    friend istream &operator>>(istream &in, Complex &z1);
};






























































































//
// vvv IMPLEMENTATIONS ABSTRACTED AWAY vvv
//

































































































//
// FUNCTION DEFINITIONS
//

double shorten(const double &d) {
    return fabs(d - round(d)) < SHORTEN_LIMIT ? round(d) : d;
}

int Complex::quadrant() const {
    // -Y-axis
    if (x == 0 and y < 0) return -4;
    // -X-axis
    else if (x < 0 and y == 0) return -3;
    // +Y-axis
    else if (x == 0 and y > 0) return -2;
    // +X-axis
    else if (x > 0 and y == 0) return -1;
    // origin
    else if (x == 0 and y == 0) return 0;
    // quadrant I
    else if (x > 0 and y > 0) return 1;
    // quadrant II
    else if (x < 0 and y > 0) return 2;
    // quadrant III
    else if (x < 0 and y < 0) return 3;
    // quadrant IV
    else return 4;
}

double Complex::magnitude() const {
    return shorten(sqrt((x * x) + (y * y)));
}

double Complex::theta() const {
    switch(quadrant())
    {
        case -4: return -piby2;
        case -3: return pi;
        case -2: return piby2;
        case -1:
        case 0: return 0;
        case 1:
        case 4: return atan(y / x);
        case 2: return pi + atan(y / x);
        case 3: return -pi + atan(y / x);
        default: return -1;
    }
}

void Complex::extract_from_string(string num) {
    if (num.find('i') == string::npos) {
        set(stod(num), 0);
    } else if (num.find("cis") != string::npos) {
        double mod, arg;
        int indexOfC = int(num.find('c'));
        int indexOfS = int(num.find('s'));
        mod = indexOfC == 0 ? 1 : (num[0] == '-' ? -1 : stod(num.substr(0, indexOfC)));
        if (mod < 0) throw runtime_error("/*** cis: MODULUS MUST BE POSITIVE ***/");
        arg = num.at(indexOfS + 1) == '(' ?
              stod(num.substr(indexOfS + 2, num.find(')') - indexOfS - 2)) :
              stod(num.substr(indexOfS + 1));
        set(mod, arg, false);
    } else {
        int indexOfI = int(num.find('i'));
        if (num.find('+') != string::npos) {
            double re, im;
            int indexOfPlus = int(num.find('+'));
            string res = num[indexOfPlus - 1] != ' ' ?
                         num.substr(0, indexOfPlus) :
                         num.substr(0, indexOfPlus - 1);
            re = stod(res);
            string ims = num[indexOfPlus + 1] != ' ' ?
                         num.substr(indexOfPlus + 1, indexOfI - indexOfPlus - 1) :
                         num.substr(indexOfPlus + 2, indexOfI - indexOfPlus - 2);
            im = ims.empty() ? 1 : stod(ims);
            set(re, im);
        } else if (num.find('-') != string::npos) {
            double re, im;
            string res, ims;
            int indexOfMinus = num[0] == '-' ?
                               (int(num.find('-', 1)) == -1 ? 0 :
                                int(num.find('-', 1))) :
                               int(num.find('-'));
            if (indexOfMinus == 0) {
                re = 0;
                ims = num.substr(1, indexOfI - 1);
            } else {
                res = num[indexOfMinus - 1] != ' ' ?
                             num.substr(0, indexOfMinus) :
                             num.substr(0, indexOfMinus - 1);
                re = stod(res);
                ims = num[indexOfMinus + 1] != ' ' ?
                             num.substr(indexOfMinus + 1, indexOfI - indexOfMinus - 1) :
                             num.substr(indexOfMinus + 2, indexOfI - indexOfMinus - 2);
            }
            im = ims.empty() ? -1 : -stod(ims);
            set(re, im);
        } else if (num.find("e^") != string::npos) {
            double mod, arg;
            int indexOfE = int(num.find('e'));
            int indexOfXor = int(num.find('^'));
            mod = indexOfE == 0 ? 1 : (num[0] == '-' ? -1 : stod(num.substr(0, indexOfE)));
            if (mod < 0) throw runtime_error("/*** e^: MODULUS MUST BE POSITIVE ***/");
            arg = num[indexOfXor + 1] == '(' and num[indexOfI - 1] == ')' ?
                  stod(num.substr(indexOfXor + 2, indexOfI - indexOfXor - 3)) :
                  stod(num.substr(indexOfXor + 1, indexOfI - indexOfXor - 1));
            set(mod, arg, false);
        } else {
            double im;
            if (indexOfI == 0) im = 1;
            else {
                string ims = num.substr(0, indexOfI);
                im = stod(ims);
            }
            set(0, im);
        }
    }
}

Complex::Complex(double num) {
    x = shorten(num);
    y = 0;
    r = magnitude();
    A = theta();
}

Complex::Complex(const char *num_s) {
    string num(num_s);
    extract_from_string(num);
}

Complex::Complex(double coord1, double coord2, const bool &rect_coords) {
    if (rect_coords) {
        x = shorten(coord1);
        y = shorten(coord2);
        r = magnitude();
        A = theta();
    } else {
        r = shorten(abs(coord1));
        A = shorten(coord2);
        x = shorten(r * cos(A));
        y = shorten(r * sin(A));
    }
}

Complex cis(double d) {
    Complex z(cos(d), sin(d));
    return z;
}

double re(const Complex &z) {
    return z.x;
}

double im(const Complex &z) {
    return z.y;
}

double mod(const Complex &z) {
    return z.r;
}

double arg(const Complex &z) {
    return z.A;
}

void Complex::set_re(double x2) {
    x = shorten(x2);
    r = magnitude();
    A = theta();
}

void Complex::set_im(double y2) {
    y = shorten(y2);
    r = magnitude();
    A = theta();
}

void Complex::set_mod(double r2) {
    r = shorten(abs(r2));
    x = shorten(r * cos(A));
    y = shorten(r * sin(A));
}

void Complex::set_arg(double A2) {
    A = shorten(A2);
    x = shorten(r * cos(A));
    y = shorten(r * sin(A));
}

void Complex::set(double coord1, double coord2, const bool &rect_coords) {
    if (rect_coords) {
        x = shorten(coord1);
        y = shorten(coord2);
        r = magnitude();
        A = theta();
    } else {
        r = shorten(abs(coord1));
        A = shorten(coord2);
        x = shorten(r * cos(A));
        y = shorten(r * sin(A));
    }
}

bool operator==(const Complex &z1, const Complex &z2) {
    return z1.x == z2.x and z1.y == z2.y;
}

bool operator!=(const Complex &z1, const Complex &z2) {
    return z1.x != z2.x or z1.y != z2.y;
}

Complex operator+(const Complex &z1, const Complex &z2) {
    Complex sum(z1.x + z2.x, z1.y + z2.y);
    return sum;
}

void operator+=(Complex &z1, const Complex &z2) {
    z1 = z1 + z2;
}

Complex operator-(const Complex &z) {
    Complex neg(-z.x, -z.y);
    return neg;
}

Complex operator-(const Complex &z1, const Complex &z2) {
    Complex diff(z1.x - z2.x, z1.y - z2.y);
    return diff;
}

void operator-=(Complex &z1, const Complex &z2) {
    z1 = z1 - z2;
}

Complex operator*(const Complex &z1, const Complex &z2) {
    Complex prod(z1.r * z2.r, z1.A + z2.A, false);
    return prod;
}

void operator*=(Complex &z1, const Complex &z2) {
    z1 = z1 * z2;
}

Complex operator/(const Complex &z1, const Complex &z2) {
    if (z2 == "0") {
        throw runtime_error("/*** DIVISION BY 0 ***/");
    }
    return (z1.r / z2.r) * cis(z1.A - z2.A);
}

void operator/=(Complex &z1, const Complex &z2) {
    z1 = z1 / z2;
}

double dist(const Complex &z1, const Complex &z2) {
    return sqrt(pow(z1.x - z2.x, 2) + pow(z1.y - z2.y, 2));
}

Complex sect(const Complex &z1, const Complex &z2, const double &k) {
    Complex sect((z1.x + (k * z2.x)) / (k + 1), (z1.y + (k * z2.y)) / (k + 1));
    return sect;
}

Complex pow(const Complex &z1, const Complex &z2) {
    switch(z2.quadrant())
    {
        case -4:
        case -2: return cis(log(pow(z1.r, z2.y))) * cis(-z1.A * z2.y);
        case -3:
        case -1: return pow(z1.r, z2.x) * cis(z1.A * z2.x);
        case 0: return 1;
        case 1:
        case 2:
        case 3:
        case 4: return pow(z1.r, z2.x) * cis(z1.A * z2.x) * cis(log(pow(z1.r, z2.y))) * cis(-z1.A * z2.y);
        default: return "0";
    }
}

double principal_arg(const Complex &z) {
    return z.theta();
}

Complex conjugate(const Complex &z) {
    Complex bar(z.x, -z.y);
    return bar;
}

bool is_real(const Complex &z) {
    return z.y == 0;
}

Complex cross_ratio(const Complex &z, const Complex &z1, const Complex &z2, const Complex &z3) {
    return ((z - z1) / (z - z3)) / ((z2 - z1) / (z2 - z3));
}

Complex rotate(const Complex &z, const Complex &about, const double &through,
            const double &fin_len, const bool &anticlockwise = true) {
    Complex vec1 = z - about;
    Complex offset;
    if (anticlockwise) {
        offset.set(fin_len / vec1.r, through, false);
    } else {
        offset.set(fin_len / vec1.r, -through, false);
    }
    return about + (vec1 * offset);
}

Complex exp(const Complex &z) {
    return exp(z.x) * cis(z.y);
}

Complex sin(const Complex &z) {
    Complex i = "i";
    return (1 / (2 * i)) * (exp(i * z) - exp(-i * z));
}

Complex cos(const Complex &z) {
    Complex i = "i";
    return (exp(i * z) + exp(-i * z)) / 2;
}

Complex tan(const Complex &z) {
    return sin(z) / cos(z);
}

Complex cosec(const Complex &z) {
    return 1 / sin(z);
}

Complex sec(const Complex &z) {
    return 1 / cos(z);
}

Complex cot(const Complex &z) {
    return cos(z) / sin(z);
}

Complex sinh(const Complex &z) {
    return (exp(z) - exp(-z)) / 2;
}

Complex cosh(const Complex &z) {
    return (exp(z) + exp(-z)) / 2;
}

Complex tanh(const Complex &z) {
    return sinh(z) / cosh(z);
}

Complex cosech(const Complex &z) {
    return 1 / sinh(z);
}

Complex sech(const Complex &z) {
    return 1 / cosh(z);
}

Complex coth(const Complex &z) {
    return cosh(z) / sinh(z);
}

void Complex::print(ostream &out, const bool &cis) const {
    if (r == 0) {
        out << "0";
    } else if (A == 0) {
        out << setprecision(4) << r;
    } else if (not cis) {
        if (r == 1) {
            out << "e^(" << setprecision(4) << A << ")i";
        } else if (r == int(r)) {
            if (A == pi) {
                out << int(r) << "e^πi";
            } else if (A == -pi) {
                out << int(r) << "e^-πi";
            } else if (A / pi == int(A / pi)) {
                out << int(r) << "e^" << int(A / pi) << "πi";
            } else {
                out << int(r) << "e^(" << setprecision(4) << A << ")i";
            }
        } else {
            if (A == pi) {
                out << "(" << setprecision(4) << r << ")e^πi";
            } else if (A == -pi) {
                out << "(" << setprecision(4) << r << ")e^-πi";
            } else if (A / pi == int(A / pi)) {
                out << "(" << setprecision(4) << r << ")e^" << int(A / pi) << "πi";
            } else {
                out << "(" << setprecision(4) << r << ")e^(" << setprecision(4) << A << ")i";
            }
        }
    } else {
        if (r == 1) {
            out << "cis(" << setprecision(4) << A << ")";
        } else if (r == int(r)) {
            if (A == pi) {
                out << int(r) << "cis(π)";
            } else if (A == -pi) {
                out << int(r) << "cis(-π)";
            } else if (A / pi == int(A / pi)) {
                out << int(r) << "cis(" << int(A / pi) << "π)";
            } else {
                out << int(r) << "cis(" << setprecision(4) << A << ")";
            }
        } else {
            if (A == pi) {
                out << "(" << setprecision(4) << r << ")cis(π)";
            } else if (A == -pi) {
                out << "(" << setprecision(4) << r << ")cis(-π)";
            } else if (A / pi == int(A / pi)) {
                out << "(" << setprecision(4) << r << ")cis(" << int(A / pi) << "π)";
            } else {
                out << "(" << setprecision(4) << r << ")cis(" << setprecision(4) << A << ")";
            }
        }
    }
}

ostream &operator<<(ostream &out, const Complex &z1) {
    switch(z1.quadrant())
    {
        case -4:
        case -2:
            if (z1.y == -1) out << "-i";
            else if (z1.y == 1) out << "i";
            else out << setprecision(4) << z1.y << "i";
            break;
        case -3:
        case -1:
            out << setprecision(4) << z1.x;
            break;
        case 0:
            out << "0";
            break;
        case 1:
        case 2:
            if (z1.y == 1) out << setprecision(4) << z1.x << " + i";
            else out << setprecision(4) << z1.x << " + " << setprecision(4) << z1.y << "i";
            break;
        case 3:
        case 4:
            if (z1.y == -1) out << setprecision(4) << z1.x << " - i";
            else out << setprecision(4) << z1.x << " - " << setprecision(4) << -z1.y << "i";
    }
    return out;
}

istream &operator>>(istream &in, Complex &z1) {
    string num;
    getline(in, num);
    z1.extract_from_string(num);
    return in;
}