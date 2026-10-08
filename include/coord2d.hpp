#ifndef _COORD2D_H_
#define _COORD2D_H_

#include <cmath>
#include <stdexcept>
#include <iostream>

template <class T>
class coord2d {
    public:
        T x{};
        T y{};

        // Constructors
        coord2d() = default;
        coord2d(const coord2d<T>& a) : x(a.x), y(a.y) { isnan(x,y); isinf(x,y); }
        coord2d(T x, T y) : x(x), y(y) { isnan(x,y); isinf(x,y); }
        coord2d(T a) : x(a), y(a) { isnan(x,y); isinf(x,y); }

        coord2d<T>& operator=(const T& a) {
            x = a;
            y = a;
            return *this;
        }

        // Arithmatic operations
        coord2d<T> operator+(const coord2d<T>& c) const { return coord2d<T>(this->x + c.x, this->y + c.y); }
        coord2d<T> operator+(T a) const { return coord2d<T>(this->x + a, this->y + a); }
        coord2d<T> operator-(const coord2d<T>& c) const { return coord2d<T>(this->x - c.x, this->y - c.y); }
        coord2d<T> operator-(T a) const { return coord2d<T>(this->x - a, this->y - a); }
        coord2d<T> operator*(T a) const { return coord2d<T>(this->x * a, this->y * a); }

        friend coord2d<T> operator+(T a, const coord2d<T>& b) { return coord2d<T>(b.x + a, b.y + a); }
        friend coord2d<T> operator*(T a, const coord2d<T>& b) { return coord2d<T>(b.x * a, b.y * a); }

        coord2d<T>& operator+=(const coord2d<T>& c) {
            this->x += c.x;
            this->y += c.y;
            return *this;
        }

        coord2d<T>& operator+=(T a) {
            this->x += a;
            this->y += a;
            return *this;
        }

        coord2d<T>& operator-=(const coord2d<T>& c) {
            this->x -= c.x;
            this->y -= c.y;
            return *this;
        }

        coord2d<T>& operator-=(T a) {
            this->x -= a;
            this->y -= a;
            return *this;
        }        

        coord2d<T>& operator*=(T a) {
            this->x *= a;
            this->y *= a;
            return *this;
        }

        // Division
        coord2d<T> operator/(const T& a) const {
            if (a == 0) {
                throw std::runtime_error("Divide by zero exception");
            }
            return coord2d<T>(this->x/a, this->y/a);
        }

        coord2d<T>& operator/=(const T& a) {
            if (a == 0) {
                throw std::runtime_error("Divide by zero exception");
            }
            this->x /= a;
            this->y /= a;
            return *this;
        }

        // Unary operations
        coord2d<T>& operator-() {
            this->x *= -1;
            this->y *= -1;
            return *this;
        }

        // Boolean operators        
        bool operator==(const coord2d<T>& c) const {
            return (this->x == c.x) && (this->y == c.y);
        }

        bool operator!=(const coord2d<T>& c) const {
            return (this->x != c.x) || (this->y != c.y);
        }

        bool operator<(const coord2d<T>& c) const {
            return normsq(*this) < normsq(c);
        }

        bool operator>(const coord2d<T>& c) const {
            return normsq(*this) > normsq(c);
        }

        bool operator<=(const coord2d<T>& c) const {
            return normsq(*this) <= normsq(c);
        }

        bool operator>=(const coord2d<T>& c) const {
            return norsq(*this) >= normsq(c);
        }

        // Stream
        friend std::ostream& operator<<(std::ostream& os, const coord2d<T>& c) {
            os << "(" << c.x << ", " << c.y << ")";
            return os;
        }

    private:
        void isnan(T x, T y) const {
            if (std::isnan(x) || std::isnan(y)) {
                throw std::runtime_error("Coord2d<T> is NaN");
            }
        }

        void isinf(T x, T y) const {
            if (std::isinf(x) || std::isinf(y))  {
                throw std::runtime_error("Coord2d<T> is Inf");
            }
        }

};

typedef coord2d<unsigned int> dim;
typedef coord2d<double> vector2d;

// Adding specific methods for vector2d
inline double dot(const vector2d& a, const vector2d& b) {
    return a.x * b.x + a.y * b.y;
}

inline double normsq(const vector2d& a) {
    return dot(a, a);
}

inline double norm(const vector2d& a) {
    return std::sqrt(normsq(a));
}

#endif