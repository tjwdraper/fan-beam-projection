#ifndef _FIELD_NEW_HPP_
#define _FIELD_NEW_HPP_

#include "coord2d.hpp"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace opticalflow {
    template <class T>
    class Field {
        public:
            // Constructors and deconstructors
            Field(dim dimin) : _dimin(dimin), 
                               _step(1, dimin.x), 
                               _size(dimin.x*dimin.y),
                               _field(new T[_size]) {}

            Field(const Field<T>& fin) : _dimin(fin.get_dimensions()), 
                                         _step(fin.get_step()), 
                                         _size(fin.get_size()),
                                         _field(new T[fin.get_size()]) {
                std::copy(fin.get_field(), fin.get_field() + _size, _field);
            }
            Field(Field<T>&& other) noexcept : _dimin(other.get_dimensions()),
                                               _step(other.get_step()),
                                               _size(other.get_size()),
                                               _field(other.get_field()) {
                other._field = nullptr;
            }
            ~Field() { delete[] _field; }

            // Getters and setters
            T* get_field() { return _field; }
            const T* get_field() const { return _field; }
            dim get_dimensions() const { return _dimin; }
            dim get_step() const { return _step; }
            std::size_t get_size() const { return _size; }



            const T& get_val(std::size_t idx) const {
                #ifdef OPTICALFLOW_DEBUG
                    Field::check_idx(idx);
                #endif
                return _field[idx];
            }

            const T& get_val(std::size_t i, std::size_t j) const {
                #ifdef OPTICALFLOW_DEBUG
                    Field::check_idx(i, j);
                #endif
                return _field[i * _step.x + j * _step.y];
            }

            void set_val(const T& val, std::size_t idx) {
                #ifdef OPTICALFLOW_DEBUG
                    Field::check_idx(idx);
                #endif
                _field[idx] = val;
            }

            void set_val(const T& val, std::size_t i, std::size_t j) {
                #ifdef OPTICALFLOW_DEBUG
                    Field::check_idx(i, j);
                #endif
                _field[i * _step.x + j * _step.y] = val;
            }

            void fill(T val) {
                for (std::size_t idx = 0; idx < _size; ++idx)
                    _field[idx] = val;
            }

            // Operator overloading
            Field<T>& operator=(const Field<T>& fin) {
                if (_dimin != fin.get_dimensions())
                    throw std::runtime_error("In Field<T>& operator=(const Field<T>&) dimensions of input and target do not match.");

                std::copy(fin.get_field(), fin.get_field() + _size, _field);
                return *this;
            }

            Field<T> operator+(const Field<T>& fin) const {
                Field<T> fout(*this);
                fout += fin;
                return fout;
            }

            Field<T> operator-(const Field<T>& fin) const {
                Field<T> fout(*this);
                fout -= fin;
                return fout;
            }

            Field<T>& operator+=(const Field<T>& fin) {
                if (_dimin != fin.get_dimensions())
                    throw std::runtime_error("In Field<T>& operator+=(const Field<T>&) dimensions of input and target do not match.");

                for (std::size_t idx = 0; idx < _size; ++idx) 
                    _field[idx] += fin.get_val(idx);
                return *this;
            }

            Field<T>& operator-=(const Field<T>& fin) {
                if (_dimin != fin.get_dimensions())
                    throw std::runtime_error("In Field<T>& operator-=(const Field<T>&) dimensions of input and target do not match.");

                for (std::size_t idx = 0; idx < _size; ++idx) 
                    _field[idx] -= fin.get_val(idx);
                return *this;
            }

            Field<T> operator*(double val) const {
                Field<T> fout(*this);
                fout *= val;
                return fout;
            }

            Field<T>& operator*=(double val) {
                for (std::size_t idx = 0; idx < _size; ++idx)
                    _field[idx] *= val;
                return *this;
            }

            Field<T> operator/(double val) const {
                if (val == 0.0)
                    throw std::runtime_error("In Field<T> operator/(double ) const, division by zero.");

                Field<T> fout(*this);
                fout /= val;
                return fout;
            }

            Field<T>& operator/=(double val) {
                if (val == 0.0)
                    throw std::runtime_error("In Field<T>& operator/=(double ), division by zero.");

                for (std::size_t idx = 0; idx < _size; ++idx)
                    _field[idx] /= val;
                return *this;
            }

            friend Field<T> operator*(double val, const Field<T>& fin) {
                return fin * val;
            }

        private:
            void check_idx(std::size_t i, std::size_t j) const {
                if (i >= _dimin.x || j >= _dimin.y)
                    throw std::runtime_error("In T Field::get_val(std::size_t, std::size_t) input indices out of bound.");
            }

            void check_idx(std::size_t idx) const {
                if (idx >= _size)
                    throw std::runtime_error("In T Field::get_val(std::size_t) input indices out of bound.");
            }

            const dim _dimin;
            const dim _step;
            const std::size_t _size;
            T* _field = nullptr;
    };

    using Image = Field<double>;

    namespace image {
        // TODO: change to mxArray* in future or separate from namespace
        inline void load_image(const double* vals, Image& image) {
            std::copy(vals, vals + image.get_size(), image.get_field());
        }

        // TODO: change to mxArray* in future or separate from namespace
        inline void save_image(double* vals, const Image& image) {
            std::copy(image.get_field(), image.get_field() + image.get_size(), vals);
        }

        inline double normsq(const Image& image) {
            double norm(0.0);
            for (std::size_t idx = 0; idx < image.get_size(); ++idx) 
                norm += std::pow(image.get_val(idx), 2);
            return norm;
        }

        inline double norm(const Image& image) {
            return std::sqrt(opticalflow::image::normsq(image));
        }

        inline double sum(const Image& image) {
            double sum(0.0);
            for (std::size_t idx = 0; idx < image.get_size(); ++idx)
                sum += image.get_val(idx);
            return sum;
        }

        inline double max(const Image& image) {
            double max(image.get_val(0));
            if (image.get_size() == 1)
                return max;

            for (std::size_t idx = 1; idx < image.get_size(); ++idx)
                if (image.get_val(idx) > max)
                    max = image.get_val(idx);
            return max;
        }

        inline double min(const Image& image) {
            double min(image.get_val(0));
            if (image.get_size() == 1)
                return min;

            for (std::size_t idx = 1; idx < image.get_size(); ++idx)
                if (image.get_val(idx) < min)
                    min = image.get_val(idx);
            return min;
        }

        inline void normalize(Image& image) {
            double low = opticalflow::image::min(image);
            double high = opticalflow::image::max(image);

            if (low == high) 
                throw std::runtime_error("In opticalflow::image::normalize(Image&) input.min() = input.max().");

            for (std::size_t idx = 0; idx < image.get_size(); ++idx) 
                image.set_val((image.get_val(idx)-low)/(high-low), idx);
        }

        inline double mse(const Image& image_a, const Image& image_b) {
            if (image_a.get_dimensions() != image_b.get_dimensions())
                throw std::runtime_error("In opticalflow::image::mse(const Image&, const Image&), input images don't share the same dimensions.");

            double val(0.0);
            for (std::size_t idx = 0; idx < image_a.get_size(); ++idx) 
                val += std::pow(image_a.get_val(idx) - image_b.get_val(idx), 2.0);
            return val / image_a.get_size();
        }
    }

    namespace gradients {
        template <typename T>
        inline T dx(const opticalflow::Field<T>& field, const std::size_t i, const std::size_t j) {
            const dim dimin = field.get_dimensions();
                if (dimin.x < 2)
                    throw std::runtime_error("In T gradients::partial_x(const Field<T>&, const std::size_t, const std::size_t), x-dimension must be at least 2.");

            if (i == dimin.x-1)
                return T(0.0);
            else
                return field.get_val(i+1, j) - field.get_val(i,j);
        }
        template <typename T>
        inline T dy(const opticalflow::Field<T>& field, const std::size_t i, const std::size_t j) {
            const dim dimin = field.get_dimensions();
                if (dimin.y < 2)
                    throw std::runtime_error("In T gradients::partial_x(const Field<T>&, const std::size_t, const std::size_t), x-dimension must be at least 2.");

            if (j == dimin.y-1)
                return T(0.0);
            else
                return field.get_val(i, j+1) - field.get_val(i,j);
        }
        template <typename T>
        inline T adjoint_dx(const opticalflow::Field<T>& field, const std::size_t i, const std::size_t j) {
            const dim dimin = field.get_dimensions();
                if (dimin.x < 2)
                    throw std::runtime_error("In T gradients::partial_x(const Field<T>&, const std::size_t, const std::size_t), x-dimension must be at least 2.");

            if (i == 0)
                return field.get_val(0,j);
            else if (i == dimin.x-1)
                return -field.get_val(dimin.x-2,j);
            else
                return field.get_val(i,j) - field.get_val(i-1,j);
        }
        template <typename T>
        inline T adjoint_dy(const opticalflow::Field<T>& field, const std::size_t i, const std::size_t j) {
            const dim dimin = field.get_dimensions();
                if (dimin.y < 2)
                    throw std::runtime_error("In T gradients::partial_x(const Field<T>&, const std::size_t, const std::size_t), x-dimension must be at least 2.");

            if (j == 0)
                return field.get_val(i,0);
            else if (j == dimin.y-1)
                return -field.get_val(i,dimin.y-2);
            else
                return field.get_val(i,j) - field.get_val(i,j-1);
        }
    }
}

#endif