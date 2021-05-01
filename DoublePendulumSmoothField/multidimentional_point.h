#ifndef _MULTIDIMENTIONAL_POINT_H_
#define _MULTIDIMENTIONAL_POINT_H_

#pragma once

#include <cstdio>
#include <array>
#include <vector>
#include <cmath>
#include <ostream>
#include <iomanip>
#include <initializer_list>

namespace dixelu {

template<size_t dims>
struct point {
	std::array<double, dims> pt;
	point() {
		for (int i = 0; i < dims; i++)
			pt[i] = 0;
	}
	explicit point(double v){
		for (int i = 0; i < dims; i++)
			pt[i] = v;
	}
	point(const point<dims>& P) {
		for (int i = 0; i < dims; i++)
			pt[i] = P.pt[i];
	}
	point(const std::initializer_list<double>& il_d) {
		auto y = il_d.begin();
		for (size_t i = 0; i < dims && y != il_d.end(); i++, y++)
			pt[i] = *y;
	}
	point(const std::initializer_list<int>& il) {
		auto y = il.begin();
		for (size_t i = 0; i < dims && y != il.end(); i++, y++)
			pt[i] = *y;
	}
	point(const std::vector<double>& il_d) {
		auto y = il_d.cbegin();
		for (size_t i = 0; i < dims && y != il_d.end(); i++, y++)
			pt[i] = *y;
	}
	point(const std::vector<int>& il) {
		auto y = il.cbegin();
		for (size_t i = 0; i < dims && y != il.end(); i++, y++)
			pt[i] = *y;
	}
	inline double get_norm2() const {
		double sum = 0;
		for (int i = 0; i < dims; i++)
			sum += pt[i] * pt[i];
		return sum;
	}
	inline double get_norm(double deg) const {
		double sum = 0;
		for (int i = 0; i < dims; i++)
			sum += std::pow(pt[i],deg);
		return std::pow(sum,1./deg);
	}
	inline double get_norm() const {
		return std::sqrt(get_norm2());
	}
	inline size_t get_dims() const {
		return dims;
	}
	inline double& operator[](size_t D) {
		return pt[D];
	}
	inline const double& operator[](size_t D) const {
		return pt[D];
	}
	inline void swap(point<dims>& P) {
		pt.swap(P.pt);
	}
	inline point<dims> operator+(const point<dims>& P) const {
		point<dims> N;
		for (size_t i = 0; i < dims; i++)
			N[i] = pt[i] + P[i];
		return N;
	}
	inline point<dims> operator+=(const point<dims>& P) {
		for (size_t i = 0; i < dims; i++)
			pt[i] += P[i];
		return *this;
	}
	inline point<dims> operator-(const point<dims>& P) const {
		point<dims> N;
		for (size_t i = 0; i < dims; i++)
			N[i] = pt[i] - P[i];
		return N;
	}
	inline point<dims> operator-=(const point<dims>& P) {
		for (size_t i = 0; i < dims; i++)
			pt[i] -= P[i];
		return *this;
	}
	inline point<dims> operator*(double M) const {
		point<dims> N;
		for (size_t i = 0; i < dims; i++)
			N[i] = pt[i] * M;
		return N;
	}
	inline point<dims> operator*=(double M) {
		for (size_t i = 0; i < dims; i++)
			pt[i] *= M;
		return *this;
	}
	inline point<dims> operator/(double M) const {
		point<dims> N;
		for (size_t i = 0; i < dims; i++)
			N[i] = pt[i] / M;
		return N;
	}
	inline point<dims> operator/=(double M) {
		return ((*this) *= (1. / M));
	}
	inline double operator*(const point<dims>& P) const {
		double sum = 0;
		for (size_t i = 0; i < dims; i++)
			sum += pt[i] * P[i];
		return sum;
	}
	inline point<dims> operator-() {
		for (size_t i = 0; i < dims; i++)
			pt[i] = 0 - pt[i];
		return *this;
	}
	inline bool operator<(const point<dims>& P) const {
		for (size_t i = 0; i < dims; i++)
			if (pt[i] >= P.pt[i])
				return false;
		return true;
	}
	inline bool operator>=(const point<dims>& P) const {
		return !(*this < P);
	}
	inline bool operator>(const point<dims>& P) const {
		return P < *this;
	}
	inline bool operator==(const point<dims>& P) const {
		for (size_t i = 0; i < dims; i++)
			if (pt[i] != P[i]) return false;
		return true;
	}
	inline bool operator!=(const point<dims>& P) const {
		return !(*this == P);
	}
	inline bool operator<=(const point<dims>& P) const {
		return !(P < *this);
	}
	inline point<dims> normalize() const {
		return (*this) / get_norm();
	}
};

template<size_t dims>
std::ostream& operator<<(std::ostream& os, const point<dims>& P) {
	os << "(";
	for (size_t i = 0; i < dims; i++) {
		os << P[i];
		if (i != dims - 1)
			os << ",";
	}
	os << ")";
	return os;
}

template<size_t dims>
inline point<dims> operator*(double M, point<dims> P) {
	point<dims> N;
	for (size_t i = 0; i < dims; i++)
		N[i] = P[i] * M;
	return N;
}

template<size_t dims>
struct sq_matrix {
	static constexpr double DOUBLE_EPSILON = 0.0000000000000001;
	double utilisation = 0;
	std::array<point<dims>, dims> ar;
	sq_matrix() {
		for (size_t i = 0; i < dims; i++)
			ar[i] = point<dims>();
	}
	explicit sq_matrix(double E_num) {
		for (size_t i = 0; i < dims; i++) {
			ar[i] = point<dims>();
			ar[i][i] = E_num;
		}
	}
	sq_matrix(std::initializer_list<point<dims>> IL) {
		size_t id = 0;
		for (auto&& p : IL) {
			if (id == dims)
				break;
			ar[id] = p;
			id++;
		}
	}
	sq_matrix(const std::array<point<dims>, dims>& ar) :ar(ar) {}
	sq_matrix(const sq_matrix<dims>& m) :ar(m.ar) {}
	inline point<dims> operator*(const point<dims>& p) const {
		point<dims> T;
		for (int i = 0; i < dims; i++) {
			T[i] = ar[i] * p;
		}
		return T;
	}
	inline sq_matrix<dims> operator*(double num) const {
		sq_matrix<dims> T;
		for (int i = 0; i < dims; i++) {
			T[i] = ar[i]*num;
		}
		return T;
	}
	inline sq_matrix<dims> operator*=(double num) {
		for (int i = 0; i < dims; i++) 
			ar[i]*=num;
		return *this;
	}
	inline sq_matrix<dims> operator/(double num) const {
		return *this * (1./num);
	}
	inline sq_matrix<dims> operator/=(double num) {
		return ((*this)*=(1./num));
	}
	inline sq_matrix<dims> operator+(const sq_matrix<dims>& p) const {
		sq_matrix<dims> T;
		for (int i = 0; i < dims; i++) {
			T[i] = ar[i] + p[i];
		}
		return T;
	}
	inline sq_matrix<dims> operator-(const sq_matrix<dims>& p) const {
		sq_matrix<dims> T;
		for (int i = 0; i < dims; i++) {
			T[i] = ar[i] - p[i];
		}
		return T;
	}
	inline sq_matrix<dims> operator+=(const sq_matrix<dims>& p) {
		for (int i = 0; i < dims; i++) {
			ar[i] += p[i];
		}
		return *this;
	}
	inline sq_matrix<dims> operator-=(const sq_matrix<dims>& p) {
		for (int i = 0; i < dims; i++) {
			ar[i] -= p[i];
		}
		return *this;
	}
	inline const point<dims>& operator[](size_t i) const {
		return ar[i];
	}
	inline point<dims>& operator[](size_t i) {
		return ar[i];
	}
	inline double& at(size_t point_id, size_t coordinate) {
		if (point_id < dims && coordinate < dims) {
			return ar[point_id][coordinate];
		}
		else return utilisation;
	}
	inline const double& at(size_t point_id, size_t coordinate) const {
		if (point_id < dims && coordinate < dims) {
			return ar[point_id][coordinate];
		}
		else return utilisation;
	}
	inline sq_matrix<dims> operator*(const sq_matrix<dims>& M) const {
		sq_matrix<dims> P;
		for (size_t y = 0; y < dims; y++) {
			for (size_t x = 0; x < dims; x++) {
				for (size_t i = 0; i < dims; i++) {
					P[y][x] += ar[y][i] * M[i][x];
				}
			}
		}
		return P;
	}
	inline sq_matrix<dims> inverse() const {
		double max_value = 0, mul = 0;
		size_t id = 0;
		sq_matrix<dims> E(1), A(*this);
		for (size_t step = 0; step < dims; step++) {
			id = step;
			max_value = 0;
			for (size_t coid = step; coid < dims; coid++) {
				if (std::abs(max_value) < std::abs(A.at(coid, step))) {
					max_value = A.at(coid, step);
					id = coid;
				}
			}
			if (std::abs(max_value) <= DOUBLE_EPSILON)
				return sq_matrix<dims>();
			if (id != step) {
				std::swap(A[id], A[step]);
				std::swap(E[id], E[step]);
			}
			for (size_t sum_id = 0; sum_id < dims; sum_id++) {
				if (sum_id == step)
					continue;
				mul = A.at(sum_id, step) / A.at(step, step);
				A[sum_id] -= (A[step] * mul);
				E[sum_id] -= (E[step] * mul);
			}
			mul = A.at(step, step);
			A[step] /= mul;
			E[step] /= mul;
		}
		return E;
	}
	inline double determinant() const {
		double determ = 1;
		double temp = 0, max_value = 0, mul = 0;
		size_t id = 0;
		sq_matrix<dims> A(*this);
		for (size_t step = 0; step < dims; step++) {
			id = step;
			max_value = 0;
			for (size_t coid = step; coid < dims; coid++) {
				if (std::abs(max_value) < std::abs(A.at(coid, step))) {
					max_value = A.at(coid, step);
					id = coid;
				}
			}
			if (std::abs(max_value) <= DOUBLE_EPSILON)
				return 0;
			if (id != step)
				std::swap(A[id], A[step]);
			for (size_t sum_id = 0; sum_id < dims; sum_id++) {
				if (sum_id == step)
					continue;
				mul = A.at(sum_id, step) / A.at(step, step);
				A[sum_id] -= (A[step] * mul);
			}
			mul = A.at(step, step);
			A[step] /= mul;
			determ *= mul;
		}
		return determ;
	}
	inline static point<dims> solve_using_eulers_method(sq_matrix<dims> A, point<dims> P) {
		double max_value = 0, mul = 0;
		size_t id = 0;
		for (size_t step = 0; step < dims; step++) {
			id = step;
			max_value = 0;
			for (size_t coid = step; coid < dims; coid++) {
				if (std::abs(max_value) < std::abs(A.at(coid, step))) {
					max_value = A.at(coid, step);
					id = coid;
				}
			}
			if (std::abs(max_value) <= DOUBLE_EPSILON)
				return point<dims>();
			if (id != step) {
				std::swap(A[id], A[step]);
				std::swap(P[id], P[step]);
			}
			for (size_t sum_id = 0; sum_id < dims; sum_id++) {
				if (sum_id == step)
					continue;
				mul = A.at(sum_id, step) / A.at(step, step);
				A[sum_id] -= A[step] * mul;
				P[sum_id] -= P[step] * mul;
			}
			mul = A.at(step, step);
			A[step] /= mul;
			P[step] /= mul;
		}
		return P;
	}
	inline sq_matrix<dims> operator^(int degree) {
		bool inverse = false;
		if (degree < 0) {
			inverse = true;
			degree = -degree;
		}
		else if (!degree)
			return sq_matrix<dims>(1.);
		auto cur_matrix = sq_matrix<dims>(1.), deg_co_matrix = *this;
		while (degree) {
			switch (degree & 1) {
			case 1:
				cur_matrix = cur_matrix * deg_co_matrix;
				deg_co_matrix = deg_co_matrix * deg_co_matrix;
				break;
			case 0:
				deg_co_matrix = deg_co_matrix * deg_co_matrix;
				break;
			}
			degree >>= 1;
		}
		return inverse ? cur_matrix.inverse() : cur_matrix;
	}
	inline sq_matrix<dims - 1> minor_matrix(const size_t& x_minor, const size_t& y_minor) const {
		auto minor_index = [](size_t x, size_t y, size_t minor_x, size_t minor_y) -> std::pair<int64_t, int64_t> {
			if (x == minor_x)
				return { -1,-1 };
			if (y == minor_y)
				return { -1,-1 };
			if (x > minor_x)
				x -= 1;
			if (y > minor_y)
				y -= 1;
			return { x,y };
		};
		sq_matrix<dims - 1> M;
		for (size_t y = 0; y < dims; y++) {
			for (size_t x = 0; x < dims; x++) {
				auto mxy = minor_index(x, y, x_minor, y_minor);
				int64_t mx = mxy.first;
				int64_t my = mxy.second;
				if (mx < 0 || my < 0)
					continue;
				M.at(my, mx) = ar[y][x];
			}
		}
		return M;
	}
};

template<size_t dims>
inline point<dims> cross_prod(const std::array<point<dims>, dims - 1>& points) {
	if (dims > 1) {
		point<dims> answer;
		std::array<point<dims>, dims> mx;
		mx[0] = point<dims>(std::vector<double>(dims, 0.));
		std::copy(points.begin(), points.end(), mx.begin() + 1);
		sq_matrix<dims> M(mx);
		for (int i = 0; i < dims; i++) {
			answer[i] = M.minor_matrix(i, 0).determinant() * ((i & 1) ? (1.) : (-1.));
		}
		return answer;
	}
	else if (dims == 1)
		return { 0 };
	else
		return {};
}

template<size_t dims>
inline std::ostream& operator<<(std::ostream& in, const sq_matrix<dims>& M) {
	for (size_t y = 0; y < dims; y++) {
		for (size_t x = 0; x < dims; x++) {
			in << std::setfill(' ') << std::setw(15) << M.at(y, x) << " ";
		}
		in << "\n";
	}
	return in;
}

} // namespace dixelu

#endif