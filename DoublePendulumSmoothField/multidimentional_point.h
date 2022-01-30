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

template<size_t dims, typename general_float_type>
struct point {
	std::array<general_float_type, dims> pt;
	using self_type = point<dims, general_float_type>;
	point() {
		for (int i = 0; i < dims; i++)
			pt[i] = general_float_type();
	}
	explicit point(general_float_type v){
		for (int i = 0; i < dims; i++)
			pt[i] = v;
	}
	point(const self_type&& P) :
		pt(std::move(P.pt))
	{
	}
	point(const self_type& P)
	{
		for (int i = 0; i < dims; i++)
			pt[i] = P.pt[i];
	}
	point(const std::initializer_list<general_float_type>& il_d) {
		auto y = il_d.begin();
		for (size_t i = 0; i < dims && y != il_d.end(); i++, y++)
			pt[i] = *y;
	}
	point(const std::initializer_list<int>& il) {
		auto y = il.begin();
		for (size_t i = 0; i < dims && y != il.end(); i++, y++)
			pt[i] = *y;
	}
	point(const std::vector<general_float_type>& il_d) {
		auto y = il_d.cbegin();
		for (size_t i = 0; i < dims && y != il_d.end(); i++, y++)
			pt[i] = *y;
	}
	point(const std::vector<int>& il) {
		auto y = il.cbegin();
		for (size_t i = 0; i < dims && y != il.end(); i++, y++)
			pt[i] = *y;
	}
	inline general_float_type get_norm2() const {
		general_float_type sum = general_float_type();
		for (int i = 0; i < dims; i++)
			sum += pt[i] * pt[i];
		return sum;
	}
	inline general_float_type get_norm(general_float_type deg) const {
		general_float_type sum = general_float_type();
		for (int i = 0; i < dims; i++)
			sum += std::pow(pt[i],deg);
		return std::pow(sum, 1./deg);
	}
	inline general_float_type get_norm() const {
		return std::sqrt(get_norm2());
	}
	inline size_t get_dims() const {
		return dims;
	}
	inline general_float_type& operator[](size_t D) {
		return pt[D];
	}
	inline const general_float_type& operator[](size_t D) const {
		return pt[D];
	}
	inline void swap(self_type& P) {
		pt.swap(P.pt);
	}
	inline self_type operator+(const self_type& P) const {
		self_type N;
		for (size_t i = 0; i < dims; i++)
			N[i] = pt[i] + P[i];
		return N;
	}
	inline self_type operator+=(const self_type& P) {
		for (size_t i = 0; i < dims; i++)
			pt[i] += P[i];
		return *this;
	}
	inline self_type operator-(const self_type& P) const {
		self_type N;
		for (size_t i = 0; i < dims; i++)
			N[i] = pt[i] - P[i];
		return N;
	}
	inline self_type operator-=(const self_type& P) {
		for (size_t i = 0; i < dims; i++)
			pt[i] -= P[i];
		return *this;
	}
	inline self_type operator*(general_float_type M) const {
		self_type N;
		for (size_t i = 0; i < dims; i++)
			N[i] = pt[i] * M;
		return N;
	}
	inline self_type operator*=(general_float_type M) {
		for (size_t i = 0; i < dims; i++)
			pt[i] *= M;
		return *this;
	}
	inline self_type operator/(general_float_type M) const {
		self_type N;
		for (size_t i = 0; i < dims; i++)
			N[i] = pt[i] / M;
		return N;
	}
	inline self_type operator/=(general_float_type M) {
		return ((*this) *= (1. / M));
	}
	inline general_float_type operator*(const self_type& P) const {
		general_float_type sum = general_float_type();
		for (size_t i = 0; i < dims; i++)
			sum += pt[i] * P[i];
		return sum;
	}
	inline self_type operator-() {
		for (size_t i = 0; i < dims; i++)
			pt[i] = 0 - pt[i];
		return *this;
	}
	inline bool operator<(const self_type& P) const {
		for (size_t i = 0; i < dims; i++)
			if (pt[i] >= P.pt[i])
				return false;
		return true;
	}
	inline bool operator>=(const self_type& P) const {
		return !(*this < P);
	}
	inline bool operator>(const self_type& P) const {
		return P < *this;
	}
	inline bool operator==(const self_type& P) const {
		for (size_t i = 0; i < dims; i++)
			if (pt[i] != P[i]) return false;
		return true;
	}
	inline bool operator!=(const self_type& P) const {
		return !(*this == P);
	}
	inline bool operator<=(const self_type& P) const {
		return !(P < *this);
	}
	inline self_type normalize() const {
		return (*this) / get_norm();
	}
	inline self_type operator=(const self_type& M) {
		for (size_t i = 0; i < dims; i++)
			pt[i] *= M;
		return *this;
	}
	inline self_type operator=(const self_type&& M) {
		pt = std::move(M.pt);
		return *this;
	}
};

template<size_t dims, typename general_float_type>
std::ostream& operator<<(std::ostream& os, const point<dims, general_float_type>& P) {
	os << "(";
	for (size_t i = 0; i < dims; i++) {
		os << P[i];
		if (i != dims - 1)
			os << ",";
	}
	os << ")";
	return os;
}

template<size_t dims, typename general_float_type>
inline point<dims, general_float_type> operator*(general_float_type M, point<dims, general_float_type> P) {
	point<dims, general_float_type> N;
	for (size_t i = 0; i < dims; i++)
		N[i] = P[i] * M;
	return N;
}

template<size_t dims, typename general_float_type, general_float_type GFLOAT_EPSILON = 0.0001>
struct sq_matrix {
	general_float_type utilisation = general_float_type();

	using point_type = point<dims, general_float_type>;
	using self_type = sq_matrix<dims, general_float_type>;

	std::array<point_type, dims> ar;

	sq_matrix() {
		for (size_t i = 0; i < dims; i++)
			ar[i] = point_type();
	}
	explicit sq_matrix(general_float_type E_num) {
		for (size_t i = 0; i < dims; i++) {
			ar[i] = point_type();
			ar[i][i] = E_num;
		}
	}
	sq_matrix(std::initializer_list<point_type> IL) {
		size_t id = 0;
		for (auto&& p : IL) {
			if (id == dims)
				break;
			ar[id] = p;
			id++;
		}
	}
	sq_matrix(const std::array<point_type, dims>& ar) :ar(ar) {}
	sq_matrix(const self_type& m) :ar(m.ar) {}
	inline point_type operator*(const point_type& p) const {
		point_type T;
		for (int i = 0; i < dims; i++) {
			T[i] = ar[i] * p;
		}
		return T;
	}
	inline self_type operator*(general_float_type num) const {
		self_type T;
		for (int i = 0; i < dims; i++) {
			T[i] = ar[i]*num;
		}
		return T;
	}
	inline self_type operator*=(general_float_type num) {
		for (int i = 0; i < dims; i++) 
			ar[i]*=num;
		return *this;
	}
	inline self_type operator/(general_float_type num) const {
		return *this * (1./num);
	}
	inline self_type operator/=(general_float_type num) {
		return ((*this)*=(1./num));
	}
	inline self_type operator+(const self_type& p) const {
		self_type T;
		for (int i = 0; i < dims; i++) {
			T[i] = ar[i] + p[i];
		}
		return T;
	}
	inline self_type operator-(const self_type& p) const {
		self_type T;
		for (int i = 0; i < dims; i++) {
			T[i] = ar[i] - p[i];
		}
		return T;
	}
	inline self_type operator+=(const self_type& p) {
		for (int i = 0; i < dims; i++) {
			ar[i] += p[i];
		}
		return *this;
	}
	inline self_type operator-=(const self_type& p) {
		for (int i = 0; i < dims; i++) {
			ar[i] -= p[i];
		}
		return *this;
	}
	inline const point_type& operator[](size_t i) const {
		return ar[i];
	}
	inline point_type& operator[](size_t i) {
		return ar[i];
	}
	inline general_float_type& at(size_t point_id, size_t coordinate) {
		if (point_id < dims && coordinate < dims) {
			return ar[point_id][coordinate];
		}
		else return utilisation;
	}
	inline const general_float_type& at(size_t point_id, size_t coordinate) const {
		if (point_id < dims && coordinate < dims) {
			return ar[point_id][coordinate];
		}
		else return utilisation;
	}
	inline self_type operator*(const self_type& M) const {
		self_type P;
		for (size_t y = 0; y < dims; y++) {
			for (size_t x = 0; x < dims; x++) {
				for (size_t i = 0; i < dims; i++) {
					P[y][x] += ar[y][i] * M[i][x];
				}
			}
		}
		return P;
	}
	inline self_type inverse() const {
		general_float_type max_value = general_float_type();
		general_float_type mul = general_float_type();
		size_t id = 0;
		self_type E(1), A(*this);
		for (size_t step = 0; step < dims; step++) {
			id = step;
			max_value = 0;
			for (size_t coid = step; coid < dims; coid++) {
				if (std::abs(max_value) < std::abs(A.at(coid, step))) {
					max_value = A.at(coid, step);
					id = coid;
				}
			}
			if (std::abs(max_value) <= GFLOAT_EPSILON)
				return self_type();
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
	inline general_float_type determinant() const {
		general_float_type determ = 1;
		general_float_type temp = 0, max_value = 0, mul = 0;
		size_t id = 0;
		self_type A(*this);
		for (size_t step = 0; step < dims; step++) {
			id = step;
			max_value = 0;
			for (size_t coid = step; coid < dims; coid++) {
				if (std::abs(max_value) < std::abs(A.at(coid, step))) {
					max_value = A.at(coid, step);
					id = coid;
				}
			}
			if (std::abs(max_value) <= GFLOAT_EPSILON)
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
	inline static point_type solve_using_eulers_method(self_type A, point_type P) {
		general_float_type max_value = 0, mul = 0;
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
			if (std::abs(max_value) <= GFLOAT_EPSILON)
				return point_type();
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
	inline self_type operator^(int degree) {
		bool inverse = false;
		if (degree < 0) {
			inverse = true;
			degree = -degree;
		}
		else if (!degree)
			return self_type(1.);
		auto cur_matrix = self_type(1.), deg_co_matrix = *this;
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
	inline sq_matrix<dims - 1, general_float_type> minor_matrix(const size_t& x_minor, const size_t& y_minor) const {
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
		sq_matrix<dims - 1, general_float_type> M;
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

template<size_t dims, typename general_float_type, general_float_type GFLOAT_EPSILON = 0.0001>
inline point<dims, general_float_type> cross_prod(const std::array<point<dims, general_float_type>, dims - 1>& points) {
	if (dims > 1) {
		point<dims, general_float_type> answer;
		std::array<point<dims, general_float_type>, dims> mx;
		mx[0] = point<dims, general_float_type>(std::vector<general_float_type>(dims, 0));
		std::copy(points.begin(), points.end(), mx.begin() + 1);
		sq_matrix<dims, general_float_type, GFLOAT_EPSILON> M(mx);
		for (int i = 0; i < dims; i++) 
			answer[i] = M.minor_matrix(i, 0).determinant() * ((i & 1) ? (1.) : (-1.));
		return answer;
	}
	else if (dims == 1)
		return { 0 };
	else
		return {};
}

template<size_t dims, typename general_float_type, general_float_type GFLOAT_EPSILON = 0.0001>
inline std::ostream& operator<<(std::ostream& in, const sq_matrix<dims, general_float_type, GFLOAT_EPSILON>& M) {
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