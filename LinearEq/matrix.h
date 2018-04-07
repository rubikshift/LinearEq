#pragma once
#include <cstddef>
#include <stdexcept>
#include <memory>
#include <tuple>
#include <iostream>

const std::size_t defaultN = 996;

class Matrix
{
	friend std::ostream& operator<<(std::ostream & os, const std::shared_ptr<Matrix>& m);
	public:
		Matrix(std::size_t rows, std::size_t cols);
		Matrix(std::size_t rows, std::size_t cols, double initValue);
		Matrix(const Matrix& src);
		Matrix(std::size_t rows, std::size_t cols, double initValues[]);
		~Matrix();

		const std::size_t cols, rows, size;

		void set(size_t row, size_t col, double value) throw(std::out_of_range);
		void set(size_t n, double value) throw(std::out_of_range);

		double at(size_t row, size_t col) const throw(std::out_of_range);
		double at(size_t n) const throw(std::out_of_range);

		std::tuple<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> getLU() const throw(std::runtime_error);

		double norm() const;

		std::shared_ptr<Matrix> allocateDotProduct(const std::shared_ptr<Matrix>& B) const throw(std::invalid_argument);
		void dot(const std::shared_ptr<Matrix>& B, std::shared_ptr<Matrix>& result) const throw(std::invalid_argument);
		void subtract(const std::shared_ptr<Matrix>& B) throw(std::invalid_argument);

		void copy(const std::shared_ptr<Matrix>& dst) const throw(std::invalid_argument);

	private:
		double* elements;
};

std::shared_ptr<Matrix> initMatrixA(std::size_t N = defaultN, double a1 = 5 + 5, double a2 = -1, double a3 = -1);
std::shared_ptr<Matrix> initVecb(std::size_t N = defaultN);

std::ostream& operator<<(std::ostream & os, const std::shared_ptr<Matrix>& m);

