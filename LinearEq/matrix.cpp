#include <cmath>
#include <iomanip>
#include <limits>
#include "matrix.h"

Matrix::Matrix(std::size_t rows, std::size_t cols) : cols(cols), rows(rows), size(cols*rows)
{
	elements = new double[size];
}

Matrix::Matrix(std::size_t rows, std::size_t cols, double initValue) : cols(cols), rows(rows), size(cols*rows)
{
	elements = new double[size];

	for (std::size_t i = 0; i < size; i++)
		elements[i] = initValue;
}

Matrix::Matrix(const Matrix & src) : cols(src.cols), rows(src.rows), size(src.size)
{
	elements = new double[size];
	memcpy(elements, src.elements, sizeof(double) * size);
}

Matrix::Matrix(std::size_t rows, std::size_t cols, double initValues[]) : cols(cols), rows(rows), size(cols*rows)
{
	elements = new double[size];
	for (std::size_t i = 0; i < rows; i++)
		for (std::size_t j = 0; j < cols; j++)
			elements[j + i * cols] = initValues[j + i * cols];
}

Matrix::~Matrix()
{
	delete[] elements;
}

void Matrix::set(size_t row, size_t col, double value)
{
	if (col >= cols || row >= rows)
		throw std::out_of_range("Nie ma takiego elementu");
	elements[col + row * cols] = value;
}

void Matrix::set(size_t n, double value)
{
	if (n > size)
		throw std::out_of_range("Nie ma takiego elementu");

	elements[n] = value;
}

double Matrix::at(std::size_t row, std::size_t col) const
{
	if (col >= cols || row >= rows)
		throw std::out_of_range("Nie ma takiego elementu");

	return elements[col + row * cols];
}

double Matrix::at(std::size_t n) const
{
	if (n > size)
		throw std::out_of_range("Nie ma takiego elementu");

	return elements[n];
}

std::tuple<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> Matrix::getLU() const
{
	if (this->cols != this->rows)
		throw std::runtime_error("Nie mozna utworzyc macierzy L i U, macierz bazowa musi byc kwadratowa");
	
	std::size_t m = this->cols;
	
	auto L = std::make_shared<Matrix>(m, m, 0);
	for (std::size_t i = 0; i < m; i++)
		L->set(i, i, 1);
	
	auto U = std::make_shared<Matrix>(*this);

	for(std::size_t k = 0; k < m - 1; k++)
		for (std::size_t j = k + 1; j < m; j++)
		{
			L->set(j, k, U->at(j, k) / U->at(k, k));
			for (std::size_t i = k; i < m; i++)
				U->set(j, i, U->at(j, i) - L->at(j, k) * U->at(k, i));
		}
	
	return std::tuple<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>>(L, U);
}

double Matrix::norm() const
{
	double n = 0.0;
	for (std::size_t i = 0; i < this->size; i++)
		n += pow(elements[i], 2);
	return sqrt(n);
}

std::shared_ptr<Matrix> Matrix::allocateDotProduct(const std::shared_ptr<Matrix>& B) const
{
	if (this->cols != B->rows)
		throw std::invalid_argument("Wymiary podanych macierzy nie pozwalaja na mnozenie");

	return std::make_shared<Matrix>(this->rows, B->cols);
}

void Matrix::dot(const std::shared_ptr<Matrix>& B, std::shared_ptr<Matrix>& product) const
{
	if (this->cols != B->rows || this->rows != product->rows || B->cols != product->cols)
		throw std::invalid_argument("Wymiary podanych macierzy nie pozwalaja na mnozenie");

	double s;
	for (std::size_t i = 0; i < this->rows; i++)
		for (std::size_t j = 0; j < product->cols; j++)
		{
			s = 0;
			for (std::size_t x = 0; x < this->cols; x++)
				s += this->at(i, x) * B->at(x, j);
			product->set(i, j, s);
		}
}
void Matrix::subtract(const std::shared_ptr<Matrix>& B) throw(std::invalid_argument)
{
	if (this->cols != B->cols && this->rows != B->rows)
		throw std::invalid_argument("Wymiary podanych macierzy nie pozwalaja na odejmowanie");
	
	double v;
	for (std::size_t i = 0; i < size; i++)
	{
		v = B->at(i);
		elements[i] -= v;
	}
}

void Matrix::copy(const std::shared_ptr<Matrix>& src) const throw(std::invalid_argument)
{
	if(this->cols != src->cols || this->rows != src->rows)
		throw std::invalid_argument("Wymiary nie pozwalaja na kopiowanie");
	
	memcpy(elements, src->elements, sizeof(double) * size);

}

std::shared_ptr<Matrix> initMatrixA(std::size_t N, double a1, double a2, double a3)
{
	auto A = std::make_shared<Matrix>(N, N, 0);
	for (std::size_t i = 0; i < N; i++)
	{
		A->set(i, i, a1);

		if (i < N - 1)
		{
			A->set(i + 1, i, a2);
			A->set(i, i + 1, a2);
		}

		if (i < N - 2)
		{
			A->set(i + 2, i, a3);
			A->set(i, i + 2, a3);
		}
	}

	return A;
}

std::shared_ptr<Matrix> initVecb(std::size_t N)
{
	auto b = std::make_shared<Matrix>(N, 1);

	for (std::size_t i = 0; i < b->rows; i++)
		b->set(i, sin(i * (5 + 1)));

	return b;
}

std::ostream & operator<<(std::ostream & os, const std::shared_ptr<Matrix>& m)
{
	os << "[";
	for (std::size_t i = 0; i < m->rows; i++)
	{
		if (i != 0)
			os << " ";
		os << "[";
		for (std::size_t j = 0; j < m->cols; j++)
		{
			os << std::setw(6) << std::setprecision(4) << m->at(i, j);
			if(j != m->cols - 1 )
				os << ", ";

		}
		os << "]";
		if (i != m->rows - 1)
			os << ",\n";
	}
	os << "]\n";
	return os;
}
