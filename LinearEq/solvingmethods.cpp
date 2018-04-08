#include <iostream>
#include "solvingmethods.h"

//#define MONITOR_RESIDUUM

#ifdef MONITOR_RESIDUUM
	#include <fstream>
	const std::size_t maxIter = 1000;
	auto jacobiMonitor = std::ofstream("zadC_jacobi.txt", std::ios::out);
	auto gaussseidlMonitor = std::ofstream("zadC_gauss.txt", std::ios::out);
#endif 


const double STOP = 1. / 1000000000;

BenchmarkResult benchmark(
		std::shared_ptr<Matrix> A,
		std::shared_ptr<Matrix> b,
		std::function<unsigned int(const std::shared_ptr<Matrix>&, const std::shared_ptr<Matrix>&, std::shared_ptr<Matrix>&)> solvingMethod)
{
	auto x = std::make_shared<Matrix>(A->rows, b->cols, 1);

	auto start = std::chrono::system_clock::now();
	auto iter = solvingMethod(A, b, x);
	auto end = std::chrono::system_clock::now();
	auto r = A->allocateDotProduct(x);
	A->dot(x, r);
	r->subtract(b);
	std::chrono::duration<double> diff = end - start;
	return BenchmarkResult(iter, diff, r->norm());
}

unsigned int jacobi(
		const std::shared_ptr<Matrix>& A, 
		const std::shared_ptr<Matrix>& b, 
		std::shared_ptr<Matrix>& x)
{
	unsigned int iter = 0;
	
	auto oldX = std::make_shared<Matrix>(*x);
	auto r = A->allocateDotProduct(x);
	A->dot(x, r);
	r->subtract(b);

	double xi;
	double norm = r->norm();
	while (norm > STOP)
	{
		oldX->copy(x);
		
		for (std::size_t i = 0; i < x->rows; i++)
		{
			xi = b->at(i);
			for (std::size_t j = 0; j < x->rows; j++)
			{
				if (j == i)
					continue;
				xi = xi - (A->at(i, j) * oldX->at(j));
			}
			xi /= A->at(i, i);
			x->set(i, xi);
		}
		
		A->dot(x, r);
		r->subtract(b);
		iter++;
		norm = r->norm();

#ifdef MONITOR_RESIDUUM
		jacobiMonitor << norm << std::endl;
		if (iter == maxIter)
			break;
#endif

	}
#ifdef MONITOR_RESIDUUM
	jacobiMonitor.close();
#endif

	return iter;
}

unsigned int gaussSeidl(
		const std::shared_ptr<Matrix>& A, 
		const std::shared_ptr<Matrix>& b, 
		std::shared_ptr<Matrix>& x)
{
	unsigned int iter = 0;
	
	auto r = A->allocateDotProduct(x);
	A->dot(x, r);
	r->subtract(b);
	double xi;
	double norm = r->norm();
	while (norm > STOP)
	{
		for (std::size_t i = 0; i < x->rows; i++)
		{
			xi = b->at(i);
			for (std::size_t j = 0; j < x->rows; j++)
			{
				if (j == i)
					continue;
				xi -= A->at(i, j) * x->at(j);
			}
			xi /= A->at(i, i);
			x->set(i, xi);
		}
		A->dot(x, r);
		iter++;
		r->subtract(b);
		norm = r->norm();
#ifdef MONITOR_RESIDUUM
		gaussseidlMonitor << norm << std::endl;
		if (iter == maxIter)
			break;
#endif
	}
#ifdef MONITOR_RESIDUUM
	gaussseidlMonitor.close();
#endif

	return iter;
}

unsigned int LUdecomposition(
		const std::shared_ptr<Matrix>& A, 
		const std::shared_ptr<Matrix>& b, 
		std::shared_ptr<Matrix>& x)
{
	unsigned int iter = 1;
	auto LU = A->getLU();
	auto y = std::make_shared<Matrix>(*x);

	double yi;
	for (std::size_t i = 0; i < y->rows; i++)
	{
		yi = b->at(i);
		for (std::size_t j = 0; j < i; j++)
			yi -= std::get<0>(LU)->at(i, j) * y->at(j);
		//yi /= std::get<0>(LU)->at(i, i); L->at(i,i) always equals 1
		y->set(i, yi);
	}

	double xi;
	for (std::size_t i = 0; i < x->rows; i++)
	{
		xi = y->at(x->rows - 1 - i);
		for (std::size_t j = 0; j < i; j++)
			xi -= std::get<1>(LU)->at(x->rows - 1 - i, x->rows - 1 - j) * x->at(x->rows - 1 - j);
		xi /= std::get<1>(LU)->at(x->rows - 1 - i, x->rows - 1 - i);
		x->set(x->rows - 1 - i, xi);
	}

	return iter;
}