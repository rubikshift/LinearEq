#pragma once
#include <tuple>
#include <chrono>
#include <functional>
#include "matrix.h"

typedef std::tuple<unsigned int, std::chrono::duration<double>, double> BenchmarkResult;


BenchmarkResult benchmark(
	std::shared_ptr<Matrix> A,
	std::shared_ptr<Matrix> b,
	std::function<unsigned int(const std::shared_ptr<Matrix>&, const std::shared_ptr<Matrix>&, std::shared_ptr<Matrix>&)> solvingMethod
);

unsigned int jacobi(
	const std::shared_ptr<Matrix>& A,
	const std::shared_ptr<Matrix>& b,
	std::shared_ptr<Matrix>& x
);

unsigned int gaussSeidl(
	const std::shared_ptr<Matrix>& A,
	const std::shared_ptr<Matrix>& b,
	std::shared_ptr<Matrix>& x
);

unsigned int LUdecomposition(
	const std::shared_ptr<Matrix>& A,
	const std::shared_ptr<Matrix>& b,
	std::shared_ptr<Matrix>& x
);
