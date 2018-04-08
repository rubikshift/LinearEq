#include <iostream>
#include <fstream>
#include "solvingmethods.h"
#include "matrix.h"

int main()
{
	std::cout << "Rozpoczynam zadanie B..." << std::endl;
	{
		auto A = initMatrixA();
		auto b = initVecb();

		auto j = benchmark(A, b, jacobi);
		auto gs = benchmark(A, b, gaussSeidl);
		auto lu = benchmark(A, b, LUdecomposition);

		auto zadB = std::ofstream("zadB.txt", std::ios::out);
		zadB << "Jacobi:" << std::endl;
		zadB << "\tIteracje:" << std::get<0>(j) << "; czas: " << std::get<1>(j).count() << "; norma residuum: " << std::get<2>(j) << std::endl;
		zadB << "Gauss-Seidl:" << std::endl;
		zadB << "\tIteracje:" << std::get<0>(gs) << "; czas: " << std::get<1>(gs).count() << "; norma residuum: " << std::get<2>(gs) << std::endl;
		zadB.close();
	}

	std::cout << "Skonczylem zadanie B..." << std::endl;


	//THIS PART DOES NOT HAVE STOP CONDITION!
	//BE CAREFUL!
	std::cout << "Rozpoczynam zadanie C..." << std::endl;
	
	{
		auto A = initMatrixA(defaultN, 3, -1, -1);
		auto b = initVecb();
		auto j = benchmark(A, b, jacobi);
		auto gs = benchmark(A, b, gaussSeidl);

		auto zadC = std::ofstream("zadC.txt", std::ios::out);
		zadC << "Jacobi:" << std::endl;
		zadC << "\tIteracje:" << std::get<0>(j) << "; czas: " << std::get<1>(j).count() << "; norma residuum: " << std::get<2>(j) << std::endl;
		zadC << "Gauss-Seidl:" << std::endl;
		zadC << "\tIteracje:" << std::get<0>(gs) << "; czas: " << std::get<1>(gs).count() << "; norma residuum: " << std::get<2>(gs) << std::endl;
		zadC.close();
	}

	std::cout << "Skonczylem zadanie C..." << std::endl;
	//YOU HAVE BEEN WARNED

	std::cout << "Rozpoczynam zadanie D..." << std::endl;

	{
		auto A = initMatrixA(defaultN, 3, -1, -1);
		auto b = initVecb();
		auto lu = benchmark(A, b, LUdecomposition);
		auto zadD = std::ofstream("zadD.txt", std::ios::out);
		zadD << "Faktoryzacja LU:" << std::endl;
		zadD << "Czas: " << std::get<1>(lu).count() << "; norma residuum: " << std::get<2>(lu) << std::endl;
		zadD.close();
	}

	std::cout << "Skonczylem zadanie D..." << std::endl;
	std::cout << "Rozpoczynam zadanie E..." << std::endl;

	{
		const std::size_t N[] = { 100, 200, 400, 800, 1600, 3200, 6400, 12800 };
		auto zadE = std::ofstream("zadE.txt", std::ios::out);
		for (std::size_t i = 0; i < sizeof(N) / sizeof(std::size_t); i++)
		{
			std::cout << "N = " << N[i] << std::endl;
			zadE << "N: " << N[i] << std::endl;
			auto A = initMatrixA(N[i]);
			auto b = initVecb(N[i]);
			std::cout << "Jacobi..." << std::endl;
			auto j = benchmark(A, b, jacobi);
			zadE << "\tJacobi: " << std::get<1>(j).count() << std::endl;
			std::cout << "Gauss-Seidl..." << std::endl;
			auto gs = benchmark(A, b, gaussSeidl);
			zadE << "\tGauss-Seidl: " << std::get<1>(gs).count() << std::endl;
			std::cout << "LU..." << std::endl;
			auto lu = benchmark(A, b, LUdecomposition);
			zadE << "\tLU: " << std::get<1>(lu).count() << std::endl;
		}
		zadE.close();
	}
	
	return 0;
}