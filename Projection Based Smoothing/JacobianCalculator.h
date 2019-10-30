#ifndef JACOBIAN_CALCULATOR
#define JACOBIAN_CALCULATOR

#include "ProjectionFunctions.h"

class JacobianCalculator {
public:

	static double calculateConditionNumber(const Triangle triangle);

private:

	static double calculateJacobianDeterminant(const Triangle triangle);
	static double calculateFrobeniusNorm(const Triangle triangle);

};

#endif
