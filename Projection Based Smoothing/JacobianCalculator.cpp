#include "JacobianCalculator.h"

double JacobianCalculator::calculateConditionNumber(const Triangle triangle) {

	// Calculate Jacobian determinant and its Frobenius norm
	double jacobianDeterminant = calculateJacobianDeterminant(triangle);
	double frobeniusNorm = calculateFrobeniusNorm(triangle);

	// Condition number formula is referenced from Shivanna paper
	return pow(frobeniusNorm, 2) / (2 * jacobianDeterminant);
}

double JacobianCalculator::calculateJacobianDeterminant(const Triangle triangle) {

	// Get the edge vectors on the projected triangle
	LineSegment left(triangle.a->coordinates, triangle.c->coordinates);
	LineSegment right(triangle.a->coordinates, triangle.b->coordinates);

	// Calculate the determinant value for x-y coordinates
	double determinant = left.getLengthX() * right.getLengthY() - right.getLengthX() * left.getLengthY();

	// Return the determinant value
	return determinant;
}

double JacobianCalculator::calculateFrobeniusNorm(const Triangle triangle) {

	// Get the edge vectors on the projected triangle
	LineSegment left(triangle.a->coordinates, triangle.c->coordinates);
	LineSegment right(triangle.a->coordinates, triangle.b->coordinates);

	// Calculate the sums of the squares of the absolute values in the Jacobian matrix
	double squareSum = pow(left.getLengthX(), 2) + pow(left.getLengthY(), 2) + pow(right.getLengthX(), 2) + pow(right.getLengthY(), 2);

	// Return the square root of the square sum
	return sqrt(squareSum);
}
