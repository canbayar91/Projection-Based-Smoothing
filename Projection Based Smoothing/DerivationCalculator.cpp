#include "DerivationCalculator.h"
#include "JacobianCalculator.h"

DerivationCalculator::DerivationCalculator(NeighborhoodVertex* vertex, Quadrilateral* face) {
	this->vertex = vertex;
	this->face = face;
}

Vector DerivationCalculator::findUntanglingVector() {

	// Calculate Jacobian determinant
	double jacobianDeterminant = calculateJacobianDeterminant();

	// Calculate the cosinus difference between normals and theta
	double cosDifference = calculateCosDifference();

	// Calculate alpha and untangling function derivative
	double alpha = cosDifference * jacobianDeterminant;
	double functionDerivative = alpha / abs(alpha) - 1;

	// Calculate the derivative for each axis (Shivanna formula 17)
	double derivativeX = deriveAlpha(jacobianDeterminant, cosDifference, X_AXIS) * functionDerivative;
	double derivativeY = deriveAlpha(jacobianDeterminant, cosDifference, Y_AXIS) * functionDerivative;
	double derivativeZ = deriveAlpha(jacobianDeterminant, cosDifference, Z_AXIS) * functionDerivative;

	// Calculate the derivative of epsilon
	double epsilonX = deriveEpsilon(X_AXIS);
	double epsilonY = deriveEpsilon(Y_AXIS);
	double epsilonZ = deriveEpsilon(Z_AXIS);

	// Calculate the derivative of eta
	double etaX = deriveEta(X_AXIS);
	double etaY = deriveEta(Y_AXIS);
	double etaZ = deriveEta(Z_AXIS);

	// Both epsilon and eta is used in the final vector
	double factorX = epsilonX + etaX;
	double factorY = epsilonY + etaY;
	double factorZ = epsilonZ + etaZ;

	// Create the final untangling vector
	Vector untanglingVector(derivativeX * factorX, derivativeY * factorY, derivativeZ * factorZ);
	GeometricFunctions::normalizeVector(untanglingVector);

	// Return the vector
	return untanglingVector;
}

Vector DerivationCalculator::findImprovementVector() {

	// Calculate Jacobian determinant
	double jacobianDeterminant = calculateJacobianDeterminant();

	// Calculate its Frobenius norm
	double frobeniusNorm = calculateFrobeniusNorm();

	// Calculate the derivative for each axis (Shivanna formula 20)
	double derivativeX = deriveConditionNumber(jacobianDeterminant, frobeniusNorm, X_AXIS);
	double derivativeY = deriveConditionNumber(jacobianDeterminant, frobeniusNorm, Y_AXIS);
	double derivativeZ = deriveConditionNumber(jacobianDeterminant, frobeniusNorm, Z_AXIS);

	// Calculate the derivative of epsilon
	double epsilonX = deriveEpsilon(X_AXIS);
	double epsilonY = deriveEpsilon(Y_AXIS);
	double epsilonZ = deriveEpsilon(Z_AXIS);

	// Calculate the derivative of eta
	double etaX = deriveEta(X_AXIS);
	double etaY = deriveEta(Y_AXIS);
	double etaZ = deriveEta(Z_AXIS);

	// Both epsilon and eta is used in the final vector
	double factorX = epsilonX + etaX;
	double factorY = epsilonY + etaY;
	double factorZ = epsilonZ + etaZ;

	// Create the final improvement vector
	Vector improvementVector(derivativeX * factorX, derivativeY * factorY, derivativeZ * factorZ);
	GeometricFunctions::normalizeVector(improvementVector);

	// Return the vector
	return improvementVector;
}

double DerivationCalculator::calculateJacobianDeterminant() {

	// Find the node that stores the given vertex 
	Node* node = face->findNode(vertex);

	// Find its next and previous vertices
	NeighborhoodVertex* p4 = node->next->data;
	NeighborhoodVertex* p2 = node->prev->data;

	// Find the Jacobian determinant for P3 triangle
	const Triangle triangle(p2, vertex, p4);
	double determinant = JacobianCalculator::calculateJacobianDeterminant(triangle);

	// Return the determinant value
	return determinant;
}

double DerivationCalculator::calculateFrobeniusNorm() {

	// Find the node that stores the given vertex 
	Node* node = face->findNode(vertex);

	// Find its next and previous vertices
	NeighborhoodVertex* p4 = node->next->data;
	NeighborhoodVertex* p2 = node->prev->data;

	// Find the frobenius norm for P3 triangle
	const Triangle triangle(p2, vertex, p4);
	double frobeniusNorm = JacobianCalculator::calculateFrobeniusNorm(triangle);

	// Return the frobenius norm value
	return frobeniusNorm;
}

double DerivationCalculator::calculateCosDifference() {

	// Find the node that stores the given vertex 
	Node* node = face->findNode(vertex);

	// Find its next and previous vertices
	NeighborhoodVertex* p4 = node->next->data;
	NeighborhoodVertex* p2 = node->prev->data;

	// Find the normal of the triangle that includes the vertex
	const LineSegment CD(vertex->coordinates, p4->coordinates);
	const LineSegment CB(vertex->coordinates, p2->coordinates);
	const Normal normalP3 = GeometricFunctions::findNormal(CD, CB);

	// Find the remaining vertix
	NeighborhoodVertex* p1 = node->next->next->data;

	// Find the normal of the triangle that doesn't include the vertex
	const LineSegment AB(p1->coordinates, p2->coordinates);
	const LineSegment AD(p1->coordinates, p4->coordinates);
	const Normal normal = GeometricFunctions::findNormal(AB, AD);

	// Find dot product of normals (cos value since normals have unit length)
	double dotProduct = GeometricFunctions::dotProduct(normal, normalP3);

	// Return difference between cos angles
	return dotProduct - cos(THETA);
}

double DerivationCalculator::deriveAlpha(double jacobianDeterminant, double cosDifference, AxisType axis) {

	// Find the node that stores the given vertex 
	Node* node = face->findNode(vertex);

	// Find its next and previous vertices
	NeighborhoodVertex* p4 = node->next->data;
	NeighborhoodVertex* p2 = node->prev->data;
	NeighborhoodVertex* p1 = node->next->next->data;

	// Find the normal of the triangle that doesn't include the vertex
	const LineSegment AB(p1->coordinates, p2->coordinates);
	const LineSegment AD(p1->coordinates, p4->coordinates);
	const Normal normal = GeometricFunctions::findNormal(AB, AD);

	// Calculate the dot product of normals
	LineSegment normalDerivative = deriveNormal(axis);
	double dotProduct = GeometricFunctions::dotProduct(normal, normalDerivative);

	// Find the derivative of alpha depending on the selected axis (Shivanna formula 18)
	double jacobianDerivative = deriveJacobianDeterminant(axis);
	return dotProduct * jacobianDeterminant + cosDifference * jacobianDerivative;
}

LineSegment DerivationCalculator::deriveNormal(AxisType axis) {

	// Find the node that stores the given vertex 
	Node* node = face->findNode(vertex);

	// Find its next and previous vertices
	// Note: p1 is actually p4 (I think)
	NeighborhoodVertex* p1 = node->next->data;
	NeighborhoodVertex* p2 = node->prev->data;

	// Return the correct vector depending on the axis (Shivanna formula 19)
	if (axis == X_AXIS) {
		Vertex start(0, p1->coordinates.z, p2->coordinates.y);
		Vertex end(0, p2->coordinates.z, p1->coordinates.y);
		return LineSegment(start, end);
	} else if (axis == Y_AXIS) {
		Vertex start(p2->coordinates.z, 0, p2->coordinates.y);
		Vertex end(p1->coordinates.z, 0, p1->coordinates.y);
		return LineSegment(start, end);
	} 

	Vertex start(p1->coordinates.y, p2->coordinates.x, 0);
	Vertex end(p2->coordinates.y, p1->coordinates.x, 0);
	return LineSegment(start, end);
}

double DerivationCalculator::deriveJacobianDeterminant(AxisType axis) {

	// Find the node that stores the given vertex 
	Node* node = face->findNode(vertex);

	// Find its next and previous vertices
	// Note: p1 is actually p4 (I think)
	NeighborhoodVertex* p1 = node->next->data;
	NeighborhoodVertex* p2 = node->prev->data;

	// Find each line segment between points
	LineSegment ls1(p1->coordinates, p2->coordinates);
	LineSegment ls2(p2->coordinates, vertex->coordinates);
	LineSegment ls3(vertex->coordinates, p1->coordinates);

	// Get their lengths
	double l1 = ls1.getLength();
	double l2 = ls2.getLength();
	double l3 = ls3.getLength();

	// Calculate x and y (Shivanna formulas 1 and 2)
	double x = (pow(l3, 2) - pow(l2, 2) + pow(l1, 2)) / (2 * l1);
	double y = sqrt(abs(pow(l3, 2) - pow(x, 2)));

	// Coordinate differences depend on the selected axis
	double diff31;
	double diff21;
	if (axis == X_AXIS) {
		diff31 = vertex->coordinates.x - p1->coordinates.x;
		diff21 = p2->coordinates.x - p1->coordinates.x;
	} else if (axis == Y_AXIS) {
		diff31 = vertex->coordinates.y - p1->coordinates.y;
		diff21 = p2->coordinates.y - p1->coordinates.y;
	} else {
		diff31 = vertex->coordinates.z - p1->coordinates.z;
		diff21 = p2->coordinates.z - p1->coordinates.z;
	}

	// Calculate and return the derivative (Shivanna formula 21)
	return (l1 * diff31 - x * diff21) / y;
}

double DerivationCalculator::deriveFrobeniusNormSquare(AxisType axis) {

	// Find the node that stores the given vertex
	Node* node = face->findNode(vertex);

	// Find its next and previous vertices
	// Note: p1 is actually p4 (I think)
	NeighborhoodVertex* p1 = node->next->data;
	NeighborhoodVertex* p2 = node->prev->data;

	// Return value depends on the selected coordinate (Shivanna formula 22 - used the last one)
	if (axis == X_AXIS) {
		return 2 * (2 * vertex->coordinates.x - (p1->coordinates.x + p2->coordinates.x));
	}

	if (axis == Y_AXIS) {
		return 2 * (2 * vertex->coordinates.y - (p1->coordinates.y + p2->coordinates.y));
	}

	return 2 * (2 * vertex->coordinates.z - (p1->coordinates.z + p2->coordinates.z));
}

double DerivationCalculator::deriveConditionNumber(double jacobianDeterminant, double frobeniusNorm, AxisType axis) {

	// Derivation values will change depending on the axis
	double jacobianDeterminantDerivation = deriveJacobianDeterminant(axis);
	double frobeniusNormSquareDerivation = deriveFrobeniusNormSquare(axis);

	// Calculate the condition number derivation (Shivanna formula 20)
	double firstPart = jacobianDeterminant * frobeniusNormSquareDerivation;
	double secondPart = jacobianDeterminantDerivation * pow(frobeniusNorm, 2);
	double denominator = 2 * pow(jacobianDeterminant, 2);

	// Return the derivation
	return (firstPart + secondPart) / denominator;
}

double DerivationCalculator::deriveEpsilon(AxisType axis) {

	// Find the node that stores the given vertex 
	Node* node = face->findNode(vertex);

	// Find its next and previous vertices
	// Note: p1 is actually p4 (I think)
	NeighborhoodVertex* p1 = node->next->data;
	NeighborhoodVertex* p2 = node->prev->data;

	// Find each line segment between points
	LineSegment ls1(p1->coordinates, p2->coordinates);
	LineSegment ls2(p2->coordinates, vertex->coordinates);
	LineSegment ls3(vertex->coordinates, p1->coordinates);

	// Get their lengths
	double l1 = ls1.getLength();
	double l2 = ls2.getLength();
	double l3 = ls3.getLength();

	// Calculate x and y (Shivanna formulas 1 and 2)
	double x = (pow(l3, 2) - pow(l2, 2) + pow(l1, 2)) / (2 * l1);
	double y = sqrt(abs(pow(l3, 2) - pow(x, 2)));

	// Temporary variables to store eta values
	double n1 = 0;
	double n2 = 0;
	double n3 = y;

	// Use the x-y-z coordinate values depending on the axis input
	double numerator = 0;
	if (axis == X_AXIS) {
		numerator = (n1 - n2) * vertex->coordinates.x + (n2 - n3) * p1->coordinates.x + (n3 - n1) * p2->coordinates.x;
	} else if (axis == Y_AXIS) {
		numerator = (n1 - n2) * vertex->coordinates.y + (n2 - n3) * p1->coordinates.y + (n3 - n1) * p2->coordinates.y;
	} else {
		numerator = (n1 - n2) * vertex->coordinates.z + (n2 - n3) * p1->coordinates.z + (n3 - n1) * p2->coordinates.z;
	}

	// Calculate the area of the triangle
	double triangleArea = heronsFormula(l1, l2, l3);

	// Calculate and return the derivative value
	return (0.5 * numerator) / triangleArea;
}

double DerivationCalculator::deriveEta(AxisType axis) {

	// Find the node that stores the given vertex 
	Node* node = face->findNode(vertex);

	// Find its next and previous vertices
	// Note: p1 is actually p4 (I think)
	NeighborhoodVertex* p1 = node->next->data;
	NeighborhoodVertex* p2 = node->prev->data;

	// Find each line segment between points
	LineSegment ls1(p1->coordinates, p2->coordinates);
	LineSegment ls2(p2->coordinates, vertex->coordinates);
	LineSegment ls3(vertex->coordinates, p1->coordinates);

	// Get their lengths
	double l1 = ls1.getLength();
	double l2 = ls2.getLength();
	double l3 = ls3.getLength();

	// Calculate x and y (Shivanna formulas 1 and 2)
	double x = (pow(l3, 2) - pow(l2, 2) + pow(l1, 2)) / (2 * l1);
	double y = sqrt(abs(pow(l3, 2) - pow(x, 2)));

	// Temporary variables to store epsilon values
	double e1 = 0;
	double e2 = l1;
	double e3 = x;

	// Use the x-y-z coordinate values depending on the axis input
	double numerator = 0;
	if (axis == X_AXIS) {
		numerator = (e1 - e2) * vertex->coordinates.x + (e2 - e3) * p1->coordinates.x + (e3 - e1) * p2->coordinates.x;
	} else if (axis == Y_AXIS) {
		numerator = (e1 - e2) * vertex->coordinates.y + (e2 - e3) * p1->coordinates.y + (e3 - e1) * p2->coordinates.y;
	} else {
		numerator = (e1 - e2) * vertex->coordinates.z + (e2 - e3) * p1->coordinates.z + (e3 - e1) * p2->coordinates.z;
	}

	// Calculate the area of the triangle
	double triangleArea = heronsFormula(l1, l2, l3);

	// Calculate and return the derivative value
	return (0.5 * numerator) / triangleArea;
}

double DerivationCalculator::heronsFormula(double l1, double l2, double l3) {

	// Used Heron's formula to find the area of the triangle using the lengths of its sides
	double s = 0.5 * (l1 + l2 + l3);
	double area = sqrt(s * (s - l1) * (s - l2) * (s - l3));
	return area;
}
