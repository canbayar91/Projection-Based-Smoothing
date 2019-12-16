#ifndef DERIVATION_CALCULATOR
#define DERIVATION_CALCULATOR

#include "Quadrilateral.h"
#include "GeometricFunctions.h"

static const Angle THETA = 90;

enum AxisType {
	X_AXIS,
	Y_AXIS,
	Z_AXIS
};

class DerivationCalculator {
public:

	DerivationCalculator(NeighborhoodVertex* vertex, Quadrilateral* face);

	Vector findUntanglingVector();
	Vector findImprovementVector();

private:

	NeighborhoodVertex* vertex;
	Quadrilateral* face;

	double calculateCosDifference();
	double calculateJacobianDeterminant();
	double calculateFrobeniusNorm();

	double deriveAlpha(double jacobianDeterminant, double cosDifference, AxisType axis);

	LineSegment deriveNormal(AxisType axis);

	double deriveJacobianDeterminant(AxisType axis);
	double deriveFrobeniusNormSquare(AxisType axis);

	double deriveConditionNumber(double jacobianDeterminant, double frobeniusNorm, AxisType axis);

	double deriveEpsilon(AxisType axis);
	double deriveEta(AxisType axis);

	double heronsFormula(double l1, double l2, double l3);
};

#endif
