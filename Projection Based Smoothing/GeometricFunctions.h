#ifndef GEOMETRIC_FUNCTIONS
#define GEOMETRIC_FUNCTIONS

#include "LineSegment.h"

typedef Vertex Point;
typedef Vertex Vector;
typedef LineSegment Edge;
typedef LineSegment Normal;
typedef double Angle;

#define PI 3.14159265

class GeometricFunctions {
public:
	
	// Calculates the dot product value of two vectors
	static double dotProduct(const LineSegment &a, const LineSegment &b);

	// Calculates the cross product value of two vectors
	static double crossProduct(const LineSegment &a, const LineSegment &b);

	// Checks whether the two line segments intersect or not
	static bool checkLineIntersection(const Edge &a, const Edge &b);

	// Finds the intersection point of two line segments
	const static Vertex findLineIntersection(const Edge &a, const Edge &b);

	// Calculates the normal vector of two vectors using the cross product
	const static Normal findNormal(const LineSegment &a, const LineSegment &b);

	// Calculates the angle between two vectors
	static Angle calculateAngle(const LineSegment &a, const LineSegment &b);

	// Converts degrees to radians
	static Angle degreesToRadians(const Angle angle);

	// Converts radians to degrees
	static Angle radiansToDegrees(const Angle angle);

	// Normalizes the given vector to a unit vector
	const static LineSegment normalizeVector(const LineSegment &vector);

	// Converts a given line segment to a vector of lengths in each coordinate
	const static Vector convertLineSegment(const LineSegment &lineSegment);

};

#endif
