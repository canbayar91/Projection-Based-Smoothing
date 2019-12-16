#include "LineSegment.h"
#include <cmath>

LineSegment::LineSegment(const Vertex &start, const Vertex &end) {
	this->start = start;
	this->end = end;
}

LineSegment::LineSegment(const LineSegment& lineSegment) {
	this->start = lineSegment.start;
	this->end = lineSegment.end;
}

const double LineSegment::getLength() const {
	double componentX = pow(getLengthX(), 2);
	double componentY = pow(getLengthY(), 2);
	double componentZ = pow(getLengthZ(), 2);
	return sqrt(componentX + componentY + componentZ);
}

const double LineSegment::getProductX() const {
	return end.x - start.x;
}

const double LineSegment::getProductY() const {
	return end.y - start.y;
}

const double LineSegment::getProductZ() const {
	return end.z - start.z;
}

const double LineSegment::getLengthX() const {
	return abs(getProductX());
}

const double LineSegment::getLengthY() const {
	return abs(getProductY());
}

const double LineSegment::getLengthZ() const {
	return abs(getProductZ());
}
