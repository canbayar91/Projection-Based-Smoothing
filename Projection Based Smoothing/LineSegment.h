#ifndef LINE_SEGMENT
#define LINE_SEGMENT

#include "Vertex.h"

class LineSegment {
public:

	Vertex start;
	Vertex end;

	LineSegment(const Vertex &start, const Vertex &end);
	LineSegment(const LineSegment& lineSegment);

	const double getLength() const;

	const double getProductX() const;
	const double getProductY() const;
	const double getProductZ() const;

	const double getLengthX() const;
	const double getLengthY() const;
	const double getLengthZ() const;
	
};

#endif
