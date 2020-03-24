#ifndef TRIANGLE
#define TRIANGLE

#include "Primitive.h"

class Triangle : public Primitive {
public:

	const static int POINT_COUNT = 3;

	NeighborhoodVertex* a;
	NeighborhoodVertex* b;
	NeighborhoodVertex* c;

	Triangle(NeighborhoodVertex* a, NeighborhoodVertex* b, NeighborhoodVertex* c);
	Triangle(Vertex a, Vertex b, Vertex c);
	~Triangle();
	
	void output(std::string filename) const;

};

#endif
