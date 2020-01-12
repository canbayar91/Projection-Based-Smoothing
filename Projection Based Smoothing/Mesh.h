#ifndef MESH
#define MESH

#include "Face.h"
#include "GeometricFunctions.h"
#include <vector>
#include <map>

class Mesh {
public:

	Mesh();
	~Mesh();

	void addVertex(unsigned int id, Vertex coordinates);
	void addFace(unsigned int id, Quadrilateral* quadrilateral);

	size_t getVertexCount();
	size_t getFaceCount();

	NeighborhoodVertex* getVertex(size_t index);
	Face* getFace(size_t index);

	void smooth();
	
	void setEdgeCount(unsigned int edgeCount);
	unsigned int getEdgeCount();

private:

	const double INITIAL_STEP_SIZE = 1;

	std::vector<NeighborhoodVertex*> vertexList;
	std::vector<Face*> faceList;

	std::map<unsigned int, std::vector<unsigned int>> vertexMapping;
	std::map<unsigned int, std::vector<unsigned int>> faceMapping;

	unsigned int edgeCount;

	void processVertex(unsigned int index);
	void lineSearch(unsigned int index, Vector optimizationPath);

	double calculateConditionNumber(unsigned int index);

};

#endif
