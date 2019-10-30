#include "Mesh.h"
#include "JacobianCalculator.h"
#include <iostream>

double max = 0;

Mesh::Mesh() {
	// Default constructor
}

Mesh::~Mesh() {
	for (size_t i = 0; i < vertexList.size(); i++) {
		delete vertexList[i];
	}
}

void Mesh::addVertex(unsigned int id, Vertex coordinates) {
	NeighborhoodVertex* vertex = new NeighborhoodVertex(id, coordinates);
	vertexList.push_back(vertex);
}

void Mesh::addFace(unsigned int id, Quadrilateral* quadrilateral) {

	NeighborhoodVertex* a = quadrilateral->a;
	NeighborhoodVertex* b = quadrilateral->b;
	NeighborhoodVertex* c = quadrilateral->c;
	NeighborhoodVertex* d = quadrilateral->d;

	vertexMapping[a->id].push_back(b->id);
	vertexMapping[b->id].push_back(c->id);
	vertexMapping[c->id].push_back(d->id);
	vertexMapping[d->id].push_back(a->id);

	vertexMapping[a->id].push_back(d->id);
	vertexMapping[b->id].push_back(a->id);
	vertexMapping[c->id].push_back(b->id);
	vertexMapping[d->id].push_back(c->id);

	Face* face = new Face(id, quadrilateral);
	faceList.push_back(face);

	faceMapping[a->id].push_back(id);
	faceMapping[b->id].push_back(id);
	faceMapping[c->id].push_back(id);
	faceMapping[d->id].push_back(id);
}

size_t Mesh::getVertexCount() {
	return vertexList.size();
}

size_t Mesh::getFaceCount() {
	return faceList.size();
}

NeighborhoodVertex* Mesh::getVertex(size_t index) {
	return vertexList[index];
}

Face* Mesh::getFace(size_t index) {
	return faceList[index];
}

void Mesh::setEdgeCount(unsigned int edgeCount) {
	this->edgeCount = edgeCount;
}

unsigned int Mesh::getEdgeCount() {
	return edgeCount;
}

void Mesh::smooth() {

	for (size_t i = 0; i < vertexList.size(); i++) {
		double conditionNumber = calculateConditionNumber(i);
		if (conditionNumber > max) {
			max = conditionNumber;
		}
	}

	std::cout << "Maximum Condition Number: " << max << std::endl;
}

double Mesh::calculateConditionNumber(unsigned int index) {

	// Total value of the Jacobian condition numbers of each adjacent face
	double totalValue = 0;

	// Get the vertex with the given index
	NeighborhoodVertex* current = vertexList[index];

	// Iterate through the adjacent faces of the given vertex
	std::vector<unsigned int> adjacentFaces = faceMapping[index];
	for (unsigned int faceId : adjacentFaces) {

		// For each face, get the quadrilateral data
		Quadrilateral* face = faceList[faceId]->getQuadrilateral();

		// Find the node that stores the given vertex 
		Node* node = face->findNode(current);
		if (node != NULL) {

			// Find its next and previous vertices
			NeighborhoodVertex* next = node->next->data;
			NeighborhoodVertex* previous = node->prev->data;

			// Form a triangle using the next and previuos vertices and align it to the z = 0 plane
			Triangle* triangle = new Triangle(current, next, previous);
			Triangle projected = ProjectionFunctions::projectTriangle(triangle);

			// Calculate the Jacobian condition number value for the triangle and add it to the total sum
			double conditionNumber = JacobianCalculator::calculateConditionNumber(projected);
			totalValue += conditionNumber;
		}
	}

	return totalValue;
}
