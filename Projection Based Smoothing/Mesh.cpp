#include "Mesh.h"
#include "JacobianCalculator.h"
#include "DerivationCalculator.h"
#include <iostream>

Mesh::Mesh() {
	// Default constructor
}

Mesh::~Mesh() {

	// Delete all created faces
	for (size_t i = 0; i < faceList.size(); i++) {
		delete faceList[i];
	}

	// Delete all created vertices
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

	// Iterate over each vertex and process them
	for (size_t i = 0; i < vertexList.size(); i++) {
		processVertex(i);
	}
}

void Mesh::processVertex(unsigned int index) {

	// Get the vertex with the given index
	NeighborhoodVertex* current = vertexList[index];

	// Iterate through the adjacent faces of the given vertex
	std::vector<unsigned int> adjacentFaces = faceMapping[index];
	for (unsigned int faceId : adjacentFaces) {

		// For each face, get the quadrilateral data
		Quadrilateral* face = faceList[faceId]->getQuadrilateral();

		// A new derivation calculator instance is created for each neighbor face
		DerivationCalculator derivationCalculator(current, face);

		// Find the node that stores the given vertex 
		Node* node = face->findNode(current);
		if (node != NULL) {

			// Find its next and previous vertices
			NeighborhoodVertex* p4 = node->next->data;
			NeighborhoodVertex* p2 = node->prev->data;

			// Find the normal of the triangle that includes the vertex
			const LineSegment CD(current->coordinates, p4->coordinates);
			const LineSegment CB(current->coordinates, p2->coordinates);
			const Normal normalP3 = GeometricFunctions::findNormal(CD, CB);

			// Find the remaining vertix
			NeighborhoodVertex* p1 = node->next->next->data;

			// Find the normal of the triangle that doesn't include the vertex
			const LineSegment AB(p1->coordinates, p2->coordinates);
			const LineSegment AD(p1->coordinates, p4->coordinates);
			const Normal normal = GeometricFunctions::findNormal(AB, AD);

			// Find the angle between triangle normals
			Angle beta = GeometricFunctions::calculateAngle(normalP3, normal);
			if (beta > THETA) {

				// Untangling vector is only calculated when the triangle normals show different directions
				Vector untanglingVector = derivationCalculator.findUntanglingVector();
				lineSearch(index, untanglingVector);
			}

			// Calculate the improvement vector along the neighbor face
			Vector improvementVector = derivationCalculator.findImprovementVector();
			lineSearch(index, improvementVector);
		}
	}
}

void Mesh::lineSearch(unsigned int index, Vector optimizationPath) {

	// Get the vertex with the given index
	NeighborhoodVertex* current = vertexList[index];
	Vertex originalCoordinates = current->coordinates;

	// Initial step size
	double stepSize = INITIAL_STEP_SIZE;

	 do {

		// Calculate the condition number before the improvement
		double originalConditionNumber = calculateConditionNumber(index);

		// Calculate the updated coordinates on the optimization path
		Vertex updatedCoordinates = originalCoordinates + optimizationPath * stepSize;

		// Calculate the condition number after the improvement
		current->coordinates = updatedCoordinates;
		double updatedConditionNumber = calculateConditionNumber(index);

		// If condition number is not improved, restore to the previous state
		if (updatedConditionNumber < originalConditionNumber) {
			originalConditionNumber = updatedConditionNumber;
			originalCoordinates = updatedCoordinates;
			stepSize *= 2;
		} else {
			current->coordinates = originalCoordinates;
			stepSize /= 2;
		}

	} while (stepSize > INITIAL_STEP_SIZE);
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
			Triangle triangle(current, next, previous);
			Triangle projected = ProjectionFunctions::projectTriangle(&triangle);

			// Calculate the Jacobian condition number value for the triangle and add it to the total sum
			double conditionNumber = JacobianCalculator::calculateConditionNumber(projected);
			totalValue += conditionNumber;
		}
	}

	return totalValue;
}
