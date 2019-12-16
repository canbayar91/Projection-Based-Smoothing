#ifndef PRIMITIVE
#define PRIMITIVE

#include "DoublyLinkedList.h"
#include "LineSegment.h"
#include <string>

class Primitive {
public:

	virtual void output(std::string filename) const = 0;

	Node* findNode(NeighborhoodVertex* vertex) {
		return searcher.find(vertex);
	}

protected:

	DoublyLinkedList searcher;

};

#endif
