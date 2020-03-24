#include "DoublyLinkedList.h"

DoublyLinkedList::DoublyLinkedList() {
	head = NULL;
	length = 0;
}

void DoublyLinkedList::insert(NeighborhoodVertex* value) {

	Node* node = new Node();
	node->data = value;

	if (head == NULL) {
		head = node;
		head->next = head;
		head->prev = head;
	} else {
		head->prev->next = node;
		node->prev = head->prev;
		head->prev = node;
		node->next = head;
	}

	length++;
}

Node* DoublyLinkedList::find(NeighborhoodVertex* value) {

	if (length == 0) {
		return NULL;
	}

	Node* temp = head;
	do {

		if (temp->data == value) {
			return temp;
		}

		temp = temp->next;

	} while (temp != head);

	return NULL;
}

unsigned int DoublyLinkedList::getLength() {
	return length;
}

void DoublyLinkedList::clear() {
	
	Node* current = head;
	while (current != current->next) {

		Node* next = current->next;
		Node* prev = current->prev;

		next->prev = prev;
		prev->next = next;

		delete current;
		current = next;
	}

	delete current;
	head = NULL;
	length = 0;
}
