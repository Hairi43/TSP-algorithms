//
// Created by Admin on 11.11.2024.
//

#ifndef HEAP_H
#define HEAP_H

#include <vector>

struct Node {
    int cost;
    std::vector<int> path;
    std::vector<bool> visited;

    Node() {
        this->cost = 0;
        this->path = {};
        this->visited = {};
    };

    Node(int _cost, std::vector<int> &_path, std::vector<bool> &_visited) {
        this->cost = _cost;
        this->path = _path;
        this->visited = _visited;
    }
};

// kopiec
struct Heap {
    std::vector<Node> heap;

    Heap() {
        heap = std::vector<Node>();
    }

    void print() {
        for (int i = 0; i < heap.size(); i++){
            std::cout << heap[i].cost << " ";
        }
        std::cout << "\n";
    }

    int parent(int i) {
        return (i - 1) / 2;
    }

    int left (int i) {
        return 2 * i + 1;
    }

    int right (int i) {
        return 2 * i + 2;
    }

    bool empty() {
        return heap.empty();
    }

    int size() {
        return heap.size();
    }

    void swap(int a, int b) {
        Node tmp = heap[a];
        heap[a] = heap[b];
        heap[b] = tmp;
    };

    void HeapifyUp(int index) {
        while (index > 0) {
            int parent_index = parent(index);
            if (heap[index].cost < heap[parent_index].cost) {
                swap(index, parent_index);
                index = parent_index;
            }
            else {
                break;
            }
        }
    }

    void HeapifyDown(int index) {
        while (index < heap.size()) {
            int left_index = left(index);
            int right_index = right(index);
            int smallest = index;

            if (left_index < heap.size() && heap[left_index].cost < heap[smallest].cost) {
                smallest = left_index;
            }
            if (right_index < heap.size() && heap[right_index].cost < heap[smallest].cost) {
                smallest = right_index;
            }
            if (smallest != index) {
                swap(index, smallest);
                index = smallest;
            }
            else {
                break;
            }
        }
    }

    void push(Node &node) {
        heap.push_back(node);
        HeapifyUp(heap.size() - 1);
    }

    void pop() {
        if (heap.empty()) {
            // throw std::out_of_range("heap is empty!");
        }
        Node min = heap[0];
        heap[0] = heap.back();
        heap.pop_back();
        if (!heap.empty()) {
            HeapifyDown(0);
        }
    }

    Node top() {
        return heap[0];
    }
};


// kolejka priorytetowa
struct Priority_queue {
    Heap heap;

    void print() {
        heap.print();
    }

    bool empty() {
        return heap.empty();
    }

    int size() {
        return heap.size();
    }

    Node top() {
        return heap.top();
    }

    void push(Node &node) {
        heap.push(node);
    }

    void pop() {
        heap.pop();
    }
};

#endif //HEAP_H
