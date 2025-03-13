#ifndef QUEUE_H
#define QUEUE_H

struct Node {
    int cost;
    int bound;
    std::vector<int> path;
    std::vector<bool> visited;
    std::vector<std::vector<int>> matrix;

    Node() {
        this->cost = 0;
        this->bound = 0;
        this->path = {};
        this->visited = {};
        this->matrix = {};
    }

    Node(int _cost, int _bound, std::vector<int> _path, std::vector<bool> _visited, std::vector<std::vector<int>> &_matrix) {
        this->cost = _cost;
        this->bound = _bound;
        this->path = _path;
        this->visited = _visited;
        this->matrix = _matrix;
    }
};

// wierzchołek do trzymania wskaźnika do danej oraz wskaźników
// na poprzedni i następny wierzchołek listy
struct QueueNode {
    Node data;
    QueueNode* next;
    QueueNode* prev;

    QueueNode(Node &node) {
        data = node;
        next = nullptr;
        prev = nullptr;
    }
};


struct Queue {
    QueueNode* _front;
    QueueNode* _back;

    Queue() {
        _front = nullptr;
        _back = nullptr;
    }

    ~Queue() {
        while (!empty()) {
            pop();
        }
    }

    // sprawdza czy kolejka jest pusta
    bool empty() {
        if (_front == nullptr){
            return true;
        }
        return false;
    }

    // zwraca pierwszy element z kolejki
    Node front() {
        if (!empty()) {
            return _front->data;
        }
    }

    // zwraca ostatni element z kolejki
    Node back() {
        if (!empty()) {
            return _back->data;
        }
    }

    // usuwa pierwszy element z kolejki
    void pop() {
        if (empty()) {
            return;
        }
        QueueNode *tmp = _front;
        // jeśli jest tylko jeden wierzchołek w kolejce
        if (_front == _back){
            _front = nullptr;
            _back = nullptr;
        }
        else {
            _front = _front->next;
            _front->prev = nullptr;
        }
        delete tmp;
    }

    // wstawia element na końcu kolejki
    void push(Node &node) {
        QueueNode *new_node = new QueueNode(node);
        if (empty()) {
            _front = new_node;
            _back = new_node;
        }
        else {
            new_node->prev = _back;
            _back->next = new_node;
            _back = new_node;
        }
    }
};


#endif QUEUE_H
