#ifndef STACK_LIBRARY_H
#define STACK_LIBRARY_H

#include <vector>

struct Node {
    int id;
    int cost;
    int bound;
    std::vector<int> path;
    std::vector<bool> visited;
    // std::vector<std::vector<int>> matrix;

    Node() {
        this->id = 0;
        this->cost = 0;
        this->bound = 0;
        this->path = {};
        this->visited = {};
        // this->matrix = {};
    }

    Node(int _id, int _cost, int _bound, std::vector<int> _path, std::vector<bool> _visited) {
        this->id = _id;
        this->cost = _cost;
        this->bound = _bound;
        this->path = _path;
        this->visited = _visited;
        // this->matrix = _matrix;
    }
};

// potrzebne do podmiany tablic w strukturze stosu
template<typename T>
void swap_values(T& a, T& b) {
    T tmp(a);
    a = b;
    b = tmp;
}

template<typename T>
struct Stack {
    int stack_top;
    int capacity;
    T* buffer;

    Stack() {
        stack_top = 0;
        capacity = 16;
        buffer = new T[capacity];
    }

    int size() {
        return stack_top;
    }

    bool empty() {
        if (stack_top == 0) {
            return true;
        }
        return false;
    }

    // zwróć wartość ze szczytu stosu
    T top() {
        if (!empty()) {
            return buffer[stack_top - 1];
        }
        throw std::out_of_range("stack is empty");
    }

    // dodanie elementu do stosu
    void push(T &n) {
        // jeśli nie można dodać do tablicy to rozszerzy ją i kopiuje dane
        if (stack_top == capacity) {
            capacity *= 2;
            T* new_buffer = new T[capacity];
            for (int i = 0; i < stack_top; i++) {
                new_buffer[i] = buffer[i];
            }
            // podmiana wartości w pamięci, żeby w strukturze stosu dalej był buffer
            swap_values(buffer, new_buffer);
            delete []new_buffer;
        }
        // dodanie elementu
        buffer[stack_top] = n;
        stack_top++;
    }

    // "usuwa ze stosu", zmniejsza szczyt stosu, a przy dodawaniu w push()
    // element "nieużywany" jest nadpisywany
    T pop() {
        if (!empty()) {
            stack_top--;
            return buffer[stack_top];
        }
        throw std::out_of_range("stack is empty");
    }


};

#endif //STACK_LIBRARY_H
