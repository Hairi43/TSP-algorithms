#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <queue>    // to jest tylko do porównania mojej implementacji z tą z biblioteki std
#include <xlnt/xlnt.hpp>

#include "heap.h"

#define MAX_EDGE 9999

struct Parameters {
    std::string file_name;
    std::string output_file_name;
    std::string type;
    int time;
    int starting_city;
    int optimal_cost;
    std::vector<int> optimal_path;
};

struct Times {
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::time_point<std::chrono::steady_clock> time_point;
    // std::chrono::duration<double> elapsed;
    long long elapsed;
    std::chrono::high_resolution_clock::time_point start_high_res;
    std::chrono::high_resolution_clock::time_point end_high_res;
};

struct NN_data {
    int min_cost;
    std::vector<int> best_path;
};

class BestFirstSearch {
    int node_num;
    int min_cost;
    int bound;
    std::vector<int> optimal_path;
    std::vector<std::vector<int>> matrix;
    Parameters parameters;
    Times times;
    NN_data nn_data;

public:
    BestFirstSearch(const std::string &file_name);

    void BFSAlgorithm();

    void CalculateNodeBound(Node &node, Node parent_node, int city);

    int ReduceMatrix(std::vector<std::vector<int>> &matrix);

    // usunąć
    bool IsVisited(std::vector<int> &path, int i);

    int NodeBound(Node &node);

    int ShortestEdgeIndex(int node, std::vector<bool> &visited);

    int GetEdgeWeight(int i, int j);

    void Upperbound();

    // funkcje do obliczenia górnej granicy za pomocą NN
    void NN_UpperBound();
    std::vector<int> FollowPath(int current_node, std::vector<bool> visited, std::vector<int> current_path);
    bool NN_CheckIfAllVisited(std::vector<int> &path);
    int NN_PathCost(const std::vector<int> &path);
    std::vector<int> NN_GetALLShortestNodes(int current_node, std::vector<bool> &visited);

    // read data input, config and write results to xlsx
    int ReadTSPFile(const std::string &file_name);
    int ReadConfigFile(const std::string &file_name);
    int WriteOutputFile();

    // show result in console
    void PrintResult();
    void PrintPath(const std::vector<int> &path);
    void PrintMatrix(const std::vector<std::vector<int>> & vector);
};

int main(int argc, char *argv[]) {

    BestFirstSearch("config.txt");

    std::cout << "\n\nPress enter to exit" << std::endl;
    std::cin.ignore();
    return 0;
}


BestFirstSearch::BestFirstSearch(const std::string &file_name) {
    if (ReadConfigFile(file_name) && ReadTSPFile(this->parameters.file_name)) {
        this->min_cost = std::numeric_limits<int>::max();

        this->times.start = std::chrono::steady_clock::now();
        // this->times.start_high_res = std::chrono::high_resolution_clock::now();

        // obliczenie upper bound
        int tmp = this->parameters.starting_city;
        this->NN_UpperBound();
        this->parameters.starting_city = tmp;
        this->bound = this->nn_data.min_cost;

        std::cout << "Upper bound: " << this->bound << std::endl;
        this->times.time_point = std::chrono::steady_clock::now();
        std::cout << "Time taken for calculating UB: " << std::chrono::duration_cast<std::chrono::duration<double>>(
            this->times.time_point - this->times.start).count() << "s" << std::endl;


        // odkomentować jeśli chcę zmierzyć czas osobno dla liczenia UB i algorytmu lowest cost
        // this->times.start = std::chrono::steady_clock::now();

        this->BFSAlgorithm();

        this->times.end = std::chrono::steady_clock::now();
        this->times.elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(this->times.end - this->times.start).count();


        // time measurment
        // this->times.end_high_res = std::chrono::high_resolution_clock::now();
        // this->times.elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(
        //     this->times.end_high_res - this->times.start_high_res);

        this->WriteOutputFile();
        this->PrintResult();
    }
}

void BestFirstSearch::BFSAlgorithm() {
    Priority_queue pq;  // implementacja w heap.h

    /*
     *  to jest tylko do sprawdzenia jak czasowo różni się to co napisałem od poprawnej implementacji
     *  wychodzi, że róznica to ok. 50% szybciej na korzyść biblioteki std
     */
    // std::priority_queue<Node, std::vector<Node>, Compare> pq;

    std::vector<int> path = {this->parameters.starting_city};
    std::vector<bool> visited(this->node_num, false);
    std::vector<std::vector<int>> _matrix = this->matrix;

    int reduced_matrix_cost = this->ReduceMatrix(_matrix);
    Node node(reduced_matrix_cost, 0, path, visited, _matrix);
    // Node tmp(0, 0, path, visited, _matrix);
    node.visited[this->parameters.starting_city] = true;
    // this->CalculateNodeBound(node, tmp, this->parameters.starting_city);

    pq.push(node);

    // int min_cost = INT_MAX;

    // this->Upperbound();
    // this->bound = this->nn_data.min_cost;
    // this->bound = 250;
    std::cout << "Upper bound: " << this->bound << std::endl;

    while (!pq.empty()) {
        // sprawdzenie czy przekrodzono czas programu z config'a
        this->times.time_point = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::minutes>(
                this->times.time_point - this->times.start).count() >= this->parameters.time) {break;}

        node = pq.top();
        pq.pop();

        // jeśli górna granica się zmniejszy to niektóre wierzchołki z kolejki nie pójdą dalej
        if (node.cost > this->bound) {
            continue;
        }

        // jeśli odwiedzone zostały wszystkie wierzchołki to sprawdzany jest powrót do miasta
        if (node.path.size() == this->node_num && this->matrix[node.path.back()][this->parameters.starting_city] != MAX_EDGE) {
            // int cost = node.cost + this->matrix[node.path.back()][this->parameters.starting_city];
            int cost = node.cost;
            node.path.push_back(this->parameters.starting_city);
            if (cost <= this->bound) {
                min_cost = cost;
                this->bound = cost;
                this->optimal_path = node.path;
                this->min_cost = min_cost;
            }
        }
        else {
            for (int i = 0; i < this->node_num; i++) {
                if(!node.visited[i] && node.matrix[node.path.back()][i] != MAX_EDGE) {
                    Node child(node.cost, 0, node.path, node.visited, node.matrix);
                    child.path.push_back(i);
                    this->CalculateNodeBound(child, node, i);
                    if (child.cost <= this->bound) {
                        child.visited[i] = true;
                        pq.push(child);
                    }
                    // std::cout << "--------------------------------------------\n";
                }
            }
            // if ()
        }
    }
    std::cout << "not found\n";
}



void BestFirstSearch::CalculateNodeBound(Node &node, Node parent_node, int city) {
    node.matrix = parent_node.matrix;

    int second_to_last = node.path[node.path.size() - 2];
    // std::cout << "second to last: " << second_to_last << "\n";
    for (int i = 0; i < node.matrix.size(); i++) {
        // blokujemy wiersz i kolumnę (ścieżki z ostatniego miasta i do nowego miasta)
        node.matrix[second_to_last][i] = MAX_EDGE;
        node.matrix[i][city] = MAX_EDGE;
    }
    node.matrix[city][this->parameters.starting_city] = MAX_EDGE;

    // this->PrintPath(node.path);
    // std::cout << "before reduction: \n";
    // this->PrintMatrix(node.matrix);

    // koszt aktualnego wierzchołka
    node.cost = parent_node.cost + parent_node.matrix[second_to_last][city];

    if (node.cost > MAX_EDGE) {
        return;
    }

    // obliczenie redukcji macierzy
    node.cost += this->ReduceMatrix(node.matrix);

    // std::cout << "reduced matrix: \n";
    // this->PrintMatrix(node.matrix);
    //
    // std::cout << "parent cost: " << parent_node.cost <<
    // "\nparent matrix cost: " << parent_node.matrix[second_to_last][city] <<
    // std::cout << "\ntotal cost: " << node.cost << std::endl;
}


// liczy koszt redukcji macierzy i ją redukuje
int BestFirstSearch::ReduceMatrix(std::vector<std::vector<int>> &matrix) {
    int cost = 0;

    // redukcja wierszy
    for (int i = 0; i < matrix.size(); i++) {
        int row_min = MAX_EDGE;
        for (int j = 0; j < matrix.size(); j++) {
            if (matrix[i][j] < row_min) {
                row_min = matrix[i][j];
            }
        }
        if (row_min != MAX_EDGE && row_min > 0) {
            cost += row_min;
            for (int j = 0; j < matrix.size(); j++) {
                if (matrix[i][j] != MAX_EDGE) {
                    matrix[i][j] -= row_min;
                }
            }
        }
    }

    // redukcja kolumn
    for (int i = 0; i < matrix.size(); i++) {
        int col_min = MAX_EDGE;
        for (int j = 0; j < matrix.size(); j++) {
            if (matrix[j][i] < col_min) {
                col_min = matrix[j][i];
            }
        }
        if (col_min != MAX_EDGE && col_min > 0) {
            cost += col_min;
            for (int j = 0; j < matrix.size(); j++) {
                if (matrix[j][i] != MAX_EDGE) {
                    matrix[j][i] -= col_min;
                }
            }
        }
    }
    // std::cout << "reduction sum: " << cost << std::endl;
    return cost;

}


/*
 *  algorytm NN do określenia Upper bound
 */

 void BestFirstSearch::Upperbound() {
    std::vector<bool> visited (this->node_num, false);
    std::vector<int> path;
    path.push_back(this->parameters.starting_city);
    visited[this->parameters.starting_city] = true;
    int current_node = 0;
    int cost = 0;

    for (int i = 0; i < this->node_num - 1; i++) {
        int distance = INT_MAX;
        int nearest_node = -1;

        for (int j = 0; j < this->node_num; j++) {
            if (!visited[j] && this->matrix[current_node][j] != 0 && this->matrix[current_node][j] < distance) {
                distance = this->matrix[current_node][j];
                nearest_node = j;
            }
        }

        cost += distance;
        visited[nearest_node] = true;
        current_node = nearest_node;
        path.push_back(current_node);
    }
    if (this->matrix[current_node][path[0]]) {
        cost += this->matrix[current_node][path[0]];
        path.push_back(path[0]);
        this->nn_data.min_cost = cost;
        this->nn_data.best_path = path;
    }
}

void BestFirstSearch::NN_UpperBound() {
    std::vector<int> path;
    std::vector<bool> visited;
    visited.resize(this->node_num, false);
    this->nn_data.min_cost = std::numeric_limits<int>::max();
    int min = std::numeric_limits<int>::max();
    std::vector<int> optimal_path;

    for (int i = 0; i < this->node_num; i++) {
        this->parameters.starting_city = i;
        std::vector<int> result_path = this->FollowPath(i, visited, path);

        if (this->nn_data.min_cost < min && this->nn_data.min_cost != 0) {
            min = this->nn_data.min_cost;
            optimal_path = this->nn_data.best_path;
        }
    }

    this->nn_data.best_path = optimal_path;
    this->nn_data.min_cost = min;
}

std::vector<int> BestFirstSearch::FollowPath(int current_node, std::vector<bool> visited,
    std::vector<int> current_path) {

    current_path.push_back(current_node);
    visited.at(current_node) = true;

    if (current_path.size() >= this->node_num) {
        visited[this->parameters.starting_city] = false;
        std::vector<int> best_paths = this->NN_GetALLShortestNodes(current_node, visited);
        // jest tylko jeden sposób na powrót do wierzchołka startowego
        // jeśli jest droga to pójdzie, a jeśli nie to zwróci pusty wektor,
        // który oznacza, że nie ma drogi (dla grafów asymetrycznych)
        if (!best_paths.empty()) {
            current_path.push_back(best_paths[0]);
            return current_path;
        }
        return current_path = {};
        // return {};
    }

    std::vector<int> best_paths = this->NN_GetALLShortestNodes(current_node, visited);

    std::vector<int> best_path;

    for (int city : best_paths) {
        std::vector<int> new_path = this->FollowPath(city, visited, current_path);

        int cost = this->NN_PathCost(new_path);

        if (cost < this->min_cost && cost != 0) {
            this->nn_data.min_cost = cost;
            this->nn_data.best_path = new_path;
            // best_path = new_path;
        }
    }

    return best_path;

}

bool BestFirstSearch::NN_CheckIfAllVisited(std::vector<int> &path) {
    bool n;
    for (int i = 0; i < path.size(); i++) {
        n = false;
        for (int j = 0; j < path.size(); j++) {
            if (path[j] == i) {
                n = true;
                break;
            }
        }
        if (!n) {
            return false;
        }
    }
    return true;
}

int BestFirstSearch::NN_PathCost(const std::vector<int> &path) {
    int cost = 0;
    for (int i = 0; i < path.size(); i++) {
        int current_node = path.at(i);
        int next_node = path.at((i + 1) % path.size());
        cost += this->matrix.at(current_node).at(next_node % path.size());
    }
    return cost;
}

std::vector<int> BestFirstSearch::NN_GetALLShortestNodes(int current_node, std::vector<bool> &visited) {
    int min = std::numeric_limits<int>::max();
    std::vector<int> nearest_nodes;

    for (int i = 0; i < this->node_num; i++) {
        if (this->matrix.at(current_node).at(i) < min && this->matrix.at(current_node).at(i) > 0 && visited.at(i) == false) {
            min = this->matrix.at(current_node).at(i);
        }
    }

    for (int i = 0; i < this->node_num; i++) {
        if (this->matrix.at(current_node).at(i) == min && visited.at(i) == false) {
            nearest_nodes.push_back(i);
        }
    }
    return nearest_nodes;
}



/*
 *  algorithms's indepentadant functions
 */

int BestFirstSearch::ReadTSPFile(const std::string &file_name) {
    std::ifstream file;
    file.open(file_name);
    if (file.is_open()) {
        std::string line;
        std::getline(file, line);
        std::string cities_n = line.substr(line.find(' ') + 1);
        this->node_num = std::stoi(cities_n);
        this->matrix.resize(this->node_num, std::vector<int>(this->node_num));

        std::getline(file, line);
        std::string type = line.substr(line.find(' ') + 1);
        this->parameters.type = type;

        for (int i = 0; i < this->node_num; i++) {
            std::getline(file, line);
            std::istringstream sline(line);
            for (int j = 0; j < this->node_num; j++) {
                std::getline(sline, line, ' ');
                int n = std::stoi(line);
                this->matrix.at(i).at(j) = n;
            }
        }

        this->parameters.optimal_cost = -1;
        this->parameters.optimal_path = std::vector<int> {};

        std::getline(file, line);
        if (line != "EOF") {
            std::string optimal_cost = line.substr(line.find(' ') + 1);
            this->parameters.optimal_cost = std::stoi(optimal_cost);

            std::getline(file, line);
            if (line != "EOF") {
                for (int i = 0; i < this->node_num; i++) {
                    std::getline(file, line);
                    int n = std::stoi(line);
                    this->parameters.optimal_path.push_back(n);
                }
            }
        }
        else this->parameters.optimal_cost = -1;

        file.close();
        return 1;
    }
    std::cout << "Error opening file " << file_name << std::endl;
    return 0;
}

int BestFirstSearch::ReadConfigFile(const std::string &file_name) {
    std::ifstream file;
    file.open(file_name);
    if (file.is_open()) {
        std::string line;

        // input file
        std::getline(file, line);
        std::string substr = line.substr(line.find(' ') + 1);
        this->parameters.file_name = substr;

        // output file
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.output_file_name = substr;

        /*
         *  wartość z TSPLIB ?
         */

        // starting city
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.starting_city = std::stoi(substr);
        this->parameters.starting_city -= 1;

        // iterations
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.time = std::stoi(substr);

        file.close();
        return 1;
    }
    std::cout << "Error opening file " << file_name << std::endl;
    return 0;
}




// zmienić wypisywanie w received path, bo problem z + 1
int BestFirstSearch::WriteOutputFile() {
    xlnt::workbook wb;
    xlnt::worksheet ws = wb.active_sheet();

    ws.cell("A1").value("File name:");
    ws.cell("B1").value(this->parameters.file_name);
    ws.cell("A2").value("Optimal cost:");
    ws.cell("B2").value(this->parameters.optimal_cost);
    ws.cell("A3").value("Optimal path");
    if (!this->parameters.optimal_path.empty()) {
        for (int i = 4; i < this->node_num + 4; i++) {
            // przeusnięte o 4 indeksy, żeby zaczęło od prawidłowej komórki w pliku xls
            ws.cell("A" + std::to_string(i)).value(this->parameters.optimal_path[i-4]);
        }
    }
    ws.cell("C2").value("Received cost:");
    ws.cell("D2").value(this->min_cost);
    ws.cell("C3").value("Received path");
    if (!this->optimal_path.empty()) {
        for (int i = 4; i < this->node_num + 4; i++) {
            // przeusnięte o 4 indeksy, żeby zaczęło od prawidłowej komórki w pliku xls
            // + 1, bo sama trasa jest obliczana od węzła 0 do n-1

            // podmienić
            ws.cell("C" + std::to_string(i)).value(this->optimal_path[i-4] + 1);
        }
    }
    ws.cell("F2").value("Time taken [ms]");
    // auto time_end_sec = std::chrono::duration_cast<std::chrono::seconds>(this->times.time_point - this->times.start).count();
    ws.cell("F3").value(this->times.elapsed);
    // ws.cell("G2").value("Time taken [min]");
    // auto time_end_min = std::chrono::duration_cast<std::chrono::minutes>(this->times.time_point - this->times.start).count();
    // ws.cell("G3").value(time_end_min);

    wb.save("output.xlsx");

    return 0;
}
// */

// tu też zmienić optimal path
void BestFirstSearch::PrintResult() {
    std::cout << "File name: " << this->parameters.file_name << std::endl;
    std::cout << std::endl;

    if (this->parameters.optimal_cost >= 0) {
        std::cout << "Optimal cost: " << this->parameters.optimal_cost << std::endl;
        std::cout << "Optimal path" << std::endl;
        for (int node : this->parameters.optimal_path) {
            std::cout << node << " ";
        }
        std::cout << "\n\n";
    }

    std::cout << "Minimal cost: " << this->min_cost << std::endl;
    std::cout << "Received path" << std::endl;
    for (int node : this->optimal_path) {
        std::cout << (node + 1) << " ";
    }
    std::cout << "\n\n";
    // auto time_end_sec = std::chrono::duration_cast<std::chrono::seconds>(this->times.time_point - this->times.start).count();
    std::cout << "Time taken: " << this->times.elapsed << "ms" << std::endl;

}


void BestFirstSearch::PrintPath(const std::vector<int> &path) {
    for (int node : path) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
};

void BestFirstSearch::PrintMatrix(const std::vector<std::vector<int>> &matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            std::cout << matrix.at(i).at(j) << " ";
        }
        std::cout << std::endl;
    }
};
