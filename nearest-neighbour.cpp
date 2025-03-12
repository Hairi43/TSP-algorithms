#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <random>


struct Parameters {
    std::string file_name;
    std::string output_file_name;
    std::string type;
    int iterations;
    int starting_city;
    int optimal_cost;
    std::vector<int> optimal_path;
};

struct Times {
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed;
};

// struct Node {
//     int id;
//     int visited;
// };

class NearestNeighbour {
    int node_num;
    int min_cost;
    std::vector<bool> visited;
    std::vector<int> best_path;
    std::vector<std::vector<int>> matrix;
    Parameters parameters;
    Times times;
public:
    NearestNeighbour(const std::string &file_name);
    int ReadTSPFile(const std::string &file_name);
    int ReadConfigFile(const std::string &file_name);
    int WriteOutputFile();
    // new method
    std::vector<int> FollowPath(int current_node, std::vector<bool> visited, std::vector<int> current_path);
    std::vector<int> GetALLShortestNodes(int current_node, std::vector<bool> &visited);
    int CalculateCostOfPath(const std::vector<int> &path);
    void PrintResult();
    void PrintPath(const std::vector<int> &path);
    void PrintMatrix();
    bool CheckIfAllVisited();
};

int main(int argc, char *argv[]) {

    NearestNeighbour("config.txt");

    std::cout << "\n\nPress enter to exit" << std::endl;
    std::cin.ignore();
    return 0;
}


NearestNeighbour::NearestNeighbour(const std::string &file_name) {
    if (ReadConfigFile(file_name) && ReadTSPFile(this->parameters.file_name)) {

        std::vector<int> path;
        this->visited.resize(this->node_num, false);
        int min_cost = std::numeric_limits<int>::max();
        std::vector<int> optimal_path;

        if (this->parameters.iterations == 0) {
            this->times.start = std::chrono::high_resolution_clock::now();

            std::vector<int> result_path = this->FollowPath(this->parameters.starting_city, this->visited, path);

            this->times.end = std::chrono::high_resolution_clock::now();
            this->times.elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(
                this->times.end - this->times.start);

            this->best_path = result_path;
            this->WriteOutputFile();
            this->PrintResult();
        }
        if (this->parameters.iterations == 1) {
            this->times.start = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < this->node_num; i++) {
                std::vector<int> result_path = this->FollowPath(i, this->visited, path);

                if (this->min_cost < min_cost) {
                    min_cost = this->min_cost;
                    optimal_path = result_path;
                }
            }
            this->times.end = std::chrono::high_resolution_clock::now();
            this->times.elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(
                this->times.end - this->times.start);
            this->best_path = optimal_path;
            this->min_cost = min_cost;
            this->WriteOutputFile();
            this->PrintResult();
        }

    }
}

int NearestNeighbour::ReadTSPFile(const std::string &file_name) {
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

        std::getline(file, line);
        if (line != "EOF") {
            std::string optimal_cost = line.substr(line.find(' ') + 1);
            this->parameters.optimal_cost = std::stoi(optimal_cost);

            std::getline(file, line);
            for (int i = 0; i < this->node_num; i++) {
                std::getline(file, line);
                int n = std::stoi(line);
                this->parameters.optimal_path.push_back(n);
            }
        }
        else this->parameters.optimal_cost = -1;

        file.close();
        return 1;
    }
    std::cout << "Error opening file " << file_name << std::endl;
    return 0;
}

int NearestNeighbour::ReadConfigFile(const std::string &file_name) {
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
        this->parameters.iterations = std::stoi(substr);

        file.close();
        return 1;
    }
    std::cout << "Error opening file " << file_name << std::endl;
    return 0;
}

int NearestNeighbour::WriteOutputFile() {
    std::ofstream file;
    file.open(this->parameters.output_file_name);
    if (file.is_open()) {
        file << "File name: " << this->parameters.file_name << std::endl;
        file << "\nOptimal cost: " << this->parameters.optimal_cost << std::endl;
        file << "Optimal Path" << std::endl;
        for (int node : this->parameters.optimal_path) {
            file << node << " ";
        }
        file << "\n\n";
        file << "Iterations: " << this->parameters.iterations << std::endl;
        file << "Result from program: " << this->min_cost << std::endl;
        file << "Result path" << std::endl;
        for (int node : this->best_path) {
            file << (node + 1) << " ";
        }
        file << "\n\n";
        file << "Time taken: " << this->times.elapsed.count() << "s";
        file.close();
        return 1;
    }
    std::cout << "Cannot write to output file " << this->parameters.output_file_name << std::endl;
    return 0;
}

std::vector<int> NearestNeighbour::FollowPath(int current_node, std::vector<bool> visited,
    std::vector<int> current_path) {

    current_path.push_back(current_node);
    visited.at(current_node) = true;

    if (current_path.size() >= this->node_num) {
        return current_path;
    }

    std::vector<int> best_paths = this->GetALLShortestNodes(current_node, visited);

    // return path if there aren't any paths with the same value
    if (best_paths.empty()) {
        return current_path;
    }

    std::vector<int> best_path;
    int minimal_cost = std::numeric_limits<int>::max();

    for (int city : best_paths) {
        std::vector<int> new_path = this->FollowPath(city, visited, current_path);

        int cost = this->CalculateCostOfPath(new_path);

        if (cost < minimal_cost) {
            minimal_cost = cost;
            best_path = new_path;
        }
    }

    this->min_cost = minimal_cost;
    return best_path;

}

bool NearestNeighbour::CheckIfAllVisited() {
    for (int j = 0; j < this->visited.size(); j++) {
        if (this->visited.at(j) != true) {
            return false;
        }
    }
    return true;
}

int NearestNeighbour::CalculateCostOfPath(const std::vector<int> &path) {
    int cost = 0;
    for (int i = 0; i < path.size(); i++) {
        int current_node = path.at(i);
        int next_node = path.at((i + 1) % path.size());
        cost += this->matrix.at(current_node).at(next_node % path.size());
    }
    return cost;
}

std::vector<int> NearestNeighbour::GetALLShortestNodes(int current_node, std::vector<bool> &visited) {
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

void NearestNeighbour::PrintResult() {
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

    std::cout << "Iterations: " << this->parameters.iterations << std::endl;
    std::cout << "Minimal cost: " << this->min_cost << std::endl;
    std::cout << "Received path" << std::endl;
    for (int node : this->best_path) {
        std::cout << (node + 1) << " ";
    }

    std::cout << "\n\nTime taken: " << this->times.elapsed.count() << "s" << std::endl;
}

void NearestNeighbour::PrintPath(const std::vector<int> &path) {
    for (int node : path) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
}

void NearestNeighbour::PrintMatrix() {
    std::cout << "nodes: " << this->node_num << std::endl;
    std::cout << "matrix size: " << this->matrix.size() << std::endl;

    for (int i = 0; i < this->matrix.size(); i++) {
        for (int j = 0; j < this->matrix.size(); j++) {
            std::cout << this->matrix.at(i).at(j) << " ";
        }
        std::cout << std::endl;
    }
}