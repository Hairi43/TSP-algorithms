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
    int optimal_cost;
    std::vector<int> optimal_path;
};

struct Times {
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed;
};

class RandomMethod {
    int node_num;
    int min_cost;
    std::vector<int> optimal_path;
    std::vector<int> current_path;
    std::vector<std::vector<int>> matrix;
    Parameters parameters;
    Times times;
public:
    RandomMethod(const std::string &file_name);
    int ReadTSPFile(const std::string &file_name);
    int ReadConfigFile(const std::string &file_name);
    int WriteOutputFile();
    int CalculateCostOfPath();
    void PrintResult();
    void PrintPath(const std::vector<int> &path);
    void PrintMatrix();
};

int main(int argc, char *argv[]) {

    RandomMethod("config.txt");

    std::cout << "\n\nPress enter to exit" << std::endl;
    std::cin.ignore();
    return 0;
}


RandomMethod::RandomMethod(const std::string &file_name) {
    if (ReadConfigFile(file_name) && ReadTSPFile(this->parameters.file_name)) {
        this->min_cost = std::numeric_limits<int>::max();
        for (int i = 0; i < this->node_num; i++) {
            this->current_path.push_back(i);
        }

        int cost = 0;

        auto rd = std::random_device();
        auto rng = std::default_random_engine {rd()};


        this->times.start = std::chrono::high_resolution_clock::now();

        if (this->parameters.type == "sym") {
            for (int i = 0; i < this->parameters.iterations; i++) {
                cost = this->CalculateCostOfPath();
                if (cost < this->min_cost) {
                    this->min_cost = cost;
                    this->optimal_path = this->current_path;
                }
                std::shuffle(this->current_path.begin() + 1, this->current_path.end(), rng);
                // std::random_shuffle(this->current_path.begin() + 1, this->current_path.end());
                // this->PrintPath(this->current_path);
            }
        }
        if (this->parameters.type == "asym") {
            for (int i = 0; i < this->parameters.iterations; i++) {
                cost = this->CalculateCostOfPath();
                if (cost < this->min_cost) {
                    this->min_cost = cost;
                    this->optimal_path = this->current_path;
                }
                std::shuffle(this->current_path.begin(), this->current_path.end(), rng);
            }
        }

        this->times.end = std::chrono::high_resolution_clock::now();
        this->times.elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(
            this->times.end - this->times.start);

        this->WriteOutputFile();
        this->PrintResult();
    }
}

int RandomMethod::ReadTSPFile(const std::string &file_name) {
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

int RandomMethod::ReadConfigFile(const std::string &file_name) {
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

int RandomMethod::WriteOutputFile() {
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
        for (int i : this->optimal_path) {
            file << (i + 1) << " ";
        }
        file << "\n\n";
        file << "Time taken: " << this->times.elapsed.count() << "s";
        file.close();
        return 1;
    }
    std::cout << "Cannot write to output file " << this->parameters.output_file_name << std::endl;
    return 0;
}

int RandomMethod::CalculateCostOfPath() {
    int cost = 0;
    for (int i = 0; i < this->current_path.size(); i++) {
        int current_node = this->current_path.at(i);
        int next_node = this->current_path.at((i + 1) % this->current_path.size());
        cost += this->matrix.at(current_node).at(next_node % this->current_path.size());
    }
    return cost;
}

void RandomMethod::PrintResult() {
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
    for (int node : this->optimal_path) {
        std::cout << (node + 1) << " ";
    }

    std::cout << "\n\nTime taken: " << this->times.elapsed.count() << "s" << std::endl;
}

void RandomMethod::PrintPath(const std::vector<int> &path) {
    for (int node : path) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
};

void RandomMethod::PrintMatrix() {
    std::cout << "nodes: " << this->node_num << std::endl;
    std::cout << "matrix size: " << this->matrix.size() << std::endl;

    for (int i = 0; i < this->matrix.size(); i++) {
        for (int j = 0; j < this->matrix.size(); j++) {
            std::cout << this->matrix.at(i).at(j) << " ";
        }
        std::cout << std::endl;
    }
};
