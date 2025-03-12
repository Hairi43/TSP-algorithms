#include <algorithm> // next_permutation
#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <xlnt/xlnt.hpp>


struct Parameters {
    std::string file_name;
    std::string output_file_name;
    std::string type;
    int time;
    int optimal_cost;
    std::vector<int> optimal_path;
};

struct Times {
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::time_point<std::chrono::steady_clock> time_point;
    std::chrono::duration<double> elapsed;
    std::chrono::high_resolution_clock::time_point start_high_res;
    std::chrono::high_resolution_clock::time_point end_high_res;
};

class BruteForce {
    int node_num;
    int min_cost;
    std::vector<int> optimal_path;
    std::vector<int> current_path;
    std::vector<std::vector<int>> matrix;
    Parameters parameters;
    Times times;
public:
    BruteForce(const std::string &file_name);
    int ReadTSPFile(const std::string &file_name);
    int ReadConfigFile(const std::string &file_name);
    int WriteOutputFile();
    int CalculateCostOfPath();
    void PrintResult();
    void PrintPath(const std::vector<int> &path);
    void PrintMatrix();
};

int main(int argc, char *argv[]) {

    BruteForce("config.txt");

    std::cout << "\n\nPress enter to exit" << std::endl;
    std::cin.ignore();
    return 0;
}


BruteForce::BruteForce(const std::string &file_name) {
    if (ReadConfigFile(file_name) && ReadTSPFile(this->parameters.file_name)) {
        this->min_cost = std::numeric_limits<int>::max();
        for (int i = 0; i < this->node_num; i++) {
            this->current_path.push_back(i);
        }

        this->times.start = std::chrono::steady_clock::now();
        this->times.start_high_res = std::chrono::high_resolution_clock::now();
        auto duration = this->times.start.time_since_epoch();


        if (this->parameters.type == "sym") {
            int cost = 0;
            do {
                // this->PrintPath(this->current_path);

                cost = this->CalculateCostOfPath();
                if (cost != -1 && cost < this->min_cost) {
                    this->min_cost = cost;
                    this->optimal_path = this->current_path;
                }


                this->times.time_point = std::chrono::steady_clock::now();
                if (std::chrono::duration_cast<std::chrono::minutes>(
                    this->times.time_point - this->times.start).count() >= this->parameters.time) {
                    break;
                }

            } while (std::next_permutation(this->current_path.begin() + 1, this->current_path.end()));
        }

        if (this->parameters.type == "asym") {
            int cost = 0;
            do {
                // this->PrintPath(this->current_path);

                cost = this->CalculateCostOfPath();
                if (cost != -1 && cost < this->min_cost) {
                    this->min_cost = cost;
                    this->optimal_path = this->current_path;
                }

                this->times.time_point = std::chrono::steady_clock::now();
                if (std::chrono::duration_cast<std::chrono::minutes>(
                    this->times.time_point - this->times.start).count() >= this->parameters.time) {
                    break;
                }

                // checks every possible permutation unlike in symmetrical one
            } while (std::next_permutation(this->current_path.begin(), this->current_path.end()));
        }


        // time measurment
        this->times.end_high_res = std::chrono::high_resolution_clock::now();
        this->times.elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(
            this->times.end_high_res - this->times.start_high_res);

        // this->WriteOutputFile();
        this->PrintResult();
    }
}

int BruteForce::ReadTSPFile(const std::string &file_name) {
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

int BruteForce::ReadConfigFile(const std::string &file_name) {
    std::ifstream file;
    file.open(file_name);
    if (file.is_open()) {
        std::string line;
        std::getline(file, line);
        std::string substr = line.substr(line.find(' ') + 1);
        this->parameters.file_name = substr;

        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.output_file_name = substr;


        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.time = std::stoi(substr);

        file.close();
        return 1;
    }
    std::cout << "Error opening file " << file_name << std::endl;
    return 0;
}

int BruteForce::WriteOutputFile() {
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
            ws.cell("C" + std::to_string(i)).value(this->optimal_path[i-4] + 1);
        }
    }
    ws.cell("F2").value("Time taken [s]");
    // auto time_end_sec = std::chrono::duration_cast<std::chrono::seconds>(this->times.time_point - this->times.start).count();
    ws.cell("F3").value(this->times.elapsed.count());
    ws.cell("G2").value("Time taken [min]");
    auto time_end_min = std::chrono::duration_cast<std::chrono::minutes>(this->times.time_point - this->times.start).count();
    ws.cell("G3").value(time_end_min);

    wb.save("output.xlsx");

    return 0;
}
// */


int BruteForce::CalculateCostOfPath() {
    int cost = 0;
    int current_node, next_node;
    for (int i = 0; i < this->current_path.size(); i++) {
        current_node = this->current_path.at(i);
        next_node = this->current_path.at((i + 1) % this->current_path.size());
        // std::cout << next_node << ", ";

        // check if path is valid. 0 means there is no connection between nodes so path is invalid
        if (this->matrix.at(current_node).at(next_node) == 0 && current_node != next_node) {
            return -1;
        }
        cost += this->matrix.at(current_node).at(next_node);
        // std::cout << "cost = " << cost << "\n";
    }
    return cost;
}

void BruteForce::PrintResult() {
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
    std::cout << "Time taken: " << this->times.elapsed.count() << "s" << std::endl;

    int x = this->parameters.optimal_cost - this->min_cost;
    std::cout << "\n\nAbsolute error: " << abs(x) << "\n";
    float y = (float(x)/float(this->parameters.optimal_cost)) * 100;
    std::cout << "Relative error: " << abs(y) << "%";
}


void BruteForce::PrintPath(const std::vector<int> &path) {
    for (int node : path) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
};

void BruteForce::PrintMatrix() {
    std::cout << "nodes: " << this->node_num << std::endl;
    std::cout << "matrix size: " << this->matrix.size() << std::endl;

    for (int i = 0; i < this->matrix.size(); i++) {
        for (int j = 0; j < this->matrix.size(); j++) {
            std::cout << this->matrix.at(i).at(j) << " ";
        }
        std::cout << std::endl;
    }
};
