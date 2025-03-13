#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
//#include <xlnt/xlnt.hpp>
#include <random>

#include <math.h>
#include <bits/ranges_algo.h>


struct Parameters {
    std::string file_name;
    std::string output_file_name;
    std::string type;
    int time;
    double alfa;    // temperatura początkowa
    double beta;
    double rho;
    bool show_prints;
    int optimal_cost;
    std::vector<int> optimal_path;
    int ants;
    int iterations_NC;
    int two_opt;
    int error_iter;
};

struct Times {
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::time_point<std::chrono::steady_clock> time_point;
    long long elapsed;
    std::chrono::high_resolution_clock::time_point start_high_res;
    std::chrono::high_resolution_clock::time_point end_high_res;
};

struct NN_data {
    int min_cost;
    std::vector<int> best_path;
};

struct City {
    int id;
    double x, y;
};

struct Ant {
    int current_city;
    double path_cost;
    std::vector<int> tabu_list;
};

class AntColony {
    int node_num;
    int optimal_cost;
    int bound; // potrzebne do nn
    std::vector<int> optimal_path;
    std::vector<std::vector<double>> feromone_intensity;
    std::vector<std::vector<double>> delta_feromone_level;
    std::vector<std::vector<int>> matrix; // dla matrix
    std::vector<City> edges; // dla EUC2D
    std::vector<Ant> ants;
    Parameters parameters;
    Times times;
    NN_data nn_data;
    std::random_device rd;
    std::mt19937 gen;

public:
    AntColony(const std::string &file_name);

    void PrintMatrix(const std::vector<std::vector<double>> & vector);

    void ACO();

    int ChooseNextTown(Ant &ant);

    bool FindIfAllowed(Ant &ant, int j);

    bool CheckIfAllAntsChooseTheSamePath(std::vector<Ant> &ants);

    bool CheckPaths(std::vector<int> first, std::vector<int> second);

    std::vector<int> GeneratePathViaTwoOpt(std::vector<int> path);

    std::vector<double> CreateDecisionTable();

    double CalculateDistanceEUC2D(City &a, City &b);

    int CalculateCostOfPath(std::vector<int> &path);

    int CostOfEdge(int &i, int &j);

    void NN_UpperBound();

    std::vector<int> FollowPath(int current_node, std::vector<bool> visited, std::vector<int> current_path);

    bool NN_CheckIfAllVisited(std::vector<int> &path);

    int NN_PathCost(const std::vector<int> &path);


    void NearestNeighbour(std::vector<City> &cities);


    std::vector<int> NN_GetALLShortestNodes(int current_node, std::vector<bool> &visited);

    int ReadTSPFile(const std::string &file_name);

    // read data input, config and write results to xlsx
    int ReadTSPFileEUC2D(const std::string &file_name);
    // int ReadTSPFileList(const std::string &file_name);
    // void AddEdge(std::vector<std::vector<int>> &adj_list, int u, int v);
    int ReadConfigFile(const std::string &file_name);
    int WriteOutputFile();

    // show result in console
    void PrintResult();
    void PrintPath(const std::vector<int> &path);

    void PrintMatrixEUC2D(std::vector<City> &edges);

    void PrintMatrix(std::vector<std::vector<int>> &matrix);
};

int main(int argc, char *argv[]) {

    srand(time(NULL));
    // zamiast rand() w algorytmie jest std::uniform_real_distribution

    AntColony("config.txt");


    std::cout << "\n\nPress enter to exit" << std::endl;
    std::cin.ignore();
    return 0;
}


AntColony::AntColony(const std::string &file_name) {
    if (ReadConfigFile(file_name) && ReadTSPFile(this->parameters.file_name)) {
        this->optimal_cost = INT_MAX;

        /// do obliczania Tau_0
        this->NN_UpperBound();
        this->bound = this->nn_data.min_cost;

        std::cout << "Koszt z NN: " << this->bound << std::endl;
        std::cout << "\n";

        std::random_device rd;
        std::mt19937 gen(rd());

        this->gen = gen;


        this->times.start = std::chrono::steady_clock::now();
        // this->times.start_high_res = std::chrono::high_resolution_clock::now();


        // std::vector<int> a = {1,2,3,4,1};
        // std::vector<int> b = {3,4,1,2,3};
        // std::vector<int> c = {1,4,2,3,1};
        // std::cout << CheckPaths(a,b);
        // std::cout << CheckPaths(a,b);

        this->ACO();

        // time measurment
        this->times.end = std::chrono::steady_clock::now();
        this->times.elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(this->times.end - this->times.start).count();

        // this->WriteOutputFile();
        this->PrintResult();
    }
}

void AntColony::ACO() {

    double tau_0 = (double) parameters.ants / (double) nn_data.min_cost;
    std::cout << "tau_0: " << tau_0 << "\n\n";

    int iter = 0;

    // 1. Inicjalizacja
    // int t = 0;
    feromone_intensity.resize(node_num, std::vector<double>(node_num, tau_0));
    delta_feromone_level.resize(node_num, std::vector<double>(node_num, 0.00000000000000000000000000000000000001));

    ants.resize(this->parameters.ants);

    int last_update = 0;



    while (iter < parameters.iterations_NC) {

        // startują z różnych wierzchołków
        for (int i = 0; i < ants.size(); i++) {
            ants.at(i).tabu_list.clear();
            ants.at(i).tabu_list.push_back(i % node_num);
        }


        // indeks w liście tabu -- jej aktualny rozmiar
        int s = 0;


        // 2. Powtarzaj dopóki lista tabu nie jest pełna

        // czy tu nie będzie problemu z ostatnim miastem ???
        while (s < node_num-1) {
            s++;
            for (int i = 0; i < ants.size(); i++) {
                int town_j = ChooseNextTown(ants.at(i));
                ants.at(i).tabu_list.push_back(town_j);
            }
        }


        // 3.
        std::vector<int> two_opt_path;
        int cost = 0;
        int two_opt_cost = 0;
        // dla każdej mrówki
        for (int i = 0; i < ants.size(); i++) {
            // =================== może tutaj dodać liczenie trasy ??? =======================

            // PrintPath(ants.at(i).tabu_list);

            cost = CalculateCostOfPath(ants.at(i).tabu_list);
            // std::cout << "cost: " << cost << "\n";

            if (parameters.two_opt == 1) {
                for (int k = 0; k < node_num; k++) {
                    // std::cout << "two_opt:\n";
                    two_opt_path = GeneratePathViaTwoOpt(ants.at(i).tabu_list);
                    // PrintPath(two_opt_path);
                    two_opt_cost = CalculateCostOfPath(two_opt_path);
                    if (two_opt_cost < cost) {
                        cost = two_opt_cost;
                        ants.at(i).tabu_list = two_opt_path;
                    }
                }
            }


            // std::cout << "ant size: " << ants.size() << "\n";
            // sprawdzenie czy znaleziono lepszą ścieżkę ...
            if (cost < this->optimal_cost) {
                last_update = 0;
                this->optimal_cost = cost;
                this->optimal_path = ants.at(i).tabu_list;
            }

            // aktualizuj feromon na drodze każdej mrówki
            if (parameters.type == "sym") {
                for (int j = 0; j < node_num - 1; j++) {
                    // delta_tau += 1 / L_k     -- 1 / długość trasy k-tej mrówki
                    delta_feromone_level.at(ants.at(i).tabu_list.at(j)).at(ants.at(i).tabu_list.at(j + 1)) += 1.0 / (double) cost;
                    delta_feromone_level.at(ants.at(i).tabu_list.at(j + 1)).at(ants.at(i).tabu_list.at(j)) += 1.0 / (double) cost;
                    // feromone_intensity.at(i).at(j) = (1 - parameters.rho) * feromone_intensity.at(i).at(j) + delta_feromone_level.at(i).at(j);

                    // std::cout << "delta: " << delta_feromone_level.at(ants.at(i).tabu_list.at(j)).at(ants.at(i).tabu_list.at(j + 1)) << "\n";
                }
                // std::cout << "delta: " << "\n";
                // PrintMatrix(delta_feromone_level);
            }
            else {
                for (int j = 0; j < node_num - 1; j++) {
                    // delta_tau += 1 / L_k     -- 1 / długość trasy k-tej mrówki
                    delta_feromone_level.at(ants.at(i).tabu_list.at(j)).at(ants.at(i).tabu_list.at(j + 1)) += 1.0 / (double) cost;
                    // feromone_intensity.at(i).at(j) = (1 - parameters.rho) * feromone_intensity.at(i).at(j) + delta_feromone_level.at(i).at(j);

                    // std::cout << "delta: " << delta_feromone_level.at(ants.at(i).tabu_list.at(j)).at(ants.at(i).tabu_list.at(j + 1)) << "\n";
                }
                // std::cout << "delta: " << "\n";
                // PrintMatrix(delta_feromone_level);
            }
        }

        //===================================//
        //  zobaczyć czy dobrze jest
        //  rozkładanie feromonów
        //
        //
        //
        //
        //
        //===================================//

        // 4.   rho  - czy nie powinno być (1 - rho)? -----------------------------------


        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < matrix.size(); j++) {
                // std::cout << "feromone: " << "\n";
                feromone_intensity.at(i).at(j) = (1 - parameters.rho) * feromone_intensity.at(i).at(j) + delta_feromone_level.at(i).at(j);
                // std::cout << feromone_intensity.at(i).at(j) << "\n";
            }
        }


        // set t = t + 1...

        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < matrix.size(); j++) {
                delta_feromone_level.at(i).at(j) = 0.0;
            }
        }

        iter++;
        last_update++;


        // ograniczenie czasowe algorytmu
        this->times.time_point = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::seconds>(
                this->times.time_point - this->times.start).count() >= this->parameters.time) {
            // iter = parameters.iterations_NC;
            break;
                }

        // jeśli znalazło opt
        if (optimal_cost == parameters.optimal_cost || last_update >= parameters.error_iter) {
            break;
        }

        if (parameters.show_prints == 1) {
            std::cout << "iter: " << iter << " cost: " << this->optimal_cost << "\n";
        }


        if (iter % 50 == 0) {
            if(CheckIfAllAntsChooseTheSamePath(ants)) {
                std::cout << "\nsame paths at iteration: " << iter << "\n\n";
                break;
            }
        }
    }
    // std::cout << "feromone matrix:\n";
    // PrintMatrix(feromone_intensity);
    // std::cout << "delta feromone matrix:\n";
    // PrintMatrix(delta_feromone_level);
    std::cout << "\nlast iteration: " << iter << "\n\n";
}


int AntColony::ChooseNextTown(Ant &ant) {
    std::vector<double> probalities_to_choose_town(node_num, 0.0);
    double sum = 0.0;
    // double epsilon = 1e-300;

    // oblicza prawdopobodieństwo wybrania każdego miasta
    for (int k = 0; k < node_num; k++) {
        if (FindIfAllowed(ant, k) == false) {
            // ============================= sprawdzić czy indeksy są dobrze =============================
            probalities_to_choose_town.at(k) = pow(feromone_intensity.at(ant.tabu_list.back()).at(k), parameters.alfa) * pow(
                 1.0 / (double) this->matrix.at(ant.tabu_list.back()).at(k), parameters.beta);
            // probalities_to_choose_town.at(k) = std::max(probalities_to_choose_town.at(k), epsilon);
            sum += probalities_to_choose_town.at(k);
            // std::cout << "feromone intens: " << feromone_intensity.at(ants.at(i+j).tabu_list.size() - 1).at(k) << " ";
            // std::cout << "tour lenght: " << this->matrix.at(ants.at(i+j).tabu_list.back()).at(k) << " ";
        }
        // std::cout << "\n";
    }

    // P_ij = ...
    for (int k = 0; k < probalities_to_choose_town.size(); k++) {
        probalities_to_choose_town.at(k) /= sum;
    }

    int index = -1;
    double max = 0.0;

    // double random = ((double) rand() / (double) (RAND_MAX));
    std::uniform_real_distribution<> dist(0.0, 1.0);
    double random = dist(gen);
    // std::cout << "random: " << random << " sum: " << sum << "\n";

    for (int m = 0; m < probalities_to_choose_town.size(); m++) {
        max += probalities_to_choose_town.at(m);
        if (max >= random) {
            index = m;
            break;
        }
    }

    // jeśli jednak suma prawdopodobieństw nie przekroczy losowej liczby
    // co jest raczej niemożliwe, bo suma prawdopodobieństw jest równa 1
    if (index == -1) {
        for (int m = 0; m < node_num; m++) {
            if (FindIfAllowed(ant, m) == false) {
                index = m;
                break;
            }
        }
    }
    return index;
}

// sprawdza czy miasto jest na liście tabu
bool AntColony::FindIfAllowed(Ant &ant, int j) {
    for (int i = 0; i < ant.tabu_list.size(); i++) {
        if (ant.tabu_list.at(i) == j) {
            return true;
        }
    }
    return false;
}

bool AntColony::CheckIfAllAntsChooseTheSamePath(std::vector<Ant> &ants) {
    std::vector<int> first = ants[0].tabu_list;

    for (int i = 1; i < ants.size(); ++i) {
        if (!CheckPaths(first, ants[i].tabu_list)) {
            return false;
        }
    }

    return true;
}

bool AntColony::CheckPaths(std::vector<int> first, std::vector<int> second) {
    int n = first.size();
    // first.pop_back();
    // second.pop_back();
    // PrintPath(first);
    // PrintPath(second);

    for (int i = 0; i < n; i++) {
        bool same = true;
        for (int j = 0; j < n; j++) {
            // std::cout << "first[i]: " << first[i] << ", sec[i+shift]: " << second[(i + j) % n] << "\n";
            if (first[i] != second[(i + j) % n]) {
                same = false;
                break;
            }
        }
        if (same) return true;
    }
    // std::cout << "end loop\n";

    // first.push_back(first.at(0));
    // second.push_back(second.at(0));

    return false;
}

std::vector<int> AntColony::GeneratePathViaTwoOpt(std::vector<int> path) {
    // inwersja (2-opt, zamiana 2-krawędziowa) - losuje dwa indeksy jako początek i koniec inwersji
    // path.pop_back();
    int i = rand() % path.size();
    int j = rand() % path.size();
    while (i == j) {
        j = rand() % path.size();
    }
    // std::cout << "i = " << i << "j = " << j << std::endl;
    if (i > j) {
        while (j < i) {
            int temp = path[j];
            path[j] = path[i];
            path[i] = temp;
            j++;
            i--;
        }
    }
    else {
        while (i < j) {
            int temp = path[i];
            path[i] = path[j];
            path[j] = temp;
            i++;
            j--;
        }
    }
    // path.push_back(path[0]);
    return path;
}

// std::vector<double> AntColony::CreateDecisionTable() {
//
// }

double AntColony::CalculateDistanceEUC2D(City &a, City &b) {
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

int AntColony::CalculateCostOfPath(std::vector<int> &path) {
    int cost = 0;
    int current_node, next_node;
    for (int i = 0; i < path.size(); i++) {
        current_node = path.at(i);
        next_node = path.at((i + 1) % path.size());
        cost += this->matrix.at(current_node).at(next_node);
    }
    return cost;
}

int AntColony::CostOfEdge(int &i, int &j) {
    return matrix.at(i).at(j);
}










void AntColony::NN_UpperBound() {
    std::vector<int> path;
    std::vector<bool> visited;
    visited.resize(this->node_num, false);
    this->nn_data.min_cost = std::numeric_limits<int>::max();
    int min = std::numeric_limits<int>::max();
    std::vector<int> optimal_path;

    // for (int i = 0; i < this->node_num; i++) {
    //     std::vector<int> result_path = this->FollowPath(i, visited, path);
    //
    //     if (this->nn_data.min_cost < min && this->nn_data.min_cost != 0) {
    //         min = this->nn_data.min_cost;
    //         optimal_path = this->nn_data.best_path;
    //     }
    // }

    std::vector<int> result_path = this->FollowPath(0, visited, path);

    // this->nn_data.best_path = optimal_path;
    // this->nn_data.min_cost = min;
}

std::vector<int> AntColony::FollowPath(int current_node, std::vector<bool> visited,
    std::vector<int> current_path) {

    current_path.push_back(current_node);
    visited.at(current_node) = true;

    if (current_path.size() >= this->node_num) {
        // jeśli doszło do liścia to punkt startowy jest ustawiany jako nieodwiedzony
        // i do niego szukana jest droga
        visited[current_path[0]] = false;
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

    // for (int city : best_paths) {
    if (!best_paths.empty()) {
        std::vector<int> new_path = this->FollowPath(best_paths.at(0), visited, current_path);

        int cost = this->NN_PathCost(new_path);

        if (cost < this->optimal_cost && cost != 0) {
            this->nn_data.min_cost = cost;
            this->nn_data.best_path = new_path;
            // best_path = new_path;
        }
    }
    // }

    return best_path;

}

bool AntColony::NN_CheckIfAllVisited(std::vector<int> &path) {
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

int AntColony::NN_PathCost(const std::vector<int> &path) {
    int cost = 0;
    for (int i = 0; i < path.size(); i++) {
        int current_node = path.at(i);
        int next_node = path.at((i + 1) % path.size());
        cost += this->matrix.at(current_node).at(next_node % path.size());
    }
    return cost;
}

std::vector<int> AntColony::NN_GetALLShortestNodes(int current_node, std::vector<bool> &visited) {
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
 *  algorithm's independent functions
 */

int AntColony::ReadTSPFile(const std::string &file_name) {
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

int AntColony::ReadTSPFileEUC2D(const std::string &file_name) {
    std::ifstream file;
    file.open(file_name);
    if (file.is_open()) {
        std::string line;
        std::getline(file, line);
        std::string cities_n = line.substr(line.find(' ') + 1);
        this->node_num = std::stoi(cities_n);
        this->edges.resize(this->node_num);

        std::getline(file, line);
        std::string type = line.substr(line.find(' ') + 1);
        this->parameters.type = type;

        for (int i = 0; i < this->node_num; i++) {
            std::getline(file, line);
            std::istringstream sline(line);
            for (int j = 0; j < 3; j++) {
                std::getline(sline, line, ' ');
                if (j == 0) {
                    int n = std::stoi(line);
                    this->edges.at(i).id = n;
                }
                else if (j == 1) {
                    double n = std::stod(line);
                    this->edges.at(i).x = n;
                }
                else if (j == 2) {
                    double n = std::stod(line);
                    this->edges.at(i).y = n;
                }
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


int AntColony::ReadConfigFile(const std::string &file_name) {
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

        // time
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.time = std::stoi(substr);

        // iterations_NC
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.iterations_NC = std::stoi(substr);

        // ants
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.ants = std::stoi(substr);

        // alfa
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.alfa = std::stod(substr);

        // beta
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.beta = std::stod(substr);

        // rho
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.rho = std::stod(substr);

        // two_opt
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.rho = std::stoi(substr);

        // show prints
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.show_prints = std::stoi(substr);

        // error iter
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.error_iter = std::stoi(substr);

        file.close();
        return 1;
    }
    std::cout << "Error opening file " << file_name << std::endl;
    return 0;
}



/*
// zmienić wypisywanie w received path, bo problem z + 1
int AntColony::WriteOutputFile() {
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
void AntColony::PrintResult() {
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

    std::cout << "Minimal cost: " << this->optimal_cost << std::endl;
    std::cout << "Received path" << std::endl;
    for (int node : this->optimal_path) {
        std::cout << (node + 1) << " ";
    }
    std::cout << "\n\n";
    // auto time_end_sec = std::chrono::duration_cast<std::chrono::seconds>(this->times.time_point - this->times.start).count();
    std::cout << "Time taken: " << this->times.elapsed << "ms" << std::endl;

    // błąd względny i bezwzględny
    int x = this->parameters.optimal_cost - this->optimal_cost;
    std::cout << "\n\nAbsolute error: " << abs(x) << "\n";
    float y = (float(x)/float(this->parameters.optimal_cost)) * 100;
    std::cout << "Relative error: " << abs(y) << "%";

}


void AntColony::PrintPath(const std::vector<int> &path) {
    for (int node : path) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
};

void AntColony::PrintMatrixEUC2D(std::vector<City> &edges) {
    std::cout << "edges size: " << edges.size() << std::endl;

    std::cout << std::setprecision(1) << std::fixed;
    for (int i = 0; i < edges.size(); i++) {
        std::cout << edges.at(i).id << " " << edges.at(i).x << " " << edges.at(i).y << "\n";
    }
};

void AntColony::PrintMatrix(std::vector<std::vector<int>> &matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            std::cout << matrix.at(i).at(j) << " ";
        }
        std::cout << "\n";
    }
};

void AntColony::PrintMatrix(const std::vector<std::vector<double>> &matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            std::cout << std::setprecision(2) << matrix.at(i).at(j) << " ";
        }
        std::cout << "\n";
    }
};