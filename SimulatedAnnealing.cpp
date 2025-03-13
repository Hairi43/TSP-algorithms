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

#define MAX_EDGE 9999
#define E 2.71828
// #define MIN_TEMP 0.000000000001

struct Parameters {
    std::string file_name;
    std::string output_file_name;
    std::string type;
    int time;
    double temp;    // temperatura początkowa
    double multiplier; // mnożnik jak szybko ma się schładzać temperatura
    double min_temp;
    int show_prints;
    int optimal_cost;
    std::vector<int> optimal_path;
    int start_with_nn;
    int cooling_scheme;
    int iterations;
    int iterations_without_improvment;
    int local_search;
    int end;
    int pocz;
    long double acceptance;
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

class SimulatedAnnealing {
    int node_num;
    int optimal_cost;
    int bound; // potrzebne do nn
    std::vector<int> optimal_path;
    std::vector<std::vector<int>> matrix;
    Parameters parameters;
    Times times;
    NN_data nn_data;

public:
    SimulatedAnnealing(const std::string &file_name);

    void SA();

    std::vector<int> GeneratePathViaTwoOpt(std::vector<int> path);

    std::vector<int> GeneratePathViaSwap(std::vector<int> path);

    std::vector<int> GeneratePathViaInsertion(std::vector<int> path);

    std::vector<int> Shuffle(std::vector<int> path);

    std::vector<int> ShuffleNTimes(std::vector<int> path);

    std::vector<int> SwapRandom(std::vector<int> path);

    int CalculateCostOfPath(std::vector<int> &path);


    // funkcje do obliczenia górnej granicy za pomocą NN
    void Upperbound();
    void NN_UpperBound();
    std::vector<int> FollowPath(int current_node, std::vector<bool> visited, std::vector<int> current_path);
    bool NN_CheckIfAllVisited(std::vector<int> &path);
    int NN_PathCost(const std::vector<int> &path);
    std::vector<int> NN_GetALLShortestNodes(int current_node, std::vector<bool> &visited);

    // read data input, config and write results to xlsx
    int ReadTSPFileArray(const std::string &file_name);
    // int ReadTSPFileList(const std::string &file_name);
    // void AddEdge(std::vector<std::vector<int>> &adj_list, int u, int v);
    int ReadConfigFile(const std::string &file_name);
    int WriteOutputFile();

    // show result in console
    void PrintResult();
    void PrintPath(const std::vector<int> &path);
    void PrintMatrix(std::vector<std::vector<int>> matrix);
};

int main(int argc, char *argv[]) {


    srand(time(NULL));
    // zamiast rand() w algorytmie jest std::uniform_real_distribution

    SimulatedAnnealing("config.txt");


    std::cout << "\n\nPress enter to exit" << std::endl;
    std::cin.ignore();
    return 0;
}


SimulatedAnnealing::SimulatedAnnealing(const std::string &file_name) {
    if (ReadConfigFile(file_name) && ReadTSPFileArray(this->parameters.file_name)) {
        this->optimal_cost = INT_MAX;

        this->times.start = std::chrono::steady_clock::now();
        // this->times.start_high_res = std::chrono::high_resolution_clock::now();

        // obliczenie startowego rozwiązania za pomocą NN, które później będzie przeszukiwane tabu search
        this->NN_UpperBound();
        this->bound = this->nn_data.min_cost;

        std::cout << "koszt trasy z NN: " << this->bound << std::endl;
        std::cout << "\n";

        // std::vector<int> tmp;
        // for (int i = 0; i < 10; i++) {
        // tmp = this->Shuffle(this->nn_data.best_path);
        //     // this->PrintPath(tmp);
        // }



        // odkomentować jeśli chcę zmierzyć czas osobno dla liczenia UB i algorytmu SimulatedAnnealing
        this->times.start = std::chrono::steady_clock::now();

        /*
         *  główna funkcja
         */


        // std::vector<int> v = {1,2,3,4,5,6,7,8,9};
        // v = this->GeneratePathViaInsertion(v);
        // this->PrintPath(v);



        this->SA();

        // time measurment
        this->times.end = std::chrono::steady_clock::now();
        this->times.elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(this->times.end - this->times.start).count();

        // this->WriteOutputFile();
        this->PrintResult();
    }
}

void SimulatedAnnealing::SA() {
    std::vector<int> current_path;
    std::vector<int> new_solution;
    std::vector<int> best_path;
    int best_cost;

    // nearest neighbour jako początkowy
    if (this->parameters.start_with_nn == 1) {
        current_path = this->nn_data.best_path;
        best_path = this->nn_data.best_path;
        best_cost = this->nn_data.min_cost;
    }
    // losowy jako początkowy
    else {
        for (int i = 0; i < this->node_num; i++) {
            current_path.push_back(i);
        }
        current_path.push_back(current_path[0]);
        current_path = this->Shuffle(current_path);
        best_path = current_path;
        best_cost = this->CalculateCostOfPath(current_path);
        std::cout << "Prosze nie zwracac uwagi na to co jest na gorze. Rozwiazanie poczatkowe jest losowo generowane" << std::endl;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);



    /*
     *
     *  temperaturę początkową liczyć jako T_0 = delta_+ / ln(x_0)
     *  a logarytmiczną to T_0 = T_0 * (wynik z NN / losową trasę)
     *
     */
    std::vector<int> path = {};
    for (int i = 0; i < this->node_num; i++) {
        path.push_back(i);
    }
    path.push_back(path[0]);
    int cost = CalculateCostOfPath(path);
    int cost2 = CalculateCostOfPath(path);
    double delta_t = 0;
    for (int i = 0; i < node_num * 3; i++) {
        path = Shuffle(path);
        cost2 = cost;
        path = GeneratePathViaSwap(path);
        cost = CalculateCostOfPath(path);
        delta_t += cost - cost2;
        // std::cout << "cost: " << cost << " cost2: " << cost2 << "\n";
    }
    // std::cout << delta_t << std::endl;

    double t_0 = (double) delta_t / (double) abs(std::log(parameters.acceptance));


    if (parameters.start_with_nn == 1) {
        std::cout << "temperatura przed T_0 = T_0 * (wynik z NN / losowa trasa): " << t_0 << std::endl;
        t_0 = t_0 * ((double) best_cost / (double) cost);
        std::cout << "temperatura po: " << t_0 << std::endl;
    }

    if (parameters.pocz == 1) {
        t_0 = parameters.temp;
    }



    // początek liczenia czasu
    // this->times.start = std::chrono::steady_clock::now();


    t_0 = abs(t_0);
    std::cout << "T_0 = " << t_0 << "\n";

    // zmienna T - temperatura
    double temperature = t_0;
    // zmienna alfa do geometrycznego schematu chłodzenia (1)
    double multiplier = this->parameters.multiplier;
    std::cout << "multiplier: " << multiplier << std::endl;

    // wybór chłodzenia - geometryczny jest pod koniec algorytmu
    if (this->parameters.cooling_scheme == 2) {
        temperature = t_0 / std::log(1);
    }
    std::cout << "temperature: " << temperature << std::endl;

    // licznik iteracji
    int k = 0;
    bool time_end = true;

    // pomaga przy liczeniu ilości epok od znalezienia lepszego rozwiązania
    bool cost_memory = true;
    int cost_counter = 0;

    int end = 0;

    // główna pętla symulowanego wyżarzania. Działa do momentu wystudzenia lub końca ustawionego czasu
    while (temperature > this->parameters.min_temp && time_end == true) {
        end++;

        int cost = this->CalculateCostOfPath(current_path);
        cost_memory = true;

        // epoka
        for (int i = 0; i < this->parameters.iterations; i++) {
            if (parameters.local_search == 0) {
                new_solution = GeneratePathViaSwap(current_path);
            }
            else if (parameters.local_search == 1) {
                new_solution = GeneratePathViaTwoOpt(current_path);
            }
            else if (parameters.local_search == 2) {
                new_solution = GeneratePathViaInsertion(current_path);
            }

            int new_cost = this->CalculateCostOfPath(new_solution);

            // int delta = cost - new_cost;
            int delta = new_cost - cost;
            // std::cout << "delta: " << delta << std::endl;
            if (delta < 0) {
                current_path = new_solution;

                if (new_cost < best_cost) {
                    // if (this->parameters.show_prints == 1) {
                    //     std::cout << "Found better path, cost: " << new_cost << "\n";
                    //     std::cout << "temperature: " << temperature << std::endl;
                    // }
                    if (this->parameters.show_prints == 3) {
                        std::cout << temperature << "," << new_cost << "\n";
                    }
                    end = 0;
                    best_cost = new_cost;
                    best_path = new_solution;
                    cost_memory = false;
                }
            }
            // kryterium metropolisa
            else {
                float random = ((float) rand() / (RAND_MAX));
                // double random = dis(gen);

                float e = pow(E, -(float)delta/(float)temperature);
                if (random < e) {
                    current_path = new_solution;
                    cost = new_cost;
                }
                if (this->parameters.show_prints == 2) {
                    std::cout << "random: " << random << std::endl;
                    std::cout << "e: " << e << "\n" <<std::endl;
                }
            }

        }


        // ograniczenie czasowe algorytmu
        this->times.time_point = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::seconds>(
                this->times.time_point - this->times.start).count() >= this->parameters.time) {
            time_end = false;
            break;
        }


        // jeśli znaleziono opt to kończy algorytm
        // if (best_cost == this->parameters.optimal_cost) {
        //     this->times.time_point = std::chrono::steady_clock::now();
        //     time_end = false;
        //     break;
        // }

        if (this->parameters.cooling_scheme == 1) {
            // jeśli przez ileś iteracji nie znaleziono rozwiązania to zwiększa temperaturę
            if (cost_memory == false) {
                cost_counter = 0;
            }
            if (cost_counter > this->parameters.iterations_without_improvment) {
                // std::cout << "temp before in counter: " << temperature << std::endl;
                temperature = temperature * 100;
                // std::cout << "temp after in counter: " << temperature << std::endl;
                cost_counter = 0;
            }
            cost_counter++;
        }

        if (parameters.end != 0) {
            if (end >= parameters.end) {
                break;
            }
        }

        // zwiększenie ilości iteracji
        k++;

        // wybór schematu chłodzenia
        if (this->parameters.cooling_scheme == 1) {
            temperature *= multiplier;
        }
        else if (this->parameters.cooling_scheme == 2) {
            temperature = t_0 / std::log(1 + k);
            if (cost_memory == false) {
                cost_counter = 0;
            }
            if (cost_counter > this->parameters.iterations_without_improvment) {
                // std::cout << "temp before in counter: " << temperature << std::endl;
                // k = k / 10;
                // std::cout << "temp after in counter: " << temperature << std::endl;
                cost_counter = 0;
            }
            cost_counter++;
        }
        else {
            std::cout << "cooling schemene not set" << std::endl;
        }
    }
    std::cout << "temperature at end: " << temperature << std::endl;
    this->optimal_path = best_path;
    this->optimal_cost = best_cost;
}

std::vector<int> SimulatedAnnealing::GeneratePathViaTwoOpt(std::vector<int> path) {
    // inwersja (2-opt, zamiana 2-krawędziowa) - losuje dwa indeksy jako początek i koniec inwersji
    path.pop_back();
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
    path.push_back(path[0]);
    return path;
}

std::vector<int> SimulatedAnnealing::GeneratePathViaSwap(std::vector<int> path) {
    path.pop_back();
    int i = rand() % path.size();
    int j = rand() % path.size();
    while (i == j) {
        j = rand() % path.size();
    }
    // std::cout << "i = " << i << "j = " << j << std::endl;
    int tmp = path[i];
    path[i] = path[j];
    path[j] = tmp;
    path.push_back(path[0]);
    return path;
}

std::vector<int> SimulatedAnnealing::GeneratePathViaInsertion(std::vector<int> path) {
    path.pop_back();
    int i = rand() % path.size();
    int j = rand() % path.size();
    // std::cout << "i = " << i << " j = " << j << std::endl;
    while (i == j) {
        j = rand() % path.size();
    }
    if (i < j) {
        int tmp = path[j];
        while (i < j) {
            path[j] = path[j-1];
            j--;
        }
        path[j] = tmp;
    }
    else {
        int tmp = path[i];
        while (j < i) {
            path[i] = path[i-1];
            i--;
        }
        path[i] = tmp;
    }
    path.push_back(path[0]);
    return path;
}


std::vector<int> SimulatedAnnealing::Shuffle(std::vector<int> path) {
    path.pop_back();
    // this->PrintPath(path);
    for (int i = path.size() - 1; i > 0; --i) {
        int random = rand() % (i + 1);
        // swap
        int tmp = path[i];
        path[i] = path[random];
        path[random] = tmp;
    }
    path.push_back(path[0]);
    // this->PrintPath(path);
    return path;
}

std::vector<int> SimulatedAnnealing::ShuffleNTimes(std::vector<int> path) {
    path.pop_back();
    // this->PrintPath(path);
    for (int i = 0; i < std::sqrt(this->node_num); i++) {
        int random_1 = rand() % (this->node_num - 1);
        // swap
        int tmp = path[random_1];
        int random_2 = rand() % (this->node_num - 1) ;
        while (random_2 == random_1) {
            random_2 = std::rand() % (node_num - 1);
        }
        path[random_1] = path[random_2];
        path[random_2] = tmp;
    }
    path.push_back(path[0]);
    // this->PrintPath(path);
    return path;
}

// to jest tak na wszelki wypadek, bo nie wiem czy sam swap z tabu będzie wystarczająco losowy
std::vector<int> SimulatedAnnealing::SwapRandom(std::vector<int> path) {
    this->PrintPath(path);
    path.pop_back();
    int random_1 = rand() % (path.size() - 1);
    // dodać while random_1 != random_2, żeby były zawsze nie takie same wartości ?? i też nie wiem czy
    // nie powinienem rand() zrobić z przedziału, bo pierwsze i ostatnie miasto
    int random_2 = rand() % (path.size() - 1);
    // swap
    int tmp = path[random_1];
    path[random_1] = path[random_2];
    path[random_2] = tmp;
    path.push_back(path[0]);
    this->PrintPath(path);
    return path;
}

// liczy koszt drogi
int SimulatedAnnealing::CalculateCostOfPath(std::vector<int> &path) {
    int cost = 0;
    int current_node, next_node;
    for (int i = 0; i < path.size(); i++) {
        current_node = path.at(i);
        next_node = path.at((i + 1) % path.size());
        cost += this->matrix.at(current_node).at(next_node);
    }
    return cost;
}

// double SimulatedAnnealing::CalculateCostOfPathInEuclideanSpace2D(std::vector<int> &path) {
//     double cost = 0.0;
//
// }


/*
 *  ====================================================================================================================
 */

/*
 *  algorytm NN do określenia Upper bound
 */

void SimulatedAnnealing::NN_UpperBound() {
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

std::vector<int> SimulatedAnnealing::FollowPath(int current_node, std::vector<bool> visited,
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

bool SimulatedAnnealing::NN_CheckIfAllVisited(std::vector<int> &path) {
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

int SimulatedAnnealing::NN_PathCost(const std::vector<int> &path) {
    int cost = 0;
    for (int i = 0; i < path.size(); i++) {
        int current_node = path.at(i);
        int next_node = path.at((i + 1) % path.size());
        cost += this->matrix.at(current_node).at(next_node % path.size());
    }
    return cost;
}

std::vector<int> SimulatedAnnealing::NN_GetALLShortestNodes(int current_node, std::vector<bool> &visited) {
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

int SimulatedAnnealing::ReadTSPFileArray(const std::string &file_name) {
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


int SimulatedAnnealing::ReadConfigFile(const std::string &file_name) {
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

        // start_with_nn
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.start_with_nn = std::stoi(substr);

        // starting_temp
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.temp = std::stof(substr);

        // min_temp
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.min_temp = std::stof(substr);

        // multiplier
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.multiplier = std::stof(substr);

        // podwyzsz_temperature_co:
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.iterations_without_improvment = std::stoi(substr);

        // ilosc_iteracji_bez_poprawy_zakoncz:
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.end = std::stoi(substr);

        // cooling_scheme
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.cooling_scheme = std::stoi(substr);

        // interations - k
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.iterations = std::stoi(substr);

        // local_search
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.local_search = std::stoi(substr);


        // show prints
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.show_prints = std::stoi(substr);

        // _poczatkowa_temperatura:
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.pocz = std::stoi(substr);

        // _poczatkowa_temperatura:
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.acceptance = std::stod(substr);

        file.close();
        return 1;
    }
    std::cout << "Error opening file " << file_name << std::endl;
    return 0;
}



/*
// zmienić wypisywanie w received path, bo problem z + 1
int SimulatedAnnealing::WriteOutputFile() {
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
void SimulatedAnnealing::PrintResult() {
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


void SimulatedAnnealing::PrintPath(const std::vector<int> &path) {
    for (int node : path) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
};

void SimulatedAnnealing::PrintMatrix(std::vector<std::vector<int>> matrix) {
    std::cout << "matrix size: " << matrix.size() << std::endl;
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            std::cout << matrix.at(i).at(j) << " ";
        }
        std::cout << std::endl;
    }
};