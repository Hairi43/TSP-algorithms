#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <xlnt/xlnt.hpp>
#include <random>

#define MAX_EDGE 9999

struct Parameters {
    std::string file_name;
    std::string output_file_name;
    std::string type;
    int time;
    int iterate;
    int optimal_cost;
    int tabu_tenure; // kadencja - na jak długo zamieniony wierzchołek ma być na liście tabu
    std::vector<int> optimal_path;
    int starting_city = 3;
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

class TabuSearch {
    int node_num;
    int optimal_cost;
    int bound; // potrzebne do nn
    std::vector<int> optimal_path;
    std::vector<std::vector<int>> matrix;
    Parameters parameters;
    Times times;
    NN_data nn_data;

public:
    TabuSearch(const std::string &file_name);

    void TS();

    std::vector<int> Shuffle(std::vector<int> path);

    std::vector<int> ShuffleNTimes(std::vector<int> path);

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

    TabuSearch("config.txt");
    

    std::cout << "\n\nPress enter to exit" << std::endl;
    std::cin.ignore();
    return 0;
}


TabuSearch::TabuSearch(const std::string &file_name) {
    if (ReadConfigFile(file_name) && ReadTSPFileArray(this->parameters.file_name)) {
        this->optimal_cost = INT_MAX;

        this->times.start = std::chrono::steady_clock::now();
        // this->times.start_high_res = std::chrono::high_resolution_clock::now();

        // obliczenie startowego rozwiązania za pomocą NN, które później będzie przeszukiwane tabu search
        this->NN_UpperBound();
        this->bound = this->nn_data.min_cost;

        std::cout << "Initialized path cost: " << this->bound << std::endl;
        std::cout << "Initialized path: \n";
        this->PrintPath(this->nn_data.best_path);
        std::cout << "\n";

        // std::vector<int> tmp;
        // for (int i = 0; i < 10; i++) {
        // tmp = this->Shuffle(this->nn_data.best_path);
        //     // this->PrintPath(tmp);
        // }

        this->times.time_point = std::chrono::steady_clock::now();
        std::cout << "Time taken for calculating UB: " << std::chrono::duration_cast<std::chrono::duration<double>>(
            this->times.time_point - this->times.start).count() << "s" << std::endl;


        // odkomentować jeśli chcę zmierzyć czas osobno dla liczenia UB i algorytmu TabuSearch
        // this->times.start = std::chrono::steady_clock::now();

        /*
         *  główna funkcja
         *
         */
        this->TS();

        // time measurment
        this->times.end = std::chrono::steady_clock::now();
        this->times.elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(this->times.end - this->times.start).count();

        // this->WriteOutputFile();
        this->PrintResult();
    }
}

void TabuSearch::TS() {
    // dodać kryterium stopu do configa - liczba iteracji it
    int it = this->parameters.iterate;

    int tenure = this->parameters.tabu_tenure; // tabu tunure - number of iteration for move will be kept tabu where
    //  t= sqrt(n)

    /*
     *  może trzeba będzie zmienić z tego z rekurencją na zwykły, żeby raz puszczać z jednego wierzchołka, a raz z innego,
     *  bo inaczej trzeba będzie losować trasę jak w losowym.
     */

    std::vector<int> current_path = this->nn_data.best_path;
    std::vector<int> best_path = this->nn_data.best_path;
    int best_cost = this->nn_data.min_cost;
    // best_cost = MAX_EDGE;

    // lista tabu
    std::vector<std::vector<int>> tabu_list;
    int counter = 0;

    // dodane do sprawdzenia w ile czasu znalazł opt jeśli jest wiadomy
    bool found = false;
    bool diversification = true;

    bool time_end = true;
    int max_tries = 0;
    // tabu search może działać w nieskończoność, więc jest tylko ograniczenie czasowe ====================
    // while (time_end && max_tries < 4) {
    while (time_end) {

        bool go_again = true;
        bool intensification = false;

        // reset licznika iteracji
        counter = 0;

        // for (int counter = 0; counter < this->parameters.iterate; counter++) {
        // while (counter < this->parameters.iterate) {
        while (go_again) {

            std::vector<int> new_path;
            int cost = MAX_EDGE;

            // tymczasowa pamięć dla danej podmiany, żeby zapamiętać zamienione miasta i ich kadencję
            std::vector<int> tabu_swap(3, 0);


            // tworzone jest lokalne sąsiedzctwo, żeby przejrzeć potencjalne rozwiązania i wybrać najlepsze (kombinacje ścieżki
            // danej przez algorytm NN)

            /*  zmienić na i = 1, bo pierwszy też permutuje... ============================================================ */
            // zmienić na i = 1, bo pierwsze miasto nie powinno podlegać permutacji, ale działa też dla i = 0 ?
            for (int i = 1; i < this->node_num - 1; i++) {
                for (int j = i + 1; j < this->node_num; j++) {
                    // std::cout << "i = " << i << std::endl;
                    // std::cout << "j = " << j << std::endl;


                    // // //  podmiana sąsiadów
                    std::vector<int> swapped_path = current_path;
                    // int tmp = swapped_path[i];
                    // swapped_path[i] = swapped_path[j];
                    // swapped_path[j] = tmp;

                    // 2-opt
                    int k = i;
                    int l = j;
                    while (k < l) {
                        int temp = swapped_path[k];
                        swapped_path[k] = swapped_path[l];
                        swapped_path[l] = temp;
                        ++k;
                        --l;
                    }

                    // this->PrintPath(swapped_path);

                    // koszt nowej ścieżki
                    int swap_cost;
                    swap_cost = this->CalculateCostOfPath(swapped_path);

                    // sprawdzenie czy ruch nie jest zakazany
                    bool tabu = false;
                    for (int k = 0; k < tabu_list.size(); k++) {
                        if ((tabu_list[k][0] == i && tabu_list[k][1]) == j || (tabu_list[k][0] == j && tabu_list[k][1] == i)) {
                            tabu = true;
                            break;
                        }
                    }

                    // jeśli rozwiązanie aspiruje ????
                    // if (swap_cost < cost) {
                    if (swap_cost < cost || (tabu && swap_cost <= best_cost)) {
                        new_path = swapped_path;
                        cost = swap_cost;
                        tabu_swap[0] = i;
                        tabu_swap[1] = j;
                    }
                }
            }
            // zapisanie najlepszej aktualnie ścieżki
            current_path = new_path;
            counter++;

            tabu_swap[2] = tenure;
            // sprawdzenie listy tabu, żeby zmniejszyć kadencję danej zmiany
            for (int i = 0; i < tabu_list.size(); i++) {
                tabu_list[i][2] -= 1;
                if (tabu_list[i][2] == 0) {
                    tabu_list.erase(tabu_list.begin() + i);
                    // i--;
                }
            }
            // dodanie ostatniego przesunięcia do listy tabu
            tabu_list.push_back(tabu_swap);

            // według Johna Knox'a lista tabu ma nie przekreczać rozmiaru 3*n
            // - str. 164 "Jak to rozwiązać, czyli nowoczesna heurystyka"
            if (tabu_list.size() > 3 * this->node_num) {
                tabu_list.erase(tabu_list.begin());
            }

            if (cost < best_cost) {
                best_path = new_path;
                best_cost = cost;
                this->PrintPath(best_path);
                std::cout << "path^ cost: " << best_cost << std::endl;

                intensification = true;
                counter = 0;
            }

            // wypisuje czas jeśli znalazło opt. Tak dla ciekawości
            if (best_cost == this->parameters.optimal_cost && found == false) {
                found = true;
                this->times.time_point = std::chrono::steady_clock::now();
                std::cout << "found opt in: " << std::chrono::duration_cast<std::chrono::duration<double>>(
            this->times.time_point - this->times.start).count() << "s" << std::endl;
            }

            // ograniczenie czasowe algorytmu
            this->times.time_point = std::chrono::steady_clock::now();
            if (std::chrono::duration_cast<std::chrono::minutes>(
                    this->times.time_point - this->times.start).count() >= this->parameters.time)
            {
                time_end = false;
                break;
            }


            if (diversification) {
                if (intensification) {
                    current_path = best_path;
                    tenure = this->parameters.tabu_tenure / 4;
                    // std::cout << "tenure: " << tenure << std::endl;
                    intensification = false;
                }
                // lepsze wyniki jeśli dopiero po wielu iteracjach przejdzie do następnego minimum lokalnego
                else if (counter >= this->parameters.iterate) {
                    // current_path = algorytm, który daje losową/nową ścieżkę
                    tenure = this->parameters.tabu_tenure;
                    intensification = false;
                    // std::cout << "completly new path" << std::endl;
                    current_path = this->ShuffleNTimes(current_path);
                    // this->PrintPath(current_path);
                }
            }
        }
        // zwiększenie licznika ogólnych iteracji
        // max_tries++;
    }
    this->optimal_path = best_path;
    this->optimal_cost = best_cost;
}

std::vector<int> TabuSearch::Shuffle(std::vector<int> path) {
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

std::vector<int> TabuSearch::ShuffleNTimes(std::vector<int> path) {
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

// liczy koszt drogi
int TabuSearch::CalculateCostOfPath(std::vector<int> &path) {
    int cost = 0;
    int current_node, next_node;
    for (int i = 0; i < path.size(); i++) {
        current_node = path.at(i);
        next_node = path.at((i + 1) % path.size());
        cost += this->matrix.at(current_node).at(next_node);
    }
    return cost;
}


/*
 *  ====================================================================================================================
 */

/*
 *  algorytm NN do określenia Upper bound
 */

 void TabuSearch::Upperbound() {
    std::vector<bool> visited (this->node_num, false);
    std::vector<int> path;
    path.push_back(0);
    visited[0] = true;
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

void TabuSearch::NN_UpperBound() {
    std::vector<int> path;
    std::vector<bool> visited;
    visited.resize(this->node_num, false);
    this->nn_data.min_cost = std::numeric_limits<int>::max();
    int min = std::numeric_limits<int>::max();
    std::vector<int> optimal_path;

    for (int i = 0; i < this->node_num; i++) {
        std::vector<int> result_path = this->FollowPath(i, visited, path);

        if (this->nn_data.min_cost < min && this->nn_data.min_cost != 0) {
            min = this->nn_data.min_cost;
            optimal_path = this->nn_data.best_path;
        }
    }

    this->nn_data.best_path = optimal_path;
    this->nn_data.min_cost = min;
}

std::vector<int> TabuSearch::FollowPath(int current_node, std::vector<bool> visited,
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

    for (int city : best_paths) {
        std::vector<int> new_path = this->FollowPath(city, visited, current_path);

        int cost = this->NN_PathCost(new_path);

        if (cost < this->optimal_cost && cost != 0) {
            this->nn_data.min_cost = cost;
            this->nn_data.best_path = new_path;
            // best_path = new_path;
        }
    }

    return best_path;

}

bool TabuSearch::NN_CheckIfAllVisited(std::vector<int> &path) {
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

int TabuSearch::NN_PathCost(const std::vector<int> &path) {
    int cost = 0;
    for (int i = 0; i < path.size(); i++) {
        int current_node = path.at(i);
        int next_node = path.at((i + 1) % path.size());
        cost += this->matrix.at(current_node).at(next_node % path.size());
    }
    return cost;
}

std::vector<int> TabuSearch::NN_GetALLShortestNodes(int current_node, std::vector<bool> &visited) {
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

int TabuSearch::ReadTSPFileArray(const std::string &file_name) {
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


int TabuSearch::ReadConfigFile(const std::string &file_name) {
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

        // iterate
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.iterate = std::stoi(substr);

        // time
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.time = std::stoi(substr);

        // tabu tenure
        std::getline(file, line);
        substr = line.substr(line.find(' ') + 1);
        this->parameters.tabu_tenure = std::stoi(substr);

        file.close();
        return 1;
    }
    std::cout << "Error opening file " << file_name << std::endl;
    return 0;
}




// zmienić wypisywanie w received path, bo problem z + 1
int TabuSearch::WriteOutputFile() {
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
    ws.cell("D2").value(this->optimal_cost);
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
void TabuSearch::PrintResult() {
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

    int x = this->parameters.optimal_cost - this->optimal_cost;
    std::cout << "\n\nAbsolute error: " << abs(x) << "\n";
    float y = (float(x)/float(this->parameters.optimal_cost)) * 100;
    std::cout << "Relative error: " << abs(y) << "%";

}


void TabuSearch::PrintPath(const std::vector<int> &path) {
    for (int node : path) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
};

void TabuSearch::PrintMatrix(std::vector<std::vector<int>> matrix) {
    std::cout << "matrix size: " << matrix.size() << std::endl;
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            std::cout << matrix.at(i).at(j) << " ";
        }
        std::cout << std::endl;
    }
};