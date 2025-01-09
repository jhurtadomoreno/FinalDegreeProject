#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>  
#include <memory>  
#include <stdexcept>  
#include <array>  
#include <omp.h>  
#include <numeric>
#include <chrono>
#include <utility> 
#include <algorithm> 

using namespace std;

std::string vector_to_json(const std::vector<int>& vec) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) {
            oss << ",";
        }
        oss << vec[i];
    }
    oss << "]";
    return oss.str();
}

std::string vector_of_vectors_to_json(const std::vector<std::vector<int>>& vec) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vector_to_json(vec[i]);
        if (i < vec.size() - 1) {
            oss << ",";
        }
    }
    oss << "]";
    return oss.str();
}

std::string vector_of_vector_of_vectors_to_json(const std::vector<std::vector<std::vector<int>>>& vec) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vector_of_vectors_to_json(vec[i]);
        if (i < vec.size() - 1) {
            oss << ",";
        }
    }
    oss << "]";
    return oss.str();
}

std::string vector_to_json(const std::vector<std::pair<float, std::vector<int>>>& vec) {
    std::ostringstream oss;
    oss << "[";

    for (size_t i = 0; i < vec.size(); ++i) {
        const auto& [coef, monomio] = vec[i];

        oss << "[" << coef << ",[";

        for (size_t j = 0; j < monomio.size(); ++j) {
            if (j > 0) {
                oss << ",";
            }
            oss << monomio[j];
        }
        oss << "]]";

        if (i < vec.size() - 1) {
            oss << ",";
        }
    }
    oss << "]";
    return oss.str();
}

std::string vector_of_vectors_to_json(const std::vector<std::vector<std::pair<float, std::vector<int>>>>& vec) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        const auto& inner_vec = vec[i];
        std::cout << "Subvector " << i << ":" << std::endl;
        for (const auto& [coef, monomio] : inner_vec) {
            std::cout << "  Coeficiente: " << coef << ", Monomio: [";
            for (size_t j = 0; j < monomio.size(); ++j) {
                if (j > 0) std::cout << ", ";
                std::cout << monomio[j];
            }
            std::cout << "]" << std::endl;       
	}

	oss << vector_to_json(inner_vec);
        if (i < vec.size() - 1) {
            oss << ",";
        }
    }
    oss << "]";
    return oss.str();
}

std::string run_command(const std::string& command) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    //std::cout << result << std::endl;
    return result;
}

std::vector<std::vector<int>> generate_monomials(int n_variables, int deg) {
    std::vector<std::vector<int>> monomials;
    std::vector<int> current(n_variables, 0);  // Inicializar con ceros
    
    while (true) {
        if (std::accumulate(current.begin(), current.end(), 0) <= deg) {
            monomials.push_back(current);
        }

        int index = n_variables - 1;
        while (index >= 0 && current[index] == deg) {
            --index;
        }

        if (index < 0) {
            break;
        }

        ++current[index];
        std::fill(current.begin() + index + 1, current.end(), 0);
    }

    return monomials;
}

std::vector<std::vector<int>> generate_combinations(int max_deg, int len) {
    std::vector<std::vector<int>> combinations;
    std::vector<int> current(len, 0);  
    
    while (true) {
        combinations.push_back(current);  

        int index = len - 1;
        while (index >= 0 && current[index] == max_deg) {
            --index;
        }

        if (index < 0) {
            break;
        }

        ++current[index];
        std::fill(current.begin() + index + 1, current.end(), 0);
    }

    return combinations;
}

std::vector<std::vector<std::vector<std::vector<int>>>> initialize_vectors(const std::vector<std::vector<int>>& combinations, int n_variables) {
    std::vector<std::vector<std::vector<std::vector<int>>>> result;

    for (const auto& comb : combinations) {
        std::vector<std::vector<std::vector<int>>> vd;

        for (int degree : comb) {
            std::vector<std::vector<int>> monomials = generate_monomials(n_variables, degree);
            vd.push_back(monomials);
        }

        result.push_back(vd);
    }

    return result;
}

void print_combinations(const std::vector<std::vector<int>>& combinations) {
    std::cout << "{";
    for (size_t i = 0; i < combinations.size(); ++i) {
        std::cout << "{";
        for (size_t j = 0; j < combinations[i].size(); ++j) {
            std::cout << combinations[i][j];
            if (j < combinations[i].size() - 1) std::cout << ", ";
        }
        std::cout << "}";
        if (i < combinations.size() - 1) std::cout << ", ";
    }
    std::cout << "}" << std::endl;
}

void print_vectors(const std::vector<std::vector<std::vector<std::vector<int>>>>& vectors) {
    for (size_t i = 0; i < vectors.size(); ++i) {
        std::cout << "vd" << (i + 1) << ":\n";
        for (const auto& vd : vectors[i]) {
            std::cout << "  [";
            for (size_t j = 0; j < vd.size(); ++j) {
                std::cout << "{";
                for (size_t k = 0; k < vd[j].size(); ++k) {
                    std::cout << vd[j][k];
                    if (k < vd[j].size() - 1) std::cout << ", ";
                }
                std::cout << "}";
                if (j < vd.size() - 1) std::cout << ", ";
            }
            std::cout << "]\n";
        }
    }
}

bool cumpleCondicion(const vector<int>& elem, const vector<int>& parametros) {
    int n = elem.size();
    for (int i = 0; i < n; ++i) {
        int producto = 2 * elem[i] + parametros[i];

        if (producto > 2) {
            bool todosMenoresOIguales = true;
            for (int j = 0; j < n; ++j) {
                if (i != j && 2 * elem[j] + parametros[j] >= producto) {
                    todosMenoresOIguales = false;
                    break;
                }
            }
            if (todosMenoresOIguales) {
                return true; 
            }
        }
    }
    return false; 
}

void filtrarVector(vector<vector<int>>& vec, const vector<int>& parametros) {
    vec.erase(remove_if(vec.begin(), vec.end(),
                        [&](const vector<int>& elem) { return cumpleCondicion(elem, parametros); }),
              vec.end());
}

int calcularGradoMonomio(const vector<int>& monomio) {
    int grado = 0;
    for (int exponente : monomio) {
        grado += exponente;
    }
    return grado;
}

int calcularGradoPolinomio(const vector<pair<float, vector<int>>>& polinomio) {
    int gradoMax = 0;
    for (const auto& termino : polinomio) {
        int gradoTermino = calcularGradoMonomio(termino.second);
        if (gradoTermino > gradoMax) {
            gradoMax = gradoTermino;
        }
    }
    return gradoMax;
}

int main() {
    /*
    std::vector<std::pair<float, std::vector<int>>> f = {
        {2, {0, 0}}, {-1, {2, 0}}, {-1, {0, 2}}
    };
    std::vector<std::vector<std::pair<float, std::vector<int>>>> g = {
        {},
        {{1, {3, 0}}, {-1, {0, 2}}},
        {{1, {0, 0}}, {-1, {1, 0}}},
    };*/
   /*
    std::vector<std::pair<float, std::vector<int>>> f = { 
        {1, {0, 0}}, {-1, {2, 0}}, {-1, {0, 2}} 
    };
    std::vector<std::vector<std::pair<float, std::vector<int>>>> g = {
        {},
        {{1, {1, 0}}},
        {{1, {0, 1}}},
        {{1, {0, 0}}, {-1, {1, 0}}, {-1, {0, 1}}},
    };*/
    /*	
    std::vector<std::pair<float, std::vector<int>>> f = {
        {2, {0, 0, 0}}, {-1, {2, 0, 0}}, {-1, {0, 2, 0}}
    };
    std::vector<std::vector<std::pair<float, std::vector<int>>>> g = {
        {},
        {{1, {0, 0, 0}}, {1, {1, 0, 2}}, {-1, {2, 0, 0}}, {-1, {0, 2, 0}}},
	{{1, {0, 0, 0}}, {-1, {1, 0, 2}}}
    };*/
    /*std::vector<std::pair<float, std::vector<int>>> f = {
        {5, {0, 0, 0}}, {-1, {2, 0, 0}}, {-1, {0, 2, 0}}, {-1, {0, 0, 2}}
    };
    std::vector<std::vector<std::pair<float, std::vector<int>>>> g = {
        {},
        {{1, {1, 0, 0}}},
        {{1, {0, 1, 0}}},
	{{1, {0, 0, 1}}}
    };*/

    /*
    std::vector<std::pair<int, std::vector<int>>> f = { 
        {1, {0, 0, 0}}, {-1, {2, 0, 0}}, {-1, {0, 2, 0}}, {-1, {0, 0, 2}} 
    };
    std::vector<std::vector<std::pair<int, std::vector<int>>>> g = {
        {},
        {{1, {1, 0, 0}}},
        {{1, {0, 1, 0}}},
	{{1, {0, 0, 1}}},
        {{1, {0, 0, 0}}, {-1, {1, 0, 0}}, {-1, {0, 1, 0}}, {-1, {0, 0, 1}}},
    };*/
    
    std::vector<std::pair<float, std::vector<int>>> f = { 
        {4, {0, 0}}, {-1, {2, 0}}, {-1, {0, 2}} 
    };
    std::vector<std::vector<std::pair<float, std::vector<int>>>> g = {
        {},
        {{1, {1, 0}}, {-1.0f / 2.0f, {0, 0}}},
        {{1, {0, 1}}, {-1.0f / 2.0f, {0, 0}}},
        {{1, {0, 0}}, {-1, {1, 1}}},
    };
    
    /*	
    std::vector<std::pair<int, std::vector<int>>> f = { 
        {17, {0, 0, 0}}, {-1, {2, 0, 0}}, {-1, {0, 2, 0}}, {-1, {0, 0, 2}} 
    };
    std::vector<std::vector<std::pair<int, std::vector<int>>>> g = {
        {},
        {{1, {1, 0, 0}}, {-1/2, {0, 0, 0}}},
        {{1, {0, 1, 0}}, {-1/2, {0, 0, 0}}},
	{{1, {0, 0, 1}}, {-1/2, {0, 0, 0}}},
        {{1, {0, 0, 0}}, {-1, {1, 1, 1}}},
    };
	*/

    int max_deg = 4;   
    int len = g.size();
    int n_variables = 2; 
    bool found = false;

    std::vector<std::vector<int>> combinations = generate_combinations(max_deg, len);
    //print_combinations(combinations);
    std::cout << combinations.size() << std::endl;

    vector<int> grados;

    for (const auto& polinomio : g) {
        int grado = polinomio.empty() ? 0 : calcularGradoPolinomio(polinomio);
        grados.push_back(grado);
    }

    cout << "Grados de los polinomios: {";
    for (size_t i = 0; i < grados.size(); ++i) {
        cout << grados[i];
        if (i < grados.size() - 1) cout << ", ";
    }
    cout << "}" << endl;

    filtrarVector(combinations, grados);
    std::cout << combinations.size() << std::endl;

    std::vector<std::vector<std::vector<std::vector<int>>>> vectors = initialize_vectors(combinations, n_variables);

    std::string f_str = vector_to_json(f);
    std::string g_str = vector_of_vectors_to_json(g);

    auto start_time = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for num_threads(5) 
    for (size_t i = 0; i < vectors.size(); ++i) {
        #pragma omp cancellation point for
        std::string vd_str = vector_of_vector_of_vectors_to_json(vectors[i]);
        std::string command = "python3 QuadraticModule.py " + f_str + " " + g_str + " " + vd_str;
        
        try {
            std::string output = run_command(command);
            int number = std::stoi(output);

            int thread_id = omp_get_thread_num();
            std::cout << "Thread " << thread_id << " Output for vd[" << i << "]: " << output << std::endl;
            
            if (number == 1){
                std::cout << "Opció vàlida" << std::endl;
		std::cout << vd_str << std::endl;
		found = true;
		#pragma omp cancel for
            }
            
        } catch (const std::runtime_error& e) {
            int thread_id = omp_get_thread_num();
            std::cerr << "Thread " << thread_id << " Encountered an error: " << e.what() << std::endl;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    std::cout << "Tiempo de ejecución: " << elapsed.count() << " segundos." << std::endl;
    std::cout << found << std::endl;

    return 0;
}
