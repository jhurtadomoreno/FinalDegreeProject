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

/*
Methods for converting vectors to JSON (these vectors represent the degree combinations for the sum of squares polynomials in the decomposition, which is why they are vectors of integers).
*/
string vector_to_json(const vector<int>& vec) {
    ostringstream oss;
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

string vector_of_vectors_to_json(const vector<vector<int>>& vec) {
    ostringstream oss;
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

string vector_of_vector_of_vectors_to_json(const vector<vector<vector<int>>>& vec) {
    ostringstream oss;
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

/*
Methods for converting vectors to JSON (these vectors represent polynomials, which is why they are vectors of coefficient-monomial tuples  in the case of polynomial f, and vector of vectors of coefficient-monomial tuples, i.e. vector of polynomials, in the case of g=g₁, ..., gₛ).
*/
string vector_to_json(const vector<pair<float, vector<int>>>& vec) {
    ostringstream oss;
    oss << "[";

    for (size_t i = 0; i < vec.size(); ++i) {
        const auto& [coef, monomial] = vec[i];

        oss << "[" << coef << ",[";

        for (size_t j = 0; j < monomial.size(); ++j) {
            if (j > 0) {
                oss << ",";
            }
            oss << monomial[j];
        }
        oss << "]]";

        if (i < vec.size() - 1) {
            oss << ",";
        }
    }
    oss << "]";
    return oss.str();
}

string vector_of_vectors_to_json(const vector<vector<pair<float, vector<int>>>>& vec) {
    ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        const auto& inner_vec = vec[i];
	
	oss << vector_to_json(inner_vec);
        if (i < vec.size() - 1) {
            oss << ",";
        }
    }
    oss << "]";
    return oss.str();
}

/*
Method for communicating C++ and Python
*/
string run_command(const string& command) {
    array<char, 128> buffer;
    string result;
    unique_ptr<FILE, decltype(&pclose)> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe) {
        throw runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    
    return result;
}

/*
Method for generating all monomials of n_variables variables with a degree less than or equal to deg.
*/
vector<vector<int>> generate_monomials(int n_variables, int deg) {
    vector<vector<int>> monomials;
    vector<int> current(n_variables, 0); 
    
    while (true) {
        if (accumulate(current.begin(), current.end(), 0) <= deg) {
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
        fill(current.begin() + index + 1, current.end(), 0);
    }

    return monomials;
}

/*
Method for generating all degree combinations of the polynomials σ₁, ..., σₛ. Combinations of as many numbers as len, with a degree less than or equal to max_deg for each individual one. Unlike the previous method, the degree restriction is individual, not based on the monomial (the sum). These could be unified into a single method with a third boolean parameter, but this approach is clearer.
*/
vector<vector<int>> generate_combinations(int len, int max_deg) {
    vector<vector<int>> combinations;
    vector<int> current(len, 0);  
    
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
        fill(current.begin() + index + 1, current.end(), 0);
    }

    return combinations;
}

/*
Method that returns the vectors vd₁, ..., vdₛ for all degree combinations σ₁, ..., σₛ. Specifically, it produces a structure of the form [a, b, c, ...], where each letter represents a vector containing the vectors vd₁, ..., vdₛ associated with a specific degree combination. Each vector vdᵢ is a vector of integer vectors representing the monomials of n variables with a degree less than or equal to a given value (e.g. vdᵢ = [[0,0], [0,1], [1,0], [1,1], ...], where each element of the outer vector corresponds to a monomial).
*/
vector<vector<vector<vector<int>>>> initialize_vectors(const vector<vector<int>>& combinations, int n_variables) {
    vector<vector<vector<vector<int>>>> result;

    for (const auto& comb : combinations) {
        vector<vector<vector<int>>> vd;

        for (int degree : comb) {
            vector<vector<int>> monomials = generate_monomials(n_variables, degree);
            vd.push_back(monomials);
        }

        result.push_back(vd);
    }

    return result;
}

/*
Methods for printing
*/
void print_combinations(const vector<vector<int>>& combinations) {
    cout << "{";
    for (size_t i = 0; i < combinations.size(); ++i) {
        cout << "{";
        for (size_t j = 0; j < combinations[i].size(); ++j) {
            cout << combinations[i][j];
            if (j < combinations[i].size() - 1) cout << ", ";
        }
        cout << "}";
        if (i < combinations.size() - 1) cout << ", ";
    }
    cout << "}" << endl;
}

void print_vectors(const vector<vector<vector<vector<int>>>>& vectors) {
    for (size_t i = 0; i < vectors.size(); ++i) {
        cout << "vd" << (i + 1) << ":\n";
        for (const auto& vd : vectors[i]) {
            cout << "  [";
            for (size_t j = 0; j < vd.size(); ++j) {
                cout << "{";
                for (size_t k = 0; k < vd[j].size(); ++k) {
                    cout << vd[j][k];
                    if (k < vd[j].size() - 1) cout << ", ";
                }
                cout << "}";
                if (j < vd.size() - 1) cout << ", ";
            }
            cout << "]\n";
        }
    }
}

/*
Methods for filtering the degree combinations of the sigmas that can be discarded
*/
bool meetsCondition(const vector<int>& combination, const vector<int>& degs_g) {
    int n = combination.size();
    for (int i = 0; i < n; ++i) {
        int prod = 2 * combination[i] + degs_g[i];

        if (prod > 2) {
            bool allBelowOrEqual = true;
            for (int j = 0; j < n; ++j) {
                if (i != j && 2 * combination[j] + degs_g[j] >= prod) {
                    allBelowOrEqual = false;
                    break;
                }
            }
            if (allBelowOrEqual) {
                return true; 
            }
        }
    }
    return false; 
}

void filterCombinations(vector<vector<int>>& combinations, const vector<int>& degs_g) {
    combinations.erase(remove_if(combinations.begin(), combinations.end(),
                        [&](const vector<int>& comb) { return meetsCondition(comb, degs_g); }),
              combinations.end());
}

/*
Methods for calculating the degree of a polynomial
*/
int calculateMonomialDegree(const vector<int>& monomial) {
    int degree = 0;
    for (int exponent : monomial) {
        degree += exponent;
    }
    return degree;
}

int calculatePolynomialDegree(const vector<pair<float, vector<int>>>& polynomial) {
    int maxDegree = 0;
    for (const auto& term : polynomial) {
        int termDegree = calculateMonomialDegree(term.second);
        if (termDegree > maxDegree) {
            maxDegree = termDegree;
        }
    }
    return maxDegree;
}

int main() {
    /* ---------------- Examples of input for the polynomials f and g₁, ..., gₛ ----------------
    (1)
    vector<pair<float, vector<int>>> f = {
        {2, {0, 0}}, {-1, {2, 0}}, {-1, {0, 2}}
    };
    vector<vector<pair<float, vector<int>>>> g = {
        {},
        {{1, {3, 0}}, {-1, {0, 2}}},
        {{1, {0, 0}}, {-1, {1, 0}}},
    };
    
    (2)
    vector<pair<float, vector<int>>> f = { 
        {1, {0, 0}}, {-1, {2, 0}}, {-1, {0, 2}} 
    };
    vector<vector<pair<float, vector<int>>>> g = {
        {},
        {{1, {1, 0}}},
        {{1, {0, 1}}},
        {{1, {0, 0}}, {-1, {1, 0}}, {-1, {0, 1}}},
    };
    
    (3)
    vector<pair<float, vector<int>>> f = {
        {2, {0, 0, 0}}, {-1, {2, 0, 0}}, {-1, {0, 2, 0}}
    };
    vector<vector<pair<float, vector<int>>>> g = {
        {},
        {{1, {0, 0, 0}}, {1, {1, 0, 2}}, {-1, {2, 0, 0}}, {-1, {0, 2, 0}}},
	{{1, {0, 0, 0}}, {-1, {1, 0, 2}}}
    };
    
    (4)
    vector<pair<float, vector<int>>> f = {
        {5, {0, 0, 0}}, {-1, {2, 0, 0}}, {-1, {0, 2, 0}}, {-1, {0, 0, 2}}
    };
    vector<vector<pair<float, vector<int>>>> g = {
        {},
        {{1, {1, 0, 0}}},
        {{1, {0, 1, 0}}},
	{{1, {0, 0, 1}}}
    };
    
    (5)
    vector<pair<float, vector<int>>> f = { 
        {1, {0, 0, 0}}, {-1, {2, 0, 0}}, {-1, {0, 2, 0}}, {-1, {0, 0, 2}} 
    };
    vector<vector<pair<float, vector<int>>>> g = {
        {},
        {{1, {1, 0, 0}}},
        {{1, {0, 1, 0}}},
	{{1, {0, 0, 1}}},
        {{1, {0, 0, 0}}, {-1, {1, 0, 0}}, {-1, {0, 1, 0}}, {-1, {0, 0, 1}}},
    };
    
    (6)
    vector<pair<float, vector<int>>> f = { 
        {17, {0, 0, 0}}, {-1, {2, 0, 0}}, {-1, {0, 2, 0}}, {-1, {0, 0, 2}} 
    };
    vector<vector<pair<float, vector<int>>>> g = {
        {},
        {{1, {1, 0, 0}}, {-1/2, {0, 0, 0}}},
        {{1, {0, 1, 0}}, {-1/2, {0, 0, 0}}},
	{{1, {0, 0, 1}}, {-1/2, {0, 0, 0}}},
        {{1, {0, 0, 0}}, {-1, {1, 1, 1}}},
    };
    */
    
    // ---------------- Input to be modified by the user ----------------
    int max_deg = 2; // Maximum degree for the vdᵢ's in the decomposition (the maximum degree of the σᵢ's would be 2 * max_deg)
    vector<pair<float, vector<int>>> f = { 
        {4, {0, 0}}, {-1, {2, 0}}, {-1, {0, 2}} 
    };
    vector<vector<pair<float, vector<int>>>> g = {
        {},
        {{1, {1, 0}}, {-1.0f / 2.0f, {0, 0}}},
        {{1, {0, 1}}, {-1.0f / 2.0f, {0, 0}}},
        {{1, {0, 0}}, {-1, {1, 1}}},
    };
    // ------------------------------------------------------
    
    int n_variables = f.empty() ? 0 : f[0].second.size();
    int len = g.size();
    bool found = false;
    vector<int> degs_g;
    string f_str = vector_to_json(f);
    string g_str = vector_of_vectors_to_json(g);
    
    // Generate all possible degree combinations for the vdᵢ
    vector<vector<int>> combinations = generate_combinations(len, max_deg);
    cout << combinations.size() << endl;
    //print_combinations(combinations);
    
    // Calculate the degrees of g₁, ..., gₛ and display them
    for (const auto& polinomio : g) {
        int grado = polinomio.empty() ? 0 : calculatePolynomialDegree(polinomio);
        degs_g.push_back(grado);
    }

    cout << "Degrees of the polynomials gᵢ: {";
    for (size_t i = 0; i < degs_g.size(); ++i) {
        cout << degs_g[i];
        if (i < degs_g.size() - 1) cout << ", ";
    }
    cout << "}" << endl;
    
    // Filter the obtained combinations
    filterCombinations(combinations, degs_g);
    cout << combinations.size() << endl;
    
    // Convert degree combinations of each vdᵢ to the monomial vectors representing them
    vector<vector<vector<vector<int>>>> vectors = initialize_vectors(combinations, n_variables);

    auto start_time = chrono::high_resolution_clock::now();
    
    // Iterate through all degree combinations to check if any allows for a sufficiently accurate representation
    #pragma omp parallel for num_threads(5) // Adjust the number of threads 
    for (size_t i = 0; i < vectors.size(); ++i) {
        #pragma omp cancellation point for
        string vd_str = vector_of_vector_of_vectors_to_json(vectors[i]);
        string command = "python3 QuadraticModule.py " + f_str + " " + g_str + " " + vd_str;
        
        try {
            string output = run_command(command); // Execute the Python code to search for a representation
            int number = stoi(output);

            int thread_id = omp_get_thread_num();
            cout << "Thread " << thread_id << " Output for vd[" << i << "]: " << output << endl;
            
            if (number == 1){ // If we have found an accurate solution
                cout << "Accurate solution" << endl;
		cout << vd_str << endl;
		found = true;
		#pragma omp cancel for
            }
            
        } catch (const runtime_error& e) {
            int thread_id = omp_get_thread_num();
            cerr << "Thread " << thread_id << " Encountered an error: " << e.what() << endl;
        }
    }

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    cout << "Execution time: " << elapsed.count() << " seconds." << endl;
    cout << found << endl;

    return 0;
}
