import itertools, time, cvxpy as cp, numpy as np, sympy as sp, symengine as se
#import mosek
from sympy import symbols, expand
x1, x2, x3, x4, x5, x6, x7 = symbols('x1 x2 x3 x4 x5 x6 x7')

"""
Method to check if a matrix is SDP (Semidefinite Positive)
"""
def sdp(matrix, threshold=1e-4): 
    if matrix.shape[0] != matrix.shape[1]: 
        return False
    
    eigenvalues, _ = np.linalg.eig(matrix)
    eigenvalues[np.abs(eigenvalues) < threshold] = 0

    return np.all(eigenvalues >= 0)

"""
Method to decompose the product of the matrix vd * vd' into the different monomials multiplied by specific matrices
"""
def decomposition(matrix):
    no_rep_monom = set()
    desc_matrices = list()

    for row in matrix:
        for element in row:
            no_rep_monom.add(tuple(element))

    for mon in no_rep_monom:
        desc_matrix = []
        for row in matrix:
            new_row = []
            for element in row:
                if tuple(element) == mon:
                    new_row.append(1)
                else:
                    new_row.append(0)
            desc_matrix.append(new_row)

        desc_matrices.append(desc_matrix)

    return list(no_rep_monom), desc_matrices

"""
Method to expand algebraic expressions in the sum of squares decomposition
"""
def expand_expression(matrix, vd):
    matrix = np.array(matrix)
    
    num_vars = len(vd[0])
    variables = sp.symbols(f'x1:{num_vars+1}')
    
    variable_terms = []
    for exponents in vd:
        if all(e == 0 for e in exponents):
            term = 1  # (0,0) is 1
        else:
            term = sp.Mul(*(variables[j]**exponents[j] for j in range(num_vars) if exponents[j] != 0))
        variable_terms.append(term)
    
    expression = 0
    for col in range(matrix.shape[1]):
        term = 0
        for row in range(matrix.shape[0]):
            term += matrix[row, col] * variable_terms[row]  
        expression += term**2  

    return se.expand(expression)
    
"""
Method to decompose a matrix Q = H * H' using Singular Value Decomposition (SVD)
"""
def prod_decomposition(matrix):
    U, s, Vt = np.linalg.svd(matrix) 
    s_sqrt = np.sqrt(s)
    S_sqrt = np.diag(s_sqrt)
    
    H = np.dot(U, S_sqrt)

    """
    # Check that Q = H * H^T
    reconstructed_Q = np.dot(H, H.T)
    
    print("\nReconstructed Q:")
    print(reconstructed_Q)
    print("\nH:")
    print(H)
    """
    
    return H

"""
Method to check if the obtained solution is sufficiently accurate
"""
def check_output(polynomial, tol = 0.00001):
    f_monomials = [term[1] for term in f]
    out_result = 1
    for term in f:
        found = False
        for element in polynomial:
            if element[1] == term[1]:
                found = True
                if abs(element[0]-term[0]) >= tol: 
                    out_result = 0
                break
        
        if out_result == 0:
            break

        if found == False:
            out_result = 0
            break
    
    for element in polynomial:
        if element[1] not in f_monomials and abs(element[0]) >= tol: 
            out_result = 0
            break

    return out_result
    
"""
Method for generating all monomials of n_variables variables with a degree less than or equal to deg.
"""
def generate_monomials(n_variables, deg):
    monomials = []
    for comb in itertools.product(range(deg + 1), repeat=n_variables):
        if sum(comb) <= deg and comb not in monomials:
            if n_variables == 1:
                monomials.append((comb[0], 0))
            else:
                monomials.append(comb)
    return monomials
 
"""
Method to search for a sum of squares decomposition
"""
def SOS(f, vd, num_threads=1):
    Q = cp.Variable((len(vd), len(vd)), symmetric=True)
    constraints = [Q >> 0]

    # Matrix of monomials resulting from vd * vd'
    prod_vd = [[tuple(a[i] + c[i] for i in range(len(vd[0]))) for c in vd] for a in vd]

    # Decompose into the matrices B_alpha
    l1, l2 = decomposition(prod_vd)

    # Set tr(Q, B_alpha) = g_alpha
    for monomial in l1:
        is_in = False
        for coef, monom in f:
            if monomial == monom:
                constraints += [cp.trace(Q @ l2[l1.index(monomial)]) == coef]
                is_in = True
                break

        if not is_in:
            constraints += [cp.trace(Q @ l2[l1.index(monomial)]) == 0]

    prob = cp.Problem(cp.Minimize(0), constraints)

    try:
        #prob.solve(solver='MOSEK', verbose=False , mosek_params = {'MSK_IPAR_NUM_THREADS':  num_threads})
        prob.solve(solver='SCS')

        if Q.value is not None:
            pass
            """
            print(sdp(Q.value))
            print("....................................................................")
            print(Q.value)
            """
        else:
            print("No sum of squares decomposition found.")

        return Q.value

    except cp.error.SolverError as e:
        print("Error with the solver:", e)
        return None
        
"""
Method responsible for the general logic
"""
def main(f):
    univariate = all(len(monomial) == 2 and monomial[1] == 0 for _, monomial in f)
    if univariate:
        num_vars = 1
    else:
        num_vars = len(f[0][1])
    deg_f = max(sum(monomial) for _, monomial in f)
    vd = generate_monomials(num_vars, deg_f // 2)

    Q = SOS(f, vd)

    # We check if we have obtained a solution and, if so, whether it is sufficiently accurate
    if Q is not None:
        H = prod_decomposition(Q)
        expr = expand_expression(np.matrix(H), vd)

        expr = sp.expand(expr)
        terms = expr.as_ordered_terms()
        variables = expr.free_symbols
        variables = sorted(variables, key=lambda x: str(x))  

        result = []
        for term in terms:
            coef, mon = term.as_coeff_Mul()  
            if mon == 1:  
                tuple_mon = (0,) * len(variables)
            else:
                tuple_mon = tuple(mon.as_powers_dict().get(var, 0) for var in variables)

            if univariate: # If we are in the univariate case, we need to properly format the polynomial structure
                tuple_mon = (tuple_mon[0], 0)
            result.append((coef, tuple_mon))


        if len(result[0][1]) == 0:
            out_result = -1  # No solution found
        else:
            out_result = check_output(result)  # 0: solution not accepted as valid, 1: accepted

        print(out_result)
        print("------------------------")
        print(sp.expand(expand_expression(np.matrix(H), vd)))
    
# ---------------- Example of input ----------------
"""
f = [(1,(2, 4, 2)), (1,(2, 0, 2)), (16,(0, 4, 0)), (96,(0, 3, 1)), (216,(0, 2, 2)), (216,(0, 1, 3)), (81,(0, 0, 4))] 
f = [(2, (4, 0)), (2, (3, 1)), (-1, (2, 2)), (5, (0, 4)) ]
f = [(1,(10, 0)), (20,(9, 0)), (180,(8, 0)), (960,(7, 0)), (3360,(6, 0)), (8064,(5, 0)), (13440,(4, 0)),(15360,(3, 0)), (11520,(2, 0)),   (5120,(1, 0)), (1024,(0, 0))]
"""

if __name__ == "__main__":
    f = [(-1,(0, 0)), (1,(1, 0)), (1,(2, 0)), (1,(3, 0)), (1,(4, 0))] 
    main(f)
# -------------------------------------------------

