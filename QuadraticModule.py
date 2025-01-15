import sys, json, itertools, time, cvxpy as cp, numpy as np, sympy as sp, symengine as se
from sympy import symbols, expand
x1, x2, x3 = symbols('x1 x2 x3')

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
    dec_matrices = list()

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

        dec_matrices.append(desc_matrix)

    return list(no_rep_monom), dec_matrices

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
def prod_decomposition(Q_matrices):
    desc = list()
    for matrix in Q_matrices:
        U, s, Vt = np.linalg.svd(matrix) 
        s_sqrt = np.sqrt(s)
        S_sqrt = np.diag(s_sqrt)
        temp_H = np.dot(U, S_sqrt)
        desc.append(temp_H)
        
        """
        # Check that Q = H * H^T
        reconstructed_Q = np.dot(temp_H, temp_H.T)
        
        print("\nReconstructed Q:")
        print(reconstructed_Q)
        print("\nH:")
        print(temp_H)
        """
    
    return desc

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
        if element[1] not in f_monomials and element[0] >= tol:
            out_result = 0
            break

    return out_result

"""
Method to search for a decomposition in the quadratic module
"""
def QM(f, g, vd):
    Q = []
    constraints = []
    # Create the unknown matrices Qᵢ and impose the SDP constraints 
    for i in range(len(g)): 
        exec(f"Q{i} = cp.Variable((len(vd[i]), len(vd[i])), symmetric=True)")
        exec(f"Q.append(Q{i})")
        exec(f"constraints.append(Q{i} >> 0)")
    
    prods_vd = list()
    for vd_monomials in vd:
        # Matrix of monomials resulting from vd * vd'
        prod_vd = [[tuple(a[i] + c[i] for i in range(len(vd_monomials[0]))) for c in vd_monomials] for a in vd_monomials] 
        prods_vd.append(prod_vd)

    final_expression = list()
    for idx, prod in enumerate(prods_vd):
        fi = list() # We construct fi = vdᵢ' * Qi * vdᵢ * gᵢ
        l1, l2 = decomposition(prod)

        if idx == 0:
            for i in range(len(l2)):
                fi.append((cp.trace(Q[idx] @ l2[i]), l1[i]))
                
        else:
            for coef, mon in g[idx]:
                for i in range(len(l2)):
                    fi.append(( coef * cp.trace(Q[idx] @ l2[i]), tuple(sum(pair) for pair in zip(l1[i], mon)) ))

        final_expression.append(fi)

    diff_monomials = []
    for prod in final_expression:
        for coef, mon in prod: 
            if mon not in diff_monomials:
                diff_monomials.append(mon)
    
    # We store the coefficients of f in the same order as the different monomials that appear in the decomposition
    f_coef = []
    for monomial in diff_monomials:
        is_in_poly = False
        for coef, mon in f:
            if monomial == mon:
                f_coef.append(coef)
                is_in_poly = True
                continue
        if is_in_poly == False: f_coef.append(0)    

    
    # For each possible monomial, we go through the expanded products to take the corresponding coefficient and impose the constraints
    for idx, monomial in enumerate(diff_monomials):
        coef_sum = 0
        for prod in final_expression: # We are taking the coefficients of each monomial and accumulating them for each gᵢ
            for coef, mon in prod: 
                if mon == monomial:
                    coef_sum += coef
                    continue    
                    
        constraints += [ coef_sum == f_coef[idx]]
    
    prob = cp.Problem(cp.Minimize(0), constraints)
    
    try:
         prob.solve(solver='SCS')

         Q_values = []
         for i, Qi in enumerate(Q):
             if Qi.value is not None:
                 #print("\n")
                 #print(sdp(Qi.value))
                 #print(f"Q{i}:")
                 #print(Qi.value)
                 #print("...................................................................")
                 Q_values.append(Qi.value)
             else:
                 pass
   
    except cp.error.SolverError as e:
        #print("Error with the solver:", e)
        Q_values = [None for i in range(len(g))]


    return Q_values

if __name__ == "__main__":
    if len(sys.argv) > 3:
    	# We take the input from the C++ code and search for a decomposition with the specific degree combination received
        f = json.loads(sys.argv[1])
        g = json.loads(sys.argv[2])
        vd = json.loads(sys.argv[3])

        f = [ [element if isinstance(element, int) else tuple(element) for element in sublist] for sublist in f ]
       	univariate = all(len(monomial) == 2 and monomial[1] == 0 for _, monomial in f)
	
        Qs = QM(f, g, vd)
        
        # We check if we have obtained a solution and, if so, whether it is sufficiently accurate
        if any(q is None for q in Qs):
            print(-2) # If there is any error with the solver
            print("\n")
        else:
            Hs = prod_decomposition(Qs)
        
            final_expression = 0 
            # ---------------- It has to be changed by the user according to the input (it could be read from g) ----------------
            """
            g_expression = [1, x1*x3**2+1-x1**2-x2**2, -x1*x3**2+1] 
            g_expression = [1, x1, x2, 1-x1-x2]
            g_expression = [1, x1, x2, x3, 1-x1-x2-x3]
            g_expression = [1, x1-1/2, x2-1/2, x3-1/2, 1-x1*x2*x3]
            g_expression = [1, x1**3-x2**2, 1-x1]
            g_expression = [1, x1, x2, x3]
            """
            g_expression = [1, x1-1/2, x2-1/2, 1-x1*x2]
		
            for idx, H in enumerate(Hs):
                final_expression += expand_expression(Hs[idx], vd[idx]) * g_expression[idx]

            expr = sp.expand(final_expression)
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
                      
                if univariate: # If we are in the univariate case, we need to properly format the polynomial structur 
                   tuple_mon = (tuple_mon[0], 0)
                result.append((coef, tuple_mon))
       
            if len(result[0][1]) == 0:
                out_result = -1 # No solution found

            else:
                out_result = check_output(result) # 0: solution not accepted as valid, 1: accepted

            print(out_result)
            print("\n")
            #print(result)
            #print("...................................................................")

    else:
        print("Not enough arguments provided") 
