import sys
import json
import itertools, time, cvxpy as cp, numpy as np, sympy as sp, symengine as se
from sympy import symbols, expand
x1, x2, x3 = symbols('x1 x2 x3')
from contextlib import redirect_stdout

def sdp(matrix, threshold=1e-4): 
    if matrix.shape[0] != matrix.shape[1]: # Verifiquem que la matriu és quadrada
        return False
    
    eigenvalues, _ = np.linalg.eig(matrix)
    
    # Filtrem els VAPS propers a 0 a partir de l'umbral
    eigenvalues[np.abs(eigenvalues) < threshold] = 0

    #print(eigenvalues)
    return np.all(eigenvalues >= 0)

def decomposition(matrix):
    no_rep_monom = set()
    desc_matrices = list()

    # Guardem els monomis únics (ja que poden haver repeticions)
    for row in matrix:
        for element in row:
            no_rep_monom.add(tuple(element))

    # Per cada monomi creem la matriu corresponent amb 1's i 0's
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

def expand_expression(matrix, vd):
    matrix = np.array(matrix)
    
    # Definim les variables x1, x2, x3, ..., xn
    num_vars = len(vd[0])
    variables = sp.symbols(f'x1:{num_vars+1}')
    
    # Creem una llista amb les variables elevades a les potències corresponents
    variable_terms = []
    for exponents in vd:
        if all(e == 0 for e in exponents):
            term = 1  # El (0,0) és 1
        else:
            term = sp.Mul(*(variables[j]**exponents[j] for j in range(num_vars) if exponents[j] != 0))
        variable_terms.append(term)
    
    # Generem l'expressió
    expression = 0
    for col in range(matrix.shape[1]):
        term = 0
        for row in range(matrix.shape[0]):
            term += matrix[row, col] * variable_terms[row]  
        expression += term**2  

    return se.expand(expression)

def prod_decomposition(Q_matrices):
    desc = list()
    for matrix in Q_matrices:
        U, s, Vt = np.linalg.svd(matrix) 
        s_sqrt = np.sqrt(s)
        S_sqrt = np.diag(s_sqrt)
        temp_H = np.dot(U, S_sqrt)
        desc.append(temp_H)
        
        """
        # Verifiquem que Q = H * H^T
        reconstructed_Q = np.dot(temp_H, temp_H.T)
        
        print("\nReconstructed Q:")
        print(reconstructed_Q)
        print("\nH:")
        print(temp_H)
        """
    
    return desc

def check_output(polynomial): 
    f_monomials = [term[1] for term in f]
    out_result = 1
    for term in f:
        found = False
        for element in polynomial:
            if element[1] == term[1]:
                found = True
                if abs(element[0]-term[0]) >= 0.00001: 
                    out_result = 0
                break
        
        if out_result == 0:
            break

        if found == False:
            out_result = 0
            break
    
    for element in polynomial:
        if element[1] not in f_monomials and element[0] >= 0.00001:
            out_result = 0
            break

    return out_result

def QM_v(f, g, vd):
    Q = []
    constraints = []
    # Creem les matrius incògnita Q_i i imposem les restriccions de SDP 
    for i in range(len(g)): 
        exec(f"Q{i} = cp.Variable((len(vd[i]), len(vd[i])), symmetric=True)")
        exec(f"Q.append(Q{i})")
        exec(f"constraints.append(Q{i} >> 0)")
    
    prods_vd = list()
    for vd_monomials in vd:
        # Matriu de monomis resultant de vd*vd'
        prod_vd = [[tuple(a[i] + c[i] for i in range(len(vd_monomials[0]))) for c in vd_monomials] for a in vd_monomials] 
        prods_vd.append(prod_vd)

    final_expression = list()
    for idx, prod in enumerate(prods_vd):
        fi = list() # Ens muntem fi = vd_i'*Qi*vd_i*gi
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
    
    # Guardem els coeficients de f en el mateix ordre dels diferents monomis que apareixen en la descomposició
    f_coef = []
    for monomial in diff_monomials:
        is_in_poly = False
        for coef, mon in f:
            if monomial == mon:
                f_coef.append(coef)
                is_in_poly = True
                continue
        if is_in_poly == False: f_coef.append(0)    

    
    # Per cada possible monomi recorrem els productes desenvolupats per agafar el coeficient corresponent i imposar les restriccions
    for idx, monomial in enumerate(diff_monomials):
        coef_sum = 0
        for prod in final_expression: # Estem agafant els coeficients de cada monomi i sumant-los per a cada gi
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
        #print("Error con el solucionador:", e)
        Q_values = [None for i in range(len(g))]


    return Q_values

if __name__ == "__main__":
    if len(sys.argv) > 3:
        f = json.loads(sys.argv[1])
        g = json.loads(sys.argv[2])
        vd = json.loads(sys.argv[3])

        f = [ [elemento if isinstance(elemento, int) else tuple(elemento) for elemento in sublista] for sublista in f ]

        Qs = QM_v(f, g, vd)
        if any(q is None for q in Qs):
            print(-2)
            print("\n")
        else:
            Hs = prod_decomposition(Qs)
        
            final_expression = 0 
            #g_expression = [1, x1*x3**2+1-x1**2-x2**2, -x1*x3**2+1] 
            #g_expression = [1, x1, x2, 1-x1-x2]
            g_expression = [1, x1-1/2, x2-1/2, 1-x1*x2]
            #g_expression = [1, x1, x2, x3, 1-x1-x2-x3]
            #g_expression = [1, x1-1/2, x2-1/2, x3-1/2, 1-x1*x2*x3]
            #g_expression = [1, x1**3-x2**2, 1-x1]
            #g_expression = [1, x1, x2, x3]

            for idx, H in enumerate(Hs):
                final_expression += expand_expression(Hs[idx], vd[idx]) * g_expression[idx]

            expr = sp.expand(final_expression)

            terms = expr.as_ordered_terms()

            variables = expr.free_symbols
            variables = sorted(variables, key=lambda x: str(x))  

            resultado = []
            for term in terms:
                coeficiente, monomio = term.as_coeff_Mul()  
                if monomio == 1:  
                   tupla_monomio = (0,) * len(variables)
                else:
                   tupla_monomio = tuple(monomio.as_powers_dict().get(var, 0) for var in variables)
                
                resultado.append((coeficiente, tupla_monomio))

       
            if len(resultado[0][1]) == 0:
                out_result = -1 # No té solució

            else:
                out_result = check_output(resultado) # 0 solució no admesa com a vàlida, 1 admesa

            print(out_result)
            print("\n")
            #print(resultado)
            #print("...................................................................")

    else:
        print("Not enough arguments provided") 
