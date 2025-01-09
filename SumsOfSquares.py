import itertools, time, cvxpy as cp, numpy as np, sympy as sp, symengine as se
import mosek
from sympy import symbols, expand
x1, x2, x3, x4, x5, x6, x7 = symbols('x1 x2 x3 x4 x5 x6 x7')
import csv

def sdp(matrix, threshold=1e-4): 
    if matrix.shape[0] != matrix.shape[1]: # Verifiquem que la matriu és quadrada
        return False
    
    eigenvalues, _ = np.linalg.eig(matrix)
    
    # Filtrem els VAPS propers a 0 a partir de l'umbral
    eigenvalues[np.abs(eigenvalues) < threshold] = 0

    print(eigenvalues)
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
    
def prod_decomposition(matrix):
    U, s, Vt = np.linalg.svd(matrix) 
    s_sqrt = np.sqrt(s)
    S_sqrt = np.diag(s_sqrt)
    
    H = np.dot(U, S_sqrt)

    """
    # Verifiquem que Q = H * H^T
    reconstructed_Q = np.dot(H, H.T)
    
    print("\nReconstructed Q:")
    print(reconstructed_Q)
    print("\nH:")
    print(H)
    """
    
    return H
    
# Retorna tots els monomis possibles de n_variables variables i grau menor o igual que deg
def generate_monomials(n_variables, deg):
    monomials = []
    for comb in itertools.product(range(deg + 1), repeat=n_variables):
        if sum(comb) <= deg and comb not in monomials:
            if n_variables == 1:
                monomials.append((comb[0], 0))
            else:
                monomials.append(comb)
    return monomials
 
def SOS(f, vd, num_threads=1):
    Q = cp.Variable((len(vd),len(vd)), symmetric=True)
    constraints = [Q >> 0]

    # Matriu de monomis resultant de vd*vd'
    prod_vd = [[tuple(a[i] + c[i] for i in range(len(vd[0]))) for c in vd] for a in vd] 
    
    # Descomposem amb les matrius B_alpha
    l1, l2 = decomposition(prod_vd)
    
    # Igualem tr(Q, B_alpha)=g_alpha
    for monomial in l1:
        is_in = False
        for coef, monom in f:
            if monomial == monom:
                constraints += [cp.trace(Q @ l2[l1.index(monomial)]) == coef]
                is_in = True
            
                break
                
        if is_in == False:
            constraints += [cp.trace(Q @ l2[l1.index(monomial)]) == 0]
        
	
    prob = cp.Problem(cp.Minimize(0), constraints)
    #prob.solve(solver='MOSEK', verbose=False , mosek_params = {'MSK_IPAR_NUM_THREADS':  num_threads})
    prob.solve(solver='SCS')
    
    if Q.value is not None:
    	pass
    else:
    	print("No se ha encontrado solución.")
	
    return Q.value

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
        if element[1] not in f_monomials and abs(element[0]) >= 0.00001: 
            out_result = 0
            print("sota")
            print(abs(element[0]))
            break

    return out_result
    

f = [(1,(2, 4, 2)), (1,(2, 0, 2)), (16,(0, 4, 0)), (96,(0, 3, 1)), (216,(0, 2, 2)), (216,(0, 1, 3)), (81,(0, 0, 4))] 
vd = generate_monomials(3, 4)
Q = SOS(f, vd)
H = prod_decomposition(Q)
 
#----------------------------------------------------------------------------------------------------------------------- 
def main_func(deg_bound, num_threads):
    Q = SOS(f, vd, num_threads)
    H = prod_decomposition(Q)

    final_expression = expand_expression(np.matrix(H), vd)

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
        #print(resultado)
        
    #print("-------------------Resultat---------------------")
    #print(out_result)
    #print(expand_expression(np.matrix(H), vd))
    return temps

