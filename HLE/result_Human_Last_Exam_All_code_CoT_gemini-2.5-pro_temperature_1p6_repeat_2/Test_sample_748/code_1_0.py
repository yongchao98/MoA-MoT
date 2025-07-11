import numpy as np
import sympy

def poly_to_string(p):
    """Converts a sympy Poly object to a human-readable string representation."""
    if p.is_number:
        return str(p.as_expr())
        
    x = p.gen
    coeffs = p.all_coeffs()
    degree = p.degree()
    
    terms = []
    for i, coeff in enumerate(coeffs):
        power = degree - i
        
        if coeff == 0:
            continue
            
        # Coefficient part
        if power == degree: # leading term
             if coeff == 1:
                 term_str = ""
             elif coeff == -1:
                 term_str = "-"
             else:
                 term_str = str(coeff)
        else:
            if coeff > 0:
                sign = " + "
            else:
                sign = " - "
            
            abs_coeff = abs(coeff)
            if abs_coeff == 1 and power != 0:
                coeff_str = ""
            else:
                coeff_str = str(abs_coeff)
            term_str = sign + coeff_str

        # Variable part
        if power == 1:
            term_str += str(x)
        elif power > 1:
            term_str += f"{x}**{power}"
        
        terms.append(term_str)
        
    return "".join(terms).lstrip(" + ")


def demonstrate_discontinuity():
    """
    Demonstrates the discontinuity of the minimal polynomial map at a derogatory matrix.
    """
    X = sympy.Symbol('X')
    
    # Define a 2x2 derogatory matrix M (the identity matrix)
    M = np.array([[1.0, 0.0], [0.0, 1.0]])
    M_sympy = sympy.Matrix(M)
    # The minimal polynomial of M = I is X - 1.
    pi_M_poly = M_sympy.minpoly(X)

    # Consider a sequence of matrices M_k converging to M as k -> infinity.
    # Let M_k = [[1, 1/k], [0, 1]]. M_k is non-derogatory for k > 0.
    # The minimal polynomial of M_k is its characteristic polynomial, (X-1)**2.
    # The limit of these polynomials is therefore (X-1)**2.
    
    # We can compute the limit polynomial by taking any matrix from the sequence (e.g., k=1)
    Mk_for_limit_calc = sympy.Matrix([[1, 1], [0, 1]]) 
    limit_pi_Mk_poly = Mk_for_limit_calc.minpoly(X)
    
    # Display the results
    print("--- Demonstration of Discontinuity ---")
    print(f"Let M be the derogatory matrix:\n{M}\n")
    
    pi_M_str = poly_to_string(pi_M_poly)
    print(f"The minimal polynomial of M is: pi_M(X) = {pi_M_str}")
    coeffs_pi_M = pi_M_poly.all_coeffs()
    print("Equation form: {}*X + ({}) = 0".format(coeffs_pi_M[0], coeffs_pi_M[1]))
    print("-" * 35)

    print("Now consider a sequence M_k -> M, for example, M_k = [[1, 1/k], [0, 1]].")
    print("For any k>0, M_k is non-derogatory.")
    print("The minimal polynomial of M_k is always (X-1)**2.")
    
    limit_pi_Mk_str = poly_to_string(limit_pi_Mk_poly)
    print(f"The limit of the minimal polynomials is: limit_pi(X) = {limit_pi_Mk_str}")
    coeffs_limit_pi = limit_pi_Mk_poly.all_coeffs()
    print("Equation form: {}*X**2 + ({})*X + {} = 0".format(coeffs_limit_pi[0], coeffs_limit_pi[1], coeffs_limit_pi[2]))

    print("-" * 35)
    print("Comparing the polynomials:")
    print(f"pi_M(X)      = {pi_M_str}")
    print(f"limit_pi(X)  = {limit_pi_Mk_str}")
    
    if pi_M_str != limit_pi_Mk_str:
        print("\nThe two polynomials are not equal. The map is discontinuous at M.")
    else:
        print("\nThis case is not expected.")

if __name__ == '__main__':
    demonstrate_discontinuity()