import numpy as np
from scipy.special import comb

def get_ehrhart_polynomial(d):
    """
    Computes the Ehrhart polynomial p(z) = (z+1)^d for a given dimension d.
    Returns the coefficients of the polynomial in decreasing order of power.
    """
    # The coefficients of (z+1)^d are given by the binomial theorem
    # (z+1)^d = sum_{k=0 to d} C(d,k) * z^k
    # We want them in order of z^d, z^{d-1}, ..., z^0
    coeffs = [comb(d, k, exact=True) for k in range(d, -1, -1)]
    return np.array(coeffs)

def solve_and_print(d):
    """
    Calculates the Ehrhart polynomial for dimension d, prints it,
    and finds its roots.
    """
    if d < 2:
        print("Dimension d must be at least 2.")
        return

    print(f"For dimension d = {d}:")
    
    # The Ehrhart polynomial is p(z) = (z+1)^d
    coeffs = get_ehrhart_polynomial(d)
    
    # Create the polynomial object
    poly = np.poly1d(coeffs)
    
    # Print the polynomial equation
    equation_str = []
    for i, c in enumerate(coeffs):
        power = d - i
        if c == 0:
            continue
        
        term = ""
        if power > 0:
            if c != 1:
                term += f"{c}*z"
            else:
                term += "z"
            if power > 1:
                term += f"^{power}"
        else:
            term += str(c)
        
        equation_str.append(term)
        
    print("The Ehrhart polynomial equation is:")
    print(" + ".join(equation_str) + " = 0")
    
    # Find and print the roots
    roots = poly.roots
    print("\nThe roots of the polynomial are:")
    print(roots)
    
    # Check the real part of the roots
    real_parts = np.real(roots)
    if np.all(np.isclose(real_parts, -1)):
        print("\nConclusion: Every root has a real part of -1.")
    else:
        print("\nConclusion: Not every root has a real part of -1.")

# Example for d=4
solve_and_print(d=4)