import math
import numpy as np

def combinations(n, k):
    """Computes the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def calculate_Sn(n):
    """Calculates the sum S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k)"""
    total_sum = 0
    for k in range(n + 1):
        term = ( (2 * k + 1)**5 * 
                 combinations(2 * k, k) * 
                 combinations(2 * n - 2 * k, n - k) )
        total_sum += term
    return total_sum

def find_closed_form():
    """Finds the polynomial P(n) such that S_n = P(n) * 4^n"""
    # Calculate S_n for n = 0, 1, 2, 3
    S_values = [calculate_Sn(n) for n in range(4)]
    
    # Calculate P(n) = S_n / 4^n
    P_values = [S_values[n] // (4**n) for n in range(4)]
    
    # We assume P(n) = An^3 + Bn^2 + Cn + D.
    # We set up a system of linear equations M * x = P_values
    # where x = [D, C, B, A]
    M = np.array([
        [1, 0, 0, 0],  # n=0
        [1, 1, 1, 1],  # n=1
        [1, 2, 4, 8],  # n=2
        [1, 3, 9, 27]  # n=3
    ])
    
    # Solve for the coefficients [D, C, B, A]
    coeffs = np.linalg.solve(M, P_values)
    D, C, B, A = [int(round(c)) for c in coeffs]
    
    print("The closed form for the sum is of the form: (A*n^3 + B*n^2 + C*n + D) * 4^n")
    print(f"The calculated coefficients are: A={A}, B={B}, C={C}, D={D}")
    print("\nSo the closed form is:")
    
    # Format the output nicely
    def format_term(coeff, power):
        if coeff == 0:
            return ""
        
        sign = "+ " if coeff > 0 else "- "
        abs_coeff = abs(coeff)
        
        if power == 0:
            return f"{sign}{abs_coeff}"
        
        coeff_str = str(abs_coeff) if abs_coeff != 1 else ""
        power_str = f"n^{power}" if power > 1 else "n"
            
        return f"{sign}{coeff_str}{power_str}"

    poly_str = f"({A}n^3 {format_term(B, 2)}{format_term(C, 1)}{format_term(D, 0)}) * 4^n"
    poly_str = poly_str.replace("(+", "(").replace(" -", "- ").replace("+", "+ ")
    print(f"S_n = {poly_str}\n")
    
    # Verification
    print("Verification for n=0, 1, 2, 3:")
    for n in range(4):
        s_val = S_values[n]
        formula_val = (A * n**3 + B * n**2 + C * n + D) * (4**n)
        print(f"For n={n}:")
        print(f"  Direct Sum S_{n} = {s_val}")
        print(f"  Formula Result = {formula_val}")
        print(f"  Match: {s_val == formula_val}")

find_closed_form()
