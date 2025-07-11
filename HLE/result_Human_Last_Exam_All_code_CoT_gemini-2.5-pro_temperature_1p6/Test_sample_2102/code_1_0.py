import numpy as np
import sympy as sp
from collections import Counter

def solve():
    """
    Solves the problem by finding the smallest n where f(n) > 10,
    and then calculating n * ||W_n||_inf.
    """
    x = sp.Symbol('x')
    
    # Define the function symbolically using sympy
    # sp.elliptic_k(x) is the complete elliptic integral of the first kind K(x)
    g_x = (2 / sp.pi) * sp.elliptic_k(x) * sp.exp(x)

    n = 1
    while True:
        # We need a polynomial of degree n, so we need n+1 terms in the series
        # The series is computed around x=0
        series = g_x.series(x, 0, n + 1)

        # Extract coefficients from the symbolic series
        # The .removeO() method removes the O(x^n) term
        poly = sp.Poly(series.removeO(), x)
        
        # sympy provides coefficients from highest degree to lowest.
        # np.roots expects coefficients in the same order.
        # We convert sympy's arbitrary precision fractions to floats.
        coeffs = [sp.N(c) for c in poly.all_coeffs()]

        # Find the roots of the polynomial. These are the eigenvalues of S_n.
        roots = np.roots(coeffs)

        # f(n) is the sum of the absolute cubes of the eigenvalues.
        f_n = np.sum(np.abs(roots)**3)
        
        # Uncomment the following line to see the progress for each n
        # print(f"Trying n={n}, f(n) = {f_n:.4f}")

        if f_n > 10:
            # We found the smallest n that satisfies the condition.
            
            # Now, calculate ||W_n||_inf.
            # This requires checking for multiplicities of eigenvalues (roots).
            # We round the roots to a certain precision to check for equality numerically.
            rounded_roots = np.round(roots, decimals=8)
            
            # Use Counter on a hashable representation of the complex numbers (their bytes)
            # to find the multiplicity of each root.
            multiplicities = Counter(r.tobytes() for r in rounded_roots)
            
            # Create a map from the byte representation back to the complex number
            unique_roots_map = {r.tobytes(): r for r in rounded_roots}

            max_norm_val = 0.0
            for root_bytes, mult in multiplicities.items():
                # Get the complex root value
                root_val = unique_roots_map[root_bytes]
                
                # The norm of a Jordan/Weyr block is |lambda| if multiplicity is 1,
                # and |lambda| + 1 if multiplicity is > 1.
                norm = np.abs(root_val)
                if mult > 1:
                    norm += 1
                
                if norm > max_norm_val:
                    max_norm_val = norm
            
            infinity_norm_Wn = max_norm_val
            
            # Calculate the final result
            result = n * infinity_norm_Wn

            print(f"The smallest integer n for which f(n) > 10 is n = {n}.")
            print(f"For this n, the infinity norm ||W_n||_inf is {infinity_norm_Wn}.")
            print(f"The final result is n * ||W_n||_inf = {n} * {infinity_norm_Wn} = {result}.")
            
            return result

# Run the solver
final_answer = solve()
# The final answer is wrapped in <<<>>> as requested.
print(f"<<<{final_answer}>>>")