import sympy
from sympy import elliptic_k, exp, pi, series
import numpy as np

def solve():
    """
    Finds the smallest n where f(n) > 10 and calculates n * ||W_n||_inf.
    """
    # Define the symbolic function g(x)
    x = sympy.Symbol('x')
    g = (2 / pi) * elliptic_k(x) * exp(x)

    n = 0
    while True:
        n += 1
        # Get the Taylor polynomial of degree n. This requires n+1 terms (from x^0 to x^n).
        poly_sympy = series(g, x, 0, n + 1).removeO()

        # Extract coefficients. The list will be [c_0, c_1, ..., c_n].
        coeffs = [float(poly_sympy.coeff(x, k)) for k in range(n + 1)]

        # numpy.roots expects coefficients from the highest power to the lowest.
        coeffs_for_roots = list(reversed(coeffs))

        # Find the roots of the polynomial. These are the eigenvalues of W_n.
        roots = np.roots(coeffs_for_roots)

        # Calculate f(n), the sum of the absolute cubes of the eigenvalues.
        f_n = sum(abs(r)**3 for r in roots)

        if f_n > 10:
            # Smallest n has been found.
            
            # The infinity norm of the Weyr form W_n is the maximum absolute eigenvalue,
            # assuming all eigenvalues are distinct (which they are in this case).
            # If there were repeated eigenvalues, the norm would be max(|lambda|) + 1.
            infinity_norm_Wn = max(abs(r) for r in roots)
            
            # Calculate the final result.
            result = n * infinity_norm_Wn
            
            # Print the final equation as requested.
            print(f"The smallest n where f(n) > 10 is n = {n}.")
            print(f"For this n, the infinity norm ||W_n||_inf is {infinity_norm_Wn}.")
            print(f"The final calculation is:")
            print(f"{n} * {infinity_norm_Wn} = {result}")
            
            # Return the final numerical answer
            return result

# Run the solver and capture the final numerical answer for the specified format.
final_answer = solve()
print(f"<<<{final_answer}>>>")
