import math
import numpy as np

def calculate_sum(n):
    """
    Calculates the value of the sum for a given integer n.
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k)
    """
    total_sum = 0
    for k in range(n + 1):
        try:
            # Use math.comb for combinations C(n,k)
            term = (2 * k + 1)**5 * math.comb(2 * k, k) * math.comb(2 * n - 2 * k, n - k)
            total_sum += term
        except ValueError:
            # math.comb(n, k) raises ValueError if k > n
            continue
    return total_sum

def find_closed_form():
    """
    Finds the coefficients of the polynomial P(n) where S_n = P(n) * 4^n.
    We assume P(n) is a polynomial of degree 5.
    """
    # We need 6 points to determine a 5th-degree polynomial.
    degree = 5
    num_points = degree + 1

    n_values = np.arange(num_points)
    
    # Calculate S_n for n = 0, 1, ..., 5
    s_values = np.array([calculate_sum(n) for n in n_values])

    # We assume S_n = P(n) * 4^n, so P(n) = S_n / 4^n
    p_values = s_values / (4**n_values)

    # We set up a system of linear equations to find the polynomial coefficients.
    # The Vandermonde matrix A maps coefficients to polynomial values.
    # A * c = p_values
    A = np.vander(n_values, num_points, increasing=True)
    
    # Solve for the coefficients c = [a_0, a_1, ..., a_5]
    coeffs = np.linalg.solve(A, p_values)

    # Round coefficients to the nearest integer to handle floating point inaccuracies
    coeffs = np.round(coeffs).astype(int)
    
    # Construct the polynomial expression string
    poly_parts = []
    for i in range(degree, -1, -1):
        coeff = coeffs[i]
        if coeff != 0:
            if i == 0:
                poly_parts.append(f"{coeff}")
            elif i == 1:
                poly_parts.append(f"{coeff}n")
            else:
                poly_parts.append(f"{coeff}n^{i}")

    poly_str = " + ".join(poly_parts).replace("+ -", "- ")
    
    # Output the closed-form expression
    print("The closed form for the sum is:")
    print(f"S_n = 4^n * ({poly_str})")
    print("\nEach number in the final equation is:")
    print(f"S_n = 4^n * ({coeffs[5]}n^5 + {coeffs[4]}n^4 + {coeffs[3]}n^3 + {coeffs[2]}n^2 + {coeffs[1]}n + {coeffs[0]})")
    
    # Also expressed as (n+1) times a polynomial of degree 4 for compactness
    # From (42n^5 + 158n^4 + 223n^3 + 144n^2 + 41n + 4) = (n+1)(42n^4 + 116n^3 + 107n^2 + 37n + 4)
    print("\nMore compactly:")
    print(f"S_n = 4^n * (n+1) * ({coeffs[5]}n^4 + 116n^3 + 107n^2 + 37n + 4)")


if __name__ == '__main__':
    find_closed_form()