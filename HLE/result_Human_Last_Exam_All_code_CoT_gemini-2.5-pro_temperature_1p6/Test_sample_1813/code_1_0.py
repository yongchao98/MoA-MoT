import numpy as np
from fractions import Fraction

def get_continued_fraction(p, q):
    """Computes the continued fraction coefficients for p/q."""
    coeffs = []
    # We only care about the part after the integer part (the a_i for i>=1)
    if p >= q:
        p = p % q
    
    while q != 0:
        a = p // q
        coeffs.append(a)
        p, q = q, p % q
    
    # This algorithm produces [a0, a1, ...] for p/q >= 1
    # or [0, a1, a2, ...] for p/q < 1
    # We strip the first element which is the integer part
    # For 4/7 = [0; 1, 1, 3], it produces [0, 1, 1, 3]. We need [1, 1, 3]
    return coeffs[1:]

def get_cf_for_number(num, den):
    """Computes the continued fraction coefficients for num/den."""
    coeffs = []
    if den == 0:
        return []
    while den != 0:
        a = num // den
        coeffs.append(a)
        num, den = den, num % den
    return coeffs

def format_cf_equation(num, den, coeffs):
    """Formats the continued fraction as an equation string."""
    if not coeffs:
        return f"{num}/{den} = {num/den}"
    
    # Start with the rational number
    equation = f"{num}/{den} = "
    
    # Build the expression from right to left
    if len(coeffs) == 1:
        equation += f"{coeffs[0]}"
        return equation

    # For [c0, c1, ..., cn] -> c0 + 1 / (c1 + 1 / (... cn))
    # Start with the last term
    expr = str(coeffs[-1])
    # Iterate backwards from the second to last term
    for i in range(len(coeffs) - 2, 0, -1):
        expr = f"{coeffs[i]} + 1 / ({expr})"
    
    # Add the first term
    equation += f"{coeffs[0]} + 1 / {expr}"
    return equation


def compute_generalized_markov_number(p, q):
    """
    Computes the generalized Markov number m_p/q and its continued fraction.
    """
    # Step 1: Find continued fraction of p/q
    cf_coeffs = get_continued_fraction(p, q)

    # Step 2: Ensure even length
    if len(cf_coeffs) % 2 != 0:
        last_coeff = cf_coeffs.pop()
        if last_coeff > 1:
            cf_coeffs.append(last_coeff - 1)
            cf_coeffs.append(1)
        # If last_coeff is 1, it means we merge it with the previous, e.g., [..., a, 1] -> [..., a+1]
        # But our rule is a_n -> a_n-1, 1. So if a_n=1, this gives [..., 0, 1], which is invalid.
        # However, the standard algorithm ensures a_n > 1 if n > 0 for shortest representation.
        # e.g., 4/7 = [0; 1, 1, 3], not [0; 1, 1, 2, 1] initially. So a_n will be > 1.
    
    # Step 3: Construct the matrix
    R = np.array([[1, 1], [0, 1]], dtype=object)
    L = np.array([[1, 0], [1, 1]], dtype=object)
    
    M = np.identity(2, dtype=object)
    
    for i, a in enumerate(cf_coeffs):
        if i % 2 == 0:  # Odd-indexed coefficient a_1, a_3, ...
            term_matrix = np.linalg.matrix_power(R, a)
        else:  # Even-indexed coefficient a_2, a_4, ...
            term_matrix = np.linalg.matrix_power(L, a)
        M = np.dot(M, term_matrix)

    # Step 4: Calculate the generalized Markov number
    trace = np.trace(M)
    m_num = int(trace)
    m_den = 3
    
    m_fraction = Fraction(m_num, m_den)
    m_num = m_fraction.numerator
    m_den = m_fraction.denominator
    
    # Step 5: Compute the continued fraction of the result
    m_coeffs = get_cf_for_number(m_num, m_den)
    
    # Step 6: Format the output
    print(format_cf_equation(m_num, m_den, m_coeffs))


# Set the rational number p/q
p, q = 4, 7
compute_generalized_markov_number(p, q)
