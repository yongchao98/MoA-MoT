import math
import numpy as np

def get_continued_fraction(p, q):
    """Computes the continued fraction coefficients of p/q."""
    coeffs = []
    while q != 0:
        a = p // q
        coeffs.append(a)
        p, q = q, p % q
    return coeffs

def compute_convergents_matrix(coeffs):
    """
    Computes the matrix M = M_1 * M_2 * ... * M_k
    where M_i = [[coeffs[i-1], 1], [1, 0]].
    The resulting matrix is [[P_k, P_{k-1}], [Q_k, Q_{k-1}]].
    """
    if not coeffs:
        # Identity matrix for an empty list of coefficients
        return np.array([[1, 0], [0, 1]])

    # Initialize with the first matrix
    m_first = np.array([[coeffs[0], 1], [1, 0]], dtype=object)
    m_total = m_first

    # Multiply by subsequent matrices
    for i in range(1, len(coeffs)):
        m_i = np.array([[coeffs[i], 1], [1, 0]], dtype=object)
        m_total = np.dot(m_total, m_i)

    return m_total

def main():
    """
    Main function to compute the continued fraction associated with m_{4/7}.
    """
    p, q = 4, 7

    # 1. Get the CF of the rational number
    cf_coeffs = get_continued_fraction(p, q)
    print(f"The continued fraction of {p}/{q} is: {cf_coeffs}")

    # The associated CF has a period from the original CF, excluding the integer part.
    integer_part = cf_coeffs[0]
    periodic_part = cf_coeffs[1:]
    print(f"The associated continued fraction is [0; overline{periodic_part}]")

    # 2. Compute the convergents matrix for the periodic part
    k = len(periodic_part)
    matrix = compute_convergents_matrix(periodic_part)
    
    Pk = matrix[0, 0]
    Pk_minus_1 = matrix[0, 1]
    Qk = matrix[1, 0]
    Qk_minus_1 = matrix[1, 1]
    
    # 3. Derive the quadratic equation for x = [0; periodic_part]
    # The equation is P_{k-1}*x^2 + (P_k - Q_{k-1})*x - Q_k = 0
    A = Pk_minus_1
    B = Pk - Qk_minus_1
    C = -Qk
    
    print(f"\nThe value x of the continued fraction is a root of the quadratic equation Ax^2 + Bx + C = 0.")
    print(f"Substituting the computed values:")
    print(f"A = P_{k-1} = {Pk_minus_1}")
    print(f"B = P_k - Q_{k-1} = {Pk} - {Qk_minus_1} = {B}")
    print(f"C = -Q_k = -{Qk}")
    
    # Normalize the equation if possible
    common_divisor = math.gcd(math.gcd(A, B), C)
    A //= common_divisor
    B //= common_divisor
    C //= common_divisor

    print(f"\nThe final equation is {A}*x^2 + {B}*x + {C} = 0")

    # 4. Solve the quadratic equation
    discriminant = B**2 - 4 * A * C
    # Since x must be positive, we take the positive root.
    x = (-B + math.sqrt(discriminant)) / (2 * A)

    print(f"\nThe value of the continued fraction is the positive root of this equation, which is: {x}")

if __name__ == "__main__":
    main()
