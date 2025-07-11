import math

def calculate_l(n: int):
    """
    Calculates the value of l(n) for n >= 5.

    The function follows these analytical steps:
    1. Define matrices M_n and P_n.
    2. Compute Q = f^(3)(P_n), which is the inverse of P_n. Q is a tridiagonal
       matrix with 2 on the diagonal and 1 on the off-diagonals.
    3. Compute X = f_M^(2)(Q), the projection of Q onto the tangent space at M_n.
       X = Q - M_n * D, where D = mdiag(M_n^T * Q).
    4. Compute l(n) = f^(1)(X), the sum of the first and last rows of X.
       l(n) = sum_rows(Q) - sum_rows(M_n * D).
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("n must be an integer greater than or equal to 5.")

    # Let's compute l(n) using the final derived formula.
    # sum_rows(Q) = 6.
    # The term sum_rows(M_n * D) is denoted as B.
    # The final simplified formula for l(n) is:
    # l(n) = 2 + 2/n^2 - (4n-2)/n^2 * sqrt(n^2 - n + 1)

    term1 = 2
    term2 = 2 / (n**2)
    sqrt_term = math.sqrt(n**2 - n + 1)
    term3 = ((4 * n - 2) / (n**2)) * sqrt_term

    l_n = term1 + term2 - term3
    
    print("This script calculates l(n) based on its derived analytical formula.")
    print(f"For n = {n}:")
    print(f"l({n}) = 2 + 2/({n}^2) - (4*{n}-2)/({n}^2) * sqrt({n}^2 - {n} + 1)")
    print(f"l({n}) = {term1} + {term2} - {term3}")
    print(f"The final exact value is: {l_n}")

# The user can specify any n >= 5. As an example, we use n=5.
n_value = 5
calculate_l(n_value)
