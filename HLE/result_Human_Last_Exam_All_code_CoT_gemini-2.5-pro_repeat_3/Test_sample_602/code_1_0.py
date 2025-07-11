import math
import sys

def solve(n):
    """
    Calculates the exact value of l(n) for n >= 5.

    The formula for l(n) is derived as:
    l(n) = f1(A) - f1(M_n * D')
    where:
    - f1(A) = 6
    - f1(M_n * D') = 2 * S
    - S = 2*a**2 + (2*n-1)*a*b + (2*n-3)*b**2
    - a = sqrt(1 - (n-1)/n**2)
    - b = 1/n
    """
    if not isinstance(n, int) or n < 5:
        print("Error: Input must be an integer n >= 5.", file=sys.stderr)
        return

    # Step 1: Define a and b from the matrix M_(n)
    a = math.sqrt(1 - (n - 1) / (n**2))
    b = 1 / n

    # Step 2: Calculate the sum S = S_1 = S_n
    # S = 2*a**2 + (2*n-1)*a*b + (2*n-3)*b**2
    term_a_sq = 2 * a**2
    term_ab = (2 * n - 1) * a * b
    term_b_sq = (2 * n - 3) * b**2
    S = term_a_sq + term_ab + term_b_sq
    
    # The term to be subtracted is S_1 + S_n = 2*S
    subtracted_term = 2 * S

    # Step 3: Calculate f^(1)(A), which is a constant
    f1_A = 6.0

    # Step 4: Calculate the final value of l(n)
    l_n = f1_A - subtracted_term
    
    # Final equation structure: l(n) = f1(A) - (S_1 + S_n)
    print(f"For n = {n}:")
    print(f"Term from A, f1(A) = {f1_A}")
    print(f"Term from projection, f1(M*D') = {subtracted_term}")
    print(f"Final value, l({n}) = {f1_A} - {subtracted_term} = {l_n}")


if __name__ == '__main__':
    # Example usage:
    # You can change this value to any integer >= 5
    n_value = 5
    solve(n_value)
