import sys

def solve_polynomial_coeffs(p, k, n):
    """
    Computes the number of coefficients not divisible by p^k in the given iterated polynomial.

    Args:
        p: An odd prime number (p >= 3).
        k: An integer (k >= 1).
        n: An integer (n >= 1).
    """

    # Based on the mathematical analysis, the problem can be split into two cases.
    # Case 1: k = 1.
    # The number of coefficients not divisible by p is 2.
    
    # Case 2: k >= 2.
    # The number of coefficients not divisible by p^k can be shown to be k+1,
    # independent of p and n.
    
    # Combining these cases, for k=1, k+1 = 2. So the formula k+1 works for all k >= 1.
    
    # The number of coefficients is k + 1.
    answer = k + 1

    print(f"For the given inputs p = {p}, k = {k}, n = {n}:")
    print("The number of coefficients in the final polynomial that are not divisible by p^k is given by the formula k + 1.")
    print(f"Result: {k} + 1 = {answer}")

if __name__ == '__main__':
    # You can change these values to test with other inputs.
    # Per the problem, p must be an odd prime, and k, n must be >= 1.
    p_val = 5
    k_val = 4
    n_val = 2
    
    if len(sys.argv) == 4:
        try:
            p_val = int(sys.argv[1])
            k_val = int(sys.argv[2])
            n_val = int(sys.argv[3])
        except ValueError:
            print("Usage: python your_script.py p k n (with integer values)")
            sys.exit(1)

    # Note: No validation for p being an odd prime is included, as the formula is independent of p.
    if k_val < 1 or n_val < 1:
        print("Error: k and n must be integers greater than or equal to 1.")
    else:
        solve_polynomial_coeffs(p_val, k_val, n_val)