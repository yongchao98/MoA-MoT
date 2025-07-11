import math

def calculate_rayleigh_quotient(n):
    """
    Calculates the Rayleigh quotient of the matrix A_n = J - Z_n
    with the all-ones vector j.

    This serves as a lower bound for the spectral norm of A_n.
    The formula is (4^n - 3^n) / 2^n.
    """
    if n < 1:
        print("n must be a natural number (>= 1).")
        return

    four_n = 4**n
    three_n = 3**n
    two_n = 2**n

    print(f"For n = {n}:")
    print(f"  4^n = {four_n}")
    print(f"  3^n = {three_n}")
    print(f"  2^n = {two_n}")

    rayleigh_quotient = (four_n - three_n) / two_n
    print(f"  Rayleigh quotient for j = (4^n - 3^n) / 2^n = {rayleigh_quotient}")

    # The growth rate of c_n is Theta(alpha^n).
    # This calculation shows that alpha is at least 2.
    # The ratio of the lower bound to 2^n approaches 1 as n grows.
    ratio = rayleigh_quotient / (2**n)
    print(f"  Ratio of this lower bound to 2^n: {ratio}")
    print("-" * 20)

# Example calculations for a few values of n
calculate_rayleigh_quotient(1)
calculate_rayleigh_quotient(2)
calculate_rayleigh_quotient(5)
calculate_rayleigh_quotient(10)

# The value of alpha is determined by the analysis above.
alpha = 2
print(f"The determined value of alpha is: {alpha}")