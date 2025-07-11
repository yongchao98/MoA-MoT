import math

def find_c():
    """
    This function calculates the maximum real number c based on the derivation.
    """
    numerator = 2
    pi = math.pi
    denominator = pi**2
    c = numerator / denominator

    print("The problem is to find the maximum real number c such that for any positive integer n and any n real numbers x_1, ..., x_n, the inequality holds.")
    print("The inequality is: sum_{i=1 to n} sum_{j=1 to n} (n - |i-j|) * x_i * x_j >= c * sum_{i=1 to n} x_i^2")
    print("\nThis value of c is the infimum of the lowest eigenvalues of the associated matrices A_n for all n.")
    print("This value can be determined by solving the continuous analogue of the problem, which gives c = 2 / pi^2.")
    print("\nLet's calculate this value.")
    print(f"The equation for c is: c = {numerator} / ({pi:.6f})^2")
    print(f"Breaking it down: c = {numerator} / {denominator:.6f}")
    print(f"The final value for c is approximately: {c:.6f}")

find_c()