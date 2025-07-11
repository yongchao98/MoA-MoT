import numpy as np
from scipy.integrate import quad
import math

def solve_borwein_problem():
    """
    This script evaluates the Borwein integral I_n for several values of n
    and provides the necessary data to answer the problem questions.

    Definition:
    I_n = integral from 0 to infinity of Product_{k=1 to n} [sin(x/k) / (x/k)] dx

    The proposition P(n) is that I_n = pi/2.
    This holds if and only if the n-th Harmonic number H_n <= 2.
    H_n = 1 + 1/2 + 1/3 + ... + 1/n.
    """

    def sinc(x):
        # Use np.sinc which is defined as sin(pi*x)/(pi*x)
        # To match sin(x)/x, we need to use sin(t)/t = sinc(t/pi)
        # Or, we can just define it directly to avoid confusion.
        if x == 0:
            return 1.0
        return np.sin(x) / x

    def integrand(x, n):
        """The function inside the Borwein integral."""
        if x == 0:
            return 1.0
        
        prod = 1.0
        for k in range(1, n + 1):
            prod *= sinc(x / k)
        return prod

    pi_half = np.pi / 2
    h_n = 0.0

    print("Evaluating Borwein Integrals I_n and checking condition H_n <= 2\n")
    print(f"Value of π/2 = {pi_half:.10f}")
    print("-" * 70)
    print(f"{'n':<3} | {'H_n':<12} | {'H_n <= 2?':<11} | {'I_n (numerical)':<18} | {'I_n - π/2':<18}")
    print("-" * 70)

    # We evaluate for n from 1 up to 7 to see the trend
    for n in range(1, 8):
        # Calculate Harmonic Number H_n
        h_n += 1.0 / n
        
        # Numerically calculate the integral I_n
        # We integrate over a large interval, as quad is efficient.
        # np.inf can be used for the upper limit.
        result, error = quad(integrand, 0, np.inf, args=(n,))

        # The phrase "output each number in the final equation" can be interpreted
        # as showing the components that lead to the result. We show n, H_n, and I_n.
        
        condition_met = "True" if h_n <= 2 else "False"
        
        # For the final output, let's print the numbers involved for each n.
        print(f"{n:<3} | {h_n:<12.5f} | {condition_met:<11} | {result:<18.10f} | {result - pi_half:<18.10f}")

solve_borwein_problem()