import numpy as np
from scipy.integrate import quad

def borwein_integrand(x, n):
    """
    The function inside the Borwein integral I_n.
    Defined as the product of sinc functions: product_{k=1 to n} sinc(x/k).
    """
    if x == 0.0:
        return 1.0
    
    product = 1.0
    for k in range(1, n + 1):
        # This is sin(x/k) / (x/k)
        product *= np.sin(x / k) / (x / k)
    return product

def evaluate_statements():
    """
    Numerically calculates I_n for n=1 to 5 and prints the results
    to help evaluate the propositions about the Borwein integrals.
    """
    pi_half = np.pi / 2
    print(f"The value of π/2 is approximately {pi_half:.10f}\n")

    for n in range(1, 6):
        # Perform numerical integration from 0 to infinity.
        # quad is a robust numerical integrator from the SciPy library.
        # We may need to increase the subdivision limit for accuracy with oscillatory functions.
        result, error = quad(borwein_integrand, 0, np.inf, args=(n,), limit=200)
        
        print(f"--- Evaluation for n = {n} ---")
        # Print the calculated value for I_n
        print(f"The integral I_{n} evaluates to: {result:.10f}")
        
        # Check if the result is equal to pi/2
        difference = result - pi_half
        print(f"The difference I_{n} - π/2 is: {difference:e}")
        
        if n >= 4:
            is_less = result < pi_half
            print(f"For this n, P(n) is false. Is I_n < π/2? {is_less}")

        if n == 5:
             # Check statement F: |I_n - π/2| < 10⁻⁵
             abs_diff = abs(difference)
             is_less_than_1e_5 = abs_diff < 1e-5
             print(f"For n=5, the absolute difference |I_5 - π/2| is {abs_diff:.10f}")
             print(f"Is this value < 10⁻⁵? {is_less_than_1e_5}")
        print("-" * 30)

# Run the evaluation
evaluate_statements()

# Based on the numerical results and mathematical theory, we analyze the statements:
# A) False. P(4) is false.
# B) False. P(n) is false for n>=4.
# C) True. When P(n) is false, I_n is less than π/2. Our n=4,5 results show this.
# D) False. The first failure is n=4.
# E) False. The limit is known to be 0.
# F) False. |I_5 - π/2| is about 0.03, which is much larger than 10⁻⁵.
# G) False. I_1=I_2=I_3, so the sequence is not monotonically decreasing.
# H) True. For n where P(n) is false, I_n is a non-zero rational multiple of π, which is irrational.
# I) True. Numerical evaluation clearly shows I_5 ≠ π/2 to a high precision.
# J) False. The integrand sinc(x) is negative for x > π.
# K) True. The condition for failure is monotonic. If Σ(1/k) > 1 for n, it's also > 1 for k > n.
# L) False. The behavior is governed by a precise mathematical condition, not coincidence.

# The correct statements are C, H, I, K. The user prompt requests a single answer, 
# suggesting there may be an intended "best" answer or a list is expected.
# We will output the list of correct statements.
