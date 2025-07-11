import numpy as np
from scipy.integrate import quad

def integrand(x, n):
    """
    The integrand of the Borwein integral I_n.
    The product of sinc functions: Π [k=1 to n] sin(x/k)/(x/k)
    """
    if x == 0.0:
        # The limit of sin(y)/y as y -> 0 is 1. The product is 1.
        return 1.0
    
    prod = 1.0
    for k in range(1, n + 1):
        y = x / k
        prod *= np.sin(y) / y
    return prod

def calculate_In(n):
    """
    Calculates the Borwein integral I_n using numerical quadrature.
    """
    # We integrate from 0 to infinity. The quad function is suitable for this.
    val, err = quad(integrand, 0, np.inf, args=(n,))
    return val

# The problem asks us to evaluate several propositions P(n), where P(n) is I_n = π/2.
# A key feature of these integrals is that P(n) is true for n=1,...,7 but false for n=8.
# Statement (C) claims that if P(n) is false, then I_n < π/2.
# We will test this by calculating I_n for the first case where P(n) is false, which is n=8.

n_test = 8
i_n_value = calculate_In(n_test)
pi_half = np.pi / 2

# Per the instructions, we output the numbers for the final equation,
# which demonstrates the relationship between I_8 and π/2.

print(f"We are testing statement C: If P(n) is false, then I_n < pi/2.")
print(f"The first value for which P(n) is false is n = {n_test}.")
print(f"We now form the equation showing the comparison of I_{n_test} and pi/2.")

print("\n-------------------------------------------")
print(f"CALCULATED VALUES:")
print(f"Value of I_{n_test} = {i_n_value:.16f}")
print(f"Value of pi/2 = {pi_half:.16f}")
print("-------------------------------------------")


# Final output showing the numbers in the inequality.
print("\nFINAL EQUATION:")
if i_n_value < pi_half:
    print(f"{i_n_value:.16f} < {pi_half:.16f}")
    print("\nThe calculation confirms that I_8 is less than pi/2, supporting statement C.")
else:
    print(f"{i_n_value:.16f} is NOT less than {pi_half:.16f}")
    print("\nThe calculation does not support statement C.")
