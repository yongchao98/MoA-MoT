import numpy as np

def solve_phi_n():
    """
    This function solves the problem for a specific value of n.
    The problem asks to calculate phi(n) for n >= 5. As n is not specified,
    we choose the smallest possible integer value, n = 5.

    The derivation shows that phi(n) = exp(2n - 4 + 2/n).
    """
    n = 5

    # The final equation for the exponent of phi(n) is 2*n - 4 + 2/n
    val_2n = 2 * n
    val_neg_4 = -4
    val_2_over_n = 2 / n
    
    exponent = val_2n + val_neg_4 + val_2_over_n
    
    # Per the instructions, we output the numbers in the final equation for the exponent.
    print(f"The numbers in the final exponent equation for n={n} are: {val_2n}, {val_neg_4}, and {val_2_over_n}")
    
    # The value of the exponent is
    print(f"The exponent is: {exponent}")
    
    # The final result is exp(exponent)
    result = np.exp(exponent)
    print(f"The final value of phi({n}) is: {result}")

solve_phi_n()
