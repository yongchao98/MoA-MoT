import numpy as np

def solve():
    """
    This function calculates the value of l(14) based on the derived formula.
    """
    # Step 1: Define the constants from the problem.
    n = 13
    p = 14

    # Step 2: Calculate R(p), the radius term.
    # R(p) = sqrt(n) * (e^(p/sqrt(n)) - 1)
    sqrt_n = np.sqrt(n)
    R_p = sqrt_n * (np.exp(p / sqrt_n) - 1)

    # Step 3: Calculate A = 1 + R(p).
    # This is the base of the powers in the integral's final formula.
    A = 1 + R_p

    # Step 4: Calculate l(p) using the simplified closed-form expression.
    # l(p) = 2 * log(2) * log(A)
    log_2 = np.log(2)
    log_A = np.log(A)
    ell_p = 2 * log_2 * log_A

    # As requested, output the numbers in the final equation.
    # The final equation is l(p) = 2 * log(2) * log(1 + R(p))
    print(f"For n = {n} and p = {p}:")
    print(f"The radius term R({p}) is: {R_p}")
    print(f"The base term A = 1 + R({p}) is: {A}")
    print(f"The final result is calculated as: 2 * log(2) * log({A})")
    print(f"The value of l({p}) is: {ell_p}")
    
    # Final answer in the required format
    print(f"<<<{ell_p}>>>")

solve()