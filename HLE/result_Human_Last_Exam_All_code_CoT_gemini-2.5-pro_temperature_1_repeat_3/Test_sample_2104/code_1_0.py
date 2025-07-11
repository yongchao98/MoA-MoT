import math
from scipy.special import beta

def solve_problem():
    """
    Solves the problem by finding n1, n2, and then calculating T((n1-1)/n2).
    """
    # Step 1 & 2: Determine n1, n2, and alpha.
    # From the analysis, the first integer n1=2 is found by applying the standard formula for the order.
    # The second integer n2=4 is found by considering the known case of order reduction.
    n1 = 2
    n2 = 4
    
    # The argument for the function T is (n1-1)/n2
    alpha_num = n1 - 1
    alpha_den = n2
    alpha = alpha_num / alpha_den

    # Step 3: Identify T(alpha) and calculate its value.
    # The function T(alpha) is interpreted as the Beta function B(alpha, 1/2).
    # We need to calculate B(1/4, 1/2).
    
    # As per the instructions, we output the numbers in the final equation.
    print(f"The equation to be solved is T(({n1} - 1) / {n2}).")
    print(f"The numbers in this expression are n1 = {n1}, 1, and n2 = {n2}.")
    print(f"The argument of T is alpha = {alpha_num}/{alpha_den} = {alpha}.")
    
    # We identify T(alpha) as the Beta function B(alpha, 1/2).
    # The value is B(1/4, 1/2).
    print("We calculate T(alpha) as the Beta function B(alpha, 1/2).")
    
    val = beta(alpha, 0.5)
    
    print(f"The value of T({alpha}) is {val}.")

solve_problem()