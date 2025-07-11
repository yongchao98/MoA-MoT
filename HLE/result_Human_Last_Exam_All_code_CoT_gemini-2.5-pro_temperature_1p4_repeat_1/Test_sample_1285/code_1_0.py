import math

def solve_b():
    """
    This function calculates the maximum number of roots for t=5.
    The formula for the maximum number of roots for a given t is 2^(t-1) - 1.
    """
    t = 5
    
    # Calculate the maximum number of roots
    # The formula is 2^(t-1) - 1
    exponent = t - 1
    power_of_2 = 2**exponent
    max_roots = power_of_2 - 1
    
    # (a) For the general case, the expression is 2^(t-1) - 1.
    # (b) For t = 5, we calculate the value.
    
    print(f"(a) 2^(t-1) - 1; (b) 2^({t}-1) - 1 = {power_of_2} - 1 = {max_roots}")

solve_b()
