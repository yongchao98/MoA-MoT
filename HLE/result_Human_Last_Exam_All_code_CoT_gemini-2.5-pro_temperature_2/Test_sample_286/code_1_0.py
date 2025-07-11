import math

def solve_m():
    """
    This function calculates the smallest integer M based on the derivation.
    M must satisfy M >= 100000 * log10(2).
    """
    N = 100000
    
    # Calculate the lower bound for M
    lower_bound = N * math.log10(2)
    
    # M must be the smallest integer greater than or equal to this value.
    # math.ceil will find this for us.
    M = math.ceil(lower_bound)
    
    log_val = math.log10(2)

    print(f"Let N = 100000.")
    print(f"The condition for M is M >= N * log10(2).")
    print(f"log10(2) is approximately {log_val}")
    print(f"N * log10(2) is approximately {lower_bound}")
    print(f"The smallest integer M satisfying this is {M}.")

solve_m()