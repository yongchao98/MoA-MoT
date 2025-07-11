import numpy as np

def calculate_b2_for_sun(n, lambda_coeffs):
    """
    Calculates the second Betti number for a coadjoint orbit of SU(n).

    Args:
        n (int): The dimension for SU(n).
        lambda_coeffs (list or tuple): A list of n real numbers representing lambda,
                                       assumed to be in the Weyl alcove (a_1 >= ... >= a_n)
                                       and sum to zero.
    """
    if len(lambda_coeffs) != n:
        print("Error: Length of lambda_coeffs must be equal to n.")
        return

    # Check if lambda is in the Weyl alcove
    if not all(lambda_coeffs[i] >= lambda_coeffs[i+1] for i in range(n - 1)):
        print(f"Error: Lambda {lambda_coeffs} is not in the Weyl alcove.")
        return
        
    # Check if the trace is zero
    if not np.isclose(sum(lambda_coeffs), 0):
        print(f"Error: The sum of coefficients for lambda {lambda_coeffs} is not zero.")
        return

    # b2 is the number of strict inequalities a_i > a_{i+1}
    b2 = 0
    for i in range(n - 1):
        if lambda_coeffs[i] > lambda_coeffs[i+1]:
            b2 += 1
            
    rank = n - 1
    
    print(f"For G = SU({n}) and lambda = {lambda_coeffs}:")
    print(f"The second Betti number is b2 = {b2}.")
    print(f"The rank is n - 1 = {rank}.")
    
    equation_str = f"{b2} = {rank}" if b2 == rank else f"{b2} != {rank}"
    
    if b2 == rank:
        print(f"In this case, b2 equals n-1. The equation is {equation_str}.")
    else:
        print(f"In this case, b2 does NOT equal n-1. The equation is {equation_str}.")
    print("-" * 30)

# --- Demonstration for SU(3) ---
n = 3

# Case 1: Regular element (in the interior of the Weyl alcove)
# For SU(3), let lambda = (2, 0, -2). Here, a1 > a2 > a3.
lambda_regular = (2, 0, -2)
calculate_b2_for_sun(n, lambda_regular)

# Case 2: Singular element (on the boundary of the Weyl alcove)
# For SU(3), let lambda = (1, 1, -2). Here, a1 = a2 > a3.
# The corresponding orbit is the projective plane CP^2, which has b2=1.
lambda_singular = (1, 1, -2)
calculate_b2_for_sun(n, lambda_singular)
