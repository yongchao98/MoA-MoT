import math

def solve_limit():
    """
    This function calculates the integer part of 10^4 * lim_{N->inf} F(N)/ln(N).
    """

    # Constants from the characteristic equations for k=13 and k=5
    k1_num, k1_rad = 13, 165
    k2_num, k2_rad = 5, 21
    
    # The dominant roots of the characteristic equations
    phi_13 = (k1_num + math.sqrt(k1_rad)) / 2
    phi_5 = (k2_num + math.sqrt(k2_rad)) / 2

    # Natural logarithms of the roots
    log_phi_13 = math.log(phi_13)
    log_phi_5 = math.log(phi_5)

    # The limit L = 2 * (1/ln(phi_13) + 1/ln(phi_5))
    limit_val = 2 * (1 / log_phi_13 + 1 / log_phi_5)
    
    # The final value is the integer part of 10^4 * L
    result = int(10**4 * limit_val)

    # Print the equation with all the numbers
    print("The limit L is calculated by the formula:")
    print(f"L = 2 * (1/ln(({k1_num} + sqrt({k1_rad}))/2) + 1/ln(({k2_num} + sqrt({k2_rad}))/2))")
    print("\nIntermediate values:")
    print(f"phi_13 = ({k1_num} + sqrt({k1_rad}))/2 = {phi_13}")
    print(f"phi_5 = ({k2_num} + sqrt({k2_rad}))/2 = {phi_5}")
    print(f"ln(phi_13) = {log_phi_13}")
    print(f"ln(phi_5) = {log_phi_5}")
    print(f"L = 2 * (1/{log_phi_13} + 1/{log_phi_5}) = {limit_val}")

    # Print the final calculation
    print(f"\nThe value of 10^4 * L is {10**4 * limit_val}")
    print(f"The integer part is {result}")

solve_limit()