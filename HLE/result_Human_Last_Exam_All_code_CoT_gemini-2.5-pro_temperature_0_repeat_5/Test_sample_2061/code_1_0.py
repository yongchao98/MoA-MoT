import math

def solve_alpha():
    """
    This function calculates the value of alpha based on the derived formula.
    """
    # Given constants
    # T = ln(10)
    # B = (1/2) * 10^20 / 99^2
    T = math.log(10)
    B = 0.5 * (10**20) / (99**2)

    # Step 1: Calculate the term (1 - e^(-3T))
    # 1 - e^(-3*ln(10)) = 1 - 10^(-3) = 1 - 0.001 = 0.999
    term_T = 1 - math.exp(-3 * T)
    
    # Step 2: Calculate sqrt(8B)
    # sqrt(8 * (1/2) * 10^20 / 99^2) = sqrt(4 * 10^20 / 99^2) = 2 * 10^10 / 99
    sqrt_8B = math.sqrt(8 * B)

    # Step 3: Calculate alpha using the derived formula
    # alpha = (1/3) * (1 - e^(-3T)) * sqrt(8B)
    alpha = (1/3) * term_T * sqrt_8B

    # Print the components of the final equation
    print("The final equation is: alpha = (1/3) * (1 - e^(-3T)) * sqrt(8B)")
    print(f"Value of (1 - e^(-3T)): {term_T}")
    print(f"Value of B: {B}")
    print(f"Value of sqrt(8B): {sqrt_8B}")
    print(f"Final value of alpha: {alpha}")

solve_alpha()