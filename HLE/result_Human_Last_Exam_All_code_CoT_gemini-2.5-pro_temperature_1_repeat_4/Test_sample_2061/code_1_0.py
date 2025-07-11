import math

def solve_alpha():
    """
    This function calculates the value of alpha based on the derived formula.
    """
    # Given constants
    T = math.log(10)
    B = 0.5 * (10**20) / (99**2)

    # From the plan (Step 5), the simplified relation is C = sqrt(8*B)
    # where C = 3 * alpha / (1 - exp(-3*T))
    
    # Calculate C
    C = math.sqrt(8 * B)
    
    # Calculate the term involving T
    # 1 - exp(-3T) = 1 - (e^T)^-3 = 1 - 10^-3 = 1 - 0.001 = 0.999
    term_T = 1 - math.exp(-3 * T)
    
    # Solve for alpha
    # alpha = C * (1 - exp(-3T)) / 3
    alpha = C * term_T / 3
    
    # Let's print the intermediate and final values for clarity
    print(f"T = {T}")
    print(f"B = {B}")
    print(f"Calculated C = sqrt(8B) = {C}")
    print(f"Term (1 - e^(-3T)) = {term_T}")
    print(f"alpha = C * (1 - e^(-3T)) / 3 = {alpha}")
    
    # The final equation to find alpha is alpha = (1 - e^(-3T)) * sqrt(8B) / 3
    # Let's show the numbers in the equation
    # We have B = 1/2 * 10^20 / 99^2
    # So 8B = 4 * 10^20 / 99^2
    # sqrt(8B) = sqrt(4 * 10^20 / 99^2) = 2 * 10^10 / 99
    sqrt_8B_val = 2 * (10**10) / 99
    
    # Final alpha = term_T * sqrt_8B_val / 3
    # alpha = (1 - e^(-3*ln(10))) * (2 * 10^10 / 99) / 3
    
    print(f"Final equation: alpha = ({term_T}) * sqrt(8 * {B}) / 3")
    print(f"alpha = ({term_T}) * {sqrt_8B_val} / 3")
    print(f"alpha = {alpha}")

solve_alpha()