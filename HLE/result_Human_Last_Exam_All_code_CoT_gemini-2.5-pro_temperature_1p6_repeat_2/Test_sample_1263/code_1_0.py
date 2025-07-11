# Maximum number of terms an expression can have in a Titan register.
TERM_LIMIT = 10

def represent_as_titan_expression(n):
    """
    Converts an integer into a Titan expression (a list of 4-bit integers).
    Returns the list of terms or raises an error if the term limit is exceeded.
    """
    if not isinstance(n, int) or n < 0:
        raise TypeError("Input must be a non-negative integer.")

    if n == 0:
        return [0]
    
    terms = []
    remaining = n
    while remaining > 0:
        term_val = min(15, remaining)
        terms.append(term_val)
        remaining -= term_val
    
    if len(terms) > TERM_LIMIT:
        raise ValueError(f"{n} requires {len(terms)} terms, which exceeds the limit of {TERM_LIMIT}.")
    
    return terms

def main():
    """
    Main function to simulate the Pandora escape velocity calculation on Titan.
    """
    print("--- Titan Architecture Feasibility Test ---")
    print("Objective: Calculate the coefficient of the v_e^2 equation.")
    print("Equation: C = (8/3) * G * pi * R^2 * rho")
    
    # Define coefficients from approximations
    # pi = 3, rho = 3, R^2 = 4, G = 13/2
    
    try:
        # Step 1: Start with 8/3, multiply by pi_coeff=3 and rho_coeff=3
        # (8/3) * 3 * 3 = 24
        val_step1 = 24
        expr_step1 = represent_as_titan_expression(val_step1)
        print(f"1. Calculating (8/3) * 3 * 3 = {val_step1}")
        print(f"   Titan expression: {expr_step1} ({len(expr_step1)} terms) -> OK")
        
        # Step 2: Multiply by R^2_coeff = 4
        # 24 * 4 = 96
        val_step2 = val_step1 * 4
        expr_step2 = represent_as_titan_expression(val_step2)
        print(f"2. Multiplying by 4 = {val_step2}")
        print(f"   Titan expression: {expr_step2} ({len(expr_step2)} terms) -> OK")

        # Step 3: Multiply by G_coeff = 13/2
        # 96 * 13/2 = 48 * 13 = 624
        val_step3 = val_step2 * 13 // 2
        print(f"3. Multiplying by 13/2 = {val_step3}")
        
        # This is the step that will fail
        expr_step3 = represent_as_titan_expression(val_step3)
        # The script will not reach here

    except ValueError as e:
        print("\n!!! CALCULATION FAILED !!!")
        print(f"Reason: {e}")
        print("\nThe intermediate products are too large for Titan's architecture.")
        print("The number of terms needed to represent the result exceeds the 10-term limit of the registers.")
        print("\nConclusion: It is not possible to calculate Pandora's escape velocity with Titan.")
        print("<<<N0>>>")

if __name__ == "__main__":
    main()