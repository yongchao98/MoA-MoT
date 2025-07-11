import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def solve_circuit():
    """
    Finds the resistor values and calculates the maximum possible current through R3.
    """
    r2 = 6
    found = False
    r1_sol, r3_sol = 0, 0

    # Search for the smallest valid (R1, R3) pair to maximize current
    # Start with the smallest prime R3 > 6
    r3_candidate = 7
    while not found:
        if is_prime(r3_candidate):
            # R1 must be > 6 and R3 - R1 > 2 => R1 < R3 - 2
            for r1_candidate in range(7, r3_candidate - 2):
                if r1_candidate == r3_candidate:
                    continue

                # Check the z(C) condition inequalities
                # |6 - mu| > |R1 - mu|  and |6 - mu| > |R3 - mu|
                # which simplifies to:
                cond1 = (r1_candidate + r3_candidate - 12) > abs(2 * r1_candidate - r3_candidate - 6)
                cond2 = (r1_candidate + r3_candidate - 12) > abs(2 * r3_candidate - r1_candidate - 6)

                if cond1 and cond2:
                    r1_sol = r1_candidate
                    r3_sol = r3_candidate
                    found = True
                    break
        if found:
            break
        r3_candidate += 1

    print(f"Found suitable resistor values:")
    print(f"R1 = {r1_sol} ohms")
    print(f"R2 = {r2} ohms")
    print(f"R3 = {r3_sol} ohms (prime)")
    print("-" * 20)

    # Calculate the source current from the failure condition
    # I_source = 26 * (R1 + R3) / (R1 * R3)
    i_source_num = 26 * (r1_sol + r3_sol)
    i_source_den = r1_sol * r3_sol
    
    # Calculate the current through R3 when R2 is intact using the current divider rule
    # I3 = I_source * (1/R3) / (1/R1 + 1/R2 + 1/R3)
    # I3 = I_source * (R1*R2) / (R1*R2 + R1*R3 + R2*R3)
    
    # Substitute I_source
    # I3 = (26 * (R1+R3)/(R1*R3)) * (R1*R2) / (R1*R2 + R1*R3 + R2*R3)
    # I3 = 26 * R2 * (R1+R3) / (R3 * (R1*R2 + R1*R3 + R2*R3))
    
    numerator = 26 * r2 * (r1_sol + r3_sol)
    denominator = r3_sol * (r1_sol * r2 + r1_sol * r3_sol + r2 * r3_sol)
    
    current_r3 = numerator / denominator

    print("Calculation of the current I_3 through R_3:")
    print(f"I_3 = (26 * {r2} * ({r1_sol} + {r3_sol})) / ({r3_sol} * ({r1_sol}*{r2} + {r1_sol}*{r3_sol} + {r2}*{r3_sol}))")
    print(f"I_3 = ({26 * r2} * {r1_sol + r3_sol}) / ({r3_sol} * ({r1_sol * r2} + {r1_sol * r3_sol} + {r2 * r3_sol}))")
    print(f"I_3 = {numerator} / {denominator}")
    print(f"I_3 = {current_r3} Amperes")
    
    # The final answer tag requires the value directly
    # print(f'<<<{current_r3}>>>')

solve_circuit()