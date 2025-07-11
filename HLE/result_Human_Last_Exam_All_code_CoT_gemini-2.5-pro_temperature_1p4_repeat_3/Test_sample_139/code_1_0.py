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
    Finds the optimal resistor values and calculates the maximum possible current in R3.
    """
    R_2 = 6
    R1_sol, R3_sol = None, None

    # Based on analysis, solutions only exist for R3 > 10.
    # We search for the smallest possible values as I3 is maximized for smaller resistances.
    r3_candidate = 11
    while R3_sol is None:
        if is_prime(r3_candidate):
            # Conditions derived from z(C) and other constraints for R1 > 6, R3 > 6:
            # R1 >= (R3 + 6)/2 and R1 < R3 - 2
            r1_min = math.ceil((r3_candidate + 6) / 2)
            r1_max = r3_candidate - 3

            if r1_min <= r1_max:
                # The smallest valid R1 gives the largest current for a fixed R3.
                R1_sol = int(r1_min)
                R3_sol = r3_candidate
        r3_candidate += 1

    print("Step 1: Finding the optimal resistor values based on the problem's constraints.")
    print(f"The search for the smallest valid resistances that satisfy all conditions yields:")
    print(f"R1 = {R1_sol} ohms")
    print(f"R2 = {R_2} ohms")
    print(f"R3 = {R3_sol} ohms (which is a prime number)")
    print("\nThese values satisfy R3 - R1 > 2 (13 - 10 = 3 > 2) and z(10, 6, 13) = 6.")

    print("\nStep 2: Calculating the maximum possible current through R3.")
    
    # Assign found values
    R1, R2, R3 = R1_sol, R_2, R3_sol
    V_fail = 26

    print("The total current (I) from the source is determined from the failure condition (R2 fails open):")
    print("I = V_fail / R_eq_parallel(R1, R3)")
    print(f"I = {V_fail} / (({R1} * {R3}) / ({R1} + {R3})) = {V_fail} * ({R1} + {R3}) / ({R1} * {R3})")
    total_I = (V_fail * (R1 + R3)) / (R1 * R3)
    print(f"I = {V_fail * (R1 + R3)} / {R1 * R3} = {total_I:.2f} A")

    print("\nWhen the circuit is intact, this total current I is divided among R1, R2, and R3.")
    print("The current through R3 (I3) is found using the current divider rule:")
    print("I3 = I * (G3 / (G1 + G2 + G3)), where Gx = 1/Rx")
    print("This can be written as:")
    print(f"I3 = I * ( (1/{R3}) / (1/{R1} + 1/{R2} + 1/{R3}) )")
    
    print("\nCombining these steps into a single final equation for I3:")
    print(f"I3 = [ V_fail * (R1+R3) / (R1*R3) ] * [ (1/R3) / (1/R1 + 1/R2 + 1/R3) ]")
    
    # Calculation
    # First term: Total Current
    term1_num = V_fail * (R1 + R3)
    term1_den = R1 * R3
    
    # Second term: Current divider ratio
    G1, G2, G3 = 1/R1, 1/R2, 1/R3
    G_total = G1 + G2 + G3
    term2_num = G3
    term2_den = G_total

    final_I3_val = (term1_num / term1_den) * (term2_num / term2_den)
    
    # For integer arithmetic and showing the equation clearly
    num_part_2 = R1 * R2
    den_part_2 = (R1 * R2) + (R1 * R3) + (R2 * R3)

    final_num = term1_num * num_part_2
    final_den = term1_den * den_part_2
    common_divisor = math.gcd(final_num, final_den)

    print("\nSubstituting the numerical values into the equation:")
    print(f"I3 = [ {V_fail} * ({R1}+{R3}) / ({R1}*{R3}) ] * [ ({R1}*{R2}) / ({R1}*{R2} + {R1}*{R3} + {R2}*{R3}) ]")
    print(f"I3 = [ {V_fail * (R1+R3)} / {R1*R3} ] * [ {R1*R2} / {den_part_2} ]")
    print(f"I3 = [ {V_fail * (R1+R3)} * {R1*R2} ] / [ {R1*R3} * {den_part_2} ]")
    print(f"I3 = {final_num} / {final_den}")
    print(f"I3 = {final_num // common_divisor} / {final_den // common_divisor} Amperes")
    
    print(f"\nThe maximum possible current through R3 is approximately {final_I3_val:.3f} A.")

solve_circuit_problem()