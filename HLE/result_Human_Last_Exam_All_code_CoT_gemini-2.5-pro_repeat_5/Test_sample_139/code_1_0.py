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
    This function solves for the maximum possible current through R3 based on the problem's constraints.
    
    1.  The condition z({R1, 6, R3}) = 6 means 6 is the element furthest from the mean.
        Given R3 - R1 > 2 and R3 is a prime, analysis shows that 6 must be the minimum value, so R1 > 6 and R3 > 6.
    
    2.  The condition that 6 is the outlier leads to the inequalities:
        R3 + 6 <= 2 * R1  and  R1 + 6 <= 2 * R3.
    
    3.  Combining with R3 - R1 > 2 (or R1 < R3 - 2), we get the search range for R1:
        (R3 + 6) / 2 <= R1 < R3 - 2.
        For a solution to exist, (R3 + 6) / 2 < R3 - 2, which simplifies to R3 > 10.
    
    4.  The current I3 is derived from circuit principles:
        When R2 fails, V_across_R3 = I_total * (R1 * R3) / (R1 + R3) = 26 V.
        So, I_total = 26 * (R1 + R3) / (R1 * R3).
        When R2 is intact, the voltage is V = I_total * R_equivalent, where R_equivalent for three parallel resistors is
        (1/R1 + 1/R2 + 1/R3)^-1. With R2=6, this is V = I_total * (6*R1*R3) / (6*R1 + R1*R3 + 6*R3).
        The current through R3 is I3 = V / R3.
        Substituting I_total and V, we get:
        I3 = 156 * (R1 + R3) / (R3 * (6*R1 + R1*R3 + 6*R3)).
    
    5.  We search for integer R1 and prime R3 that satisfy the conditions and maximize I3.
        The current I3 generally decreases as resistances increase. Thus, the maximum value is expected for the
        smallest possible R1 and R3.
    """
    
    R2 = 6
    max_i3 = 0.0
    best_r1 = None
    best_r3 = None

    # We search for primes R3 > 10. A reasonable upper limit is sufficient to show the trend.
    for r3_candidate in range(11, 200):
        if not is_prime(r3_candidate):
            continue
        
        # Determine the valid range for R1
        r1_min = math.ceil((r3_candidate + 6) / 2.0)
        r1_max = r3_candidate - 3 # R1 < R3 - 2 means R1 <= R3 - 3 for integers
        
        for r1_candidate in range(r1_min, r1_max + 1):
            # R1 must be an integer and distinct from R3 and R2
            if r1_candidate == r3_candidate or r1_candidate == R2:
                continue

            # Calculate I3 for this valid pair (r1_candidate, r3_candidate)
            R1 = r1_candidate
            R3 = r3_candidate
            
            numerator = 156 * (R1 + R3)
            denominator = R3 * (6 * R1 + R1 * R3 + 6 * R3)
            current_i3 = numerator / denominator

            if current_i3 > max_i3:
                max_i3 = current_i3
                best_r1 = R1
                best_r3 = R3

    # Output the results
    if best_r1 is not None:
        R1 = best_r1
        R3 = best_r3
        
        print("The maximum current is achieved with the following resistor values:")
        print(f"R1 = {R1} ohms")
        print(f"R2 = {R2} ohms")
        print(f"R3 = {R3} ohms (which is a prime number)")
        print("\nVerification of conditions:")
        print(f"R3 - R1 = {R3 - R1} > 2. Condition met.")
        mean = (R1 + R2 + R3) / 3
        dist1 = abs(R1 - mean)
        dist2 = abs(R2 - mean)
        dist3 = abs(R3 - mean)
        print(f"Distances from mean ({mean:.2f}): |R1-μ|={dist1:.2f}, |R2-μ|={dist2:.2f}, |R3-μ|={dist3:.2f}.")
        print(f"R2=6 is the element furthest from the mean. Condition z=6 met.")

        print("\nCalculation of the maximum current through R3:")
        num_val = 156 * (R1 + R3)
        den_val_part1 = 6 * R1 + R1 * R3 + 6 * R3
        den_val = R3 * den_val_part1
        
        print(f"I_3 = (156 * (R1 + R3)) / (R3 * (6*R1 + R1*R3 + 6*R3))")
        print(f"I_3 = (156 * ({R1} + {R3})) / ({R3} * (6*{R1} + {R1}*{R3} + 6*{R3}))")
        print(f"I_3 = {num_val} / ({R3} * ({den_val_part1}))")
        print(f"I_3 = {num_val} / {den_val}")
        
        # For an exact fraction
        common_divisor = math.gcd(num_val, den_val)
        print(f"I_3 = {num_val // common_divisor} / {den_val // common_divisor} Amperes")
        
        print(f"\nMaximum possible current through R3 is: {max_i3:.6f} A")

    else:
        print("No solution found within the search range.")

solve_circuit_problem()