import math

def solve_circuit_problem():
    """
    Finds the maximum possible current through R3 by checking all valid resistor
    combinations and two possible circuit models.
    """
    r2 = 6
    max_current = 0.0
    best_params = {}

    def is_prime(n):
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

    # We expect smaller resistances to yield higher currents, so a limited search range is sufficient.
    for r3_candidate in range(2, 100):
        if not is_prime(r3_candidate):
            continue

        for r1_candidate in range(1, 100):
            # Condition: R1, R2, R3 must be distinct integers
            if r1_candidate == r3_candidate or r1_candidate == r2 or r3_candidate == r2:
                continue

            # Condition: R3 - R1 > 2
            if not (r3_candidate - r1_candidate > 2):
                continue
            
            # Condition: R2=6 is an extremum (min or max)
            is_extremum = (r1_candidate < r2 and r3_candidate < r2) or \
                          (r1_candidate > r2 and r3_candidate > r2)
            if not is_extremum:
                continue

            # This is a valid (R1, R3) pair. Now, calculate I3 for both circuit models.
            r1, r3 = r1_candidate, r3_candidate

            # Model A: Parallel Circuit
            # V_R3_failed = I * (r1 * r3) / (r1 + r3) = 26
            # I = 26 * (r1 + r3) / (r1 * r3)
            # I3_intact = I * (r1 * r2) / (r1*r2 + r1*r3 + r2*r3)
            current_a = (26 * (r1 + r3) / (r1 * r3)) * \
                        (r1 * r2) / (r1*r2 + r1*r3 + r2*r3)
            
            if current_a > max_current:
                max_current = current_a
                best_params = {
                    'model': 'A', 'r1': r1, 'r3': r3, 'current': current_a
                }

            # Model B: Series-Parallel Circuit (R1 in series with R2 || R3)
            # V_R3_failed = I * r3 = 26
            # I = 26 / r3
            # I3_intact = I * r2 / (r2 + r3)
            current_b = (26 / r3) * (r2 / (r2 + r3))

            if current_b > max_current:
                max_current = current_b
                best_params = {
                    'model': 'B', 'r1': r1, 'r3': r3, 'current': current_b
                }
    
    # After checking all combinations, print the calculation for the maximum case found.
    r1 = best_params['r1']
    r3 = best_params['r3']

    print("The maximum possible current occurs in a parallel circuit configuration.")
    print(f"The resistor values that maximize the current are R1 = {r1} ohms, R2 = {r2} ohms, and R3 = {r3} ohms.")
    print("\nIn this model, when R2 fails, the circuit simplifies to R1 and R3 in parallel with the current source.")
    print(f"The voltage across R3 is given as 26 V, so: I * (R1 * R3) / (R1 + R3) = 26.")
    print("\nWhen the circuit is intact, the current through R3 (I3) is calculated as follows:")
    
    # Calculate intermediate values for clarity
    r1_plus_r3 = r1 + r3
    r1_times_r3 = r1 * r3
    
    term1_I = f"26 * ({r1} + {r3}) / ({r1} * {r3})"
    
    r1_r2 = r1 * r2
    r1_r3_term = r1 * r3
    r2_r3 = r2 * r3
    denominator_sum = r1_r2 + r1_r3_term + r2_r3

    term2_dist = f"({r1} * {r2}) / ({r1}*{r2} + {r1}*{r3} + {r2}*{r3})"

    print(f"I3 = [Source Current I] * [Current Divider Ratio]")
    print(f"I3 = [{term1_I}] * [{term2_dist}]")
    print(f"I3 = [26 * {r1_plus_r3} / {r1_times_r3}] * [{r1_r2} / ({r1_r2} + {r1_r3_term} + {r2_r3})]")
    print(f"I3 = [26 * {r1_plus_r3} / {r1_times_r3}] * [{r1_r2} / {denominator_sum}]")

    final_numerator = 26 * r1_plus_r3 * r1_r2
    final_denominator = r1_times_r3 * denominator_sum
    
    print(f"I3 = {final_numerator} / {final_denominator}")
    print(f"I3 = {best_params['current']}")


solve_circuit_problem()
<<<936/205>>>