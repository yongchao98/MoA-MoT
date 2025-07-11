import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def find_resistor_values():
    """
    Finds the values of R1 and R3 based on the problem constraints.
    """
    solutions = []
    r2 = 6

    # Case A: R2 is the minimum value (R1 > 6, R3 > 6)
    # The condition z(C)=6 implies 6 is furthest from the mean.
    # If 6 is the min, this gives the inequality: 2*R1 - 6 > R3
    for r3_candidate in range(r2 + 1, 100): # Search for R3 in a reasonable range
        if is_prime(r3_candidate):
            # Constraint: R3 - R1 > 2  => R1 < R3 - 2
            for r1_candidate in range(r2 + 1, r3_candidate - 2):
                # Constraint from z(C)=6 for Case A
                if 2 * r1_candidate - 6 > r3_candidate:
                    # All constraints met for Case A
                    solutions.append((r1_candidate, r3_candidate))

    # Case B: R2 is the maximum value (R1 < 6, R3 < 6)
    # The condition z(C)=6 implies 6 is furthest from the mean.
    # If 6 is the max, this gives the inequalities: R3+6 > 2*R1 and R1+6 > 2*R3
    for r3_candidate in range(2, r2):
        if is_prime(r3_candidate):
            # Constraint: R3 - R1 > 2 => R1 < R3 - 2
            for r1_candidate in range(1, r3_candidate - 2):
                # Constraint from z(C)=6 for Case B
                if (r3_candidate + 6 > 2 * r1_candidate) and (r1_candidate + 6 > 2 * r3_candidate):
                    solutions.append((r1_candidate, r3_candidate))

    return solutions

def calculate_max_current():
    """
    Calculates the maximum possible current through R3.
    """
    print("Step 1 & 2: Finding resistor values based on constraints...")
    
    possible_solutions = find_resistor_values()
    
    if not possible_solutions:
        print("No valid resistor values found.")
        return

    # We expect a unique solution from the problem statement
    r1, r3 = possible_solutions[0]
    r2 = 6
    
    print(f"Found a unique solution for the resistor values: R1 = {r1} ohms, R2 = {r2} ohms, R3 = {r3} ohms.")
    print("\nStep 3 & 4: Analyzing circuit scenarios and calculating currents.")

    v3_failed = 26.0

    # Scenario 1: R1 in series with a parallel combination of R2 and R3
    # When R2 fails (open), V3_failed = I_source * R3
    i_source_1 = v3_failed / r3
    # When intact, I3 = I_source * (R2 / (R2 + R3)) (current divider)
    i3_scenario1 = i_source_1 * (r2 / (r2 + r3))
    
    print("\nScenario 1: R1 in series with (R2 || R3)")
    print(f"Source Current (I) = V3_failed / R3 = {v3_failed} / {r3} = {i_source_1:.4f} A")
    print(f"Current through R3 (I3) = I * (R2 / (R2 + R3)) = {i_source_1:.4f} * ({r2} / ({r2} + {r3})) = {i3_scenario1:.4f} A")

    # Scenario 2: R1, R2, R3 in parallel
    # When R2 fails, V_parallel = V3_failed = 26V.
    # I_source = V_parallel/R1 + V_parallel/R3
    i_source_2 = (v3_failed / r1) + (v3_failed / r3)
    # When intact, the total resistance is Req = 1 / (1/R1 + 1/R2 + 1/R3)
    # The voltage is V = I_source * Req
    # The current through R3 is I3 = V / R3
    req_intact_inv = (1/r1) + (1/r2) + (1/r3)
    v_intact = i_source_2 / req_intact_inv
    i3_scenario2 = v_intact / r3
    
    print("\nScenario 2: R1, R2, and R3 in parallel")
    print(f"Source Current (I) = V3_failed/R1 + V3_failed/R3 = {v3_failed}/{r1} + {v3_failed}/{r3} = {i_source_2:.4f} A")
    print(f"Voltage across intact circuit (V) = I / (1/R1 + 1/R2 + 1/R3) = {i_source_2:.4f} / ({1/r1:.4f} + {1/r2:.4f} + {1/r3:.4f}) = {v_intact:.4f} V")
    print(f"Current through R3 (I3) = V / R3 = {v_intact:.4f} / {r3} = {i3_scenario2:.4f} A")
    
    # Step 5: Determine the maximum possible current
    max_i3 = max(i3_scenario1, i3_scenario2)
    
    print("\nStep 5: Determining the maximum possible current.")
    print(f"Comparing the two scenarios: {i3_scenario1:.4f} A vs {i3_scenario2:.4f} A.")
    print(f"\nThe maximum possible current through R3 when R2 is intact is {max_i3:.4f} A.")

if __name__ == '__main__':
    calculate_max_current()