import math

def solve_circuit_problem():
    """
    Finds the maximum possible current through R3 under the given conditions.
    """
    R2 = 6
    V3_fail = 26
    max_I3 = 0.0
    best_params = {}

    def is_prime(n):
        if n <= 1:
            return False
        for i in range(2, int(math.sqrt(n)) + 1):
            if n % i == 0:
                return False
        return True

    def z_is_6(r1, r2, r3):
        # For a 3-element set, the element furthest from the mean is the min or max.
        is_min = r2 < r1 and r2 < r3
        is_max = r2 > r1 and r2 > r3
        return is_min or is_max

    # We search a reasonable range for R1 and R3.
    # The current tends to decrease as resistance values increase,
    # so the maximum will be found with small resistor values.
    for R1 in range(1, 100):
        for R3 in range(1, 100):
            # Condition 1: Distinct integer values
            if R1 == R2 or R3 == R2 or R1 == R3:
                continue

            # Condition 2: R3 is a prime number
            if not is_prime(R3):
                continue

            # Condition 3: R3 - R1 > 2
            if not (R3 - R1 > 2):
                continue

            # Condition 4: z(R1, R2, R3) = 6
            if not z_is_6(R1, R2, R3):
                continue

            # If conditions are met, we have a valid set of resistors.
            # Now, calculate I3 for the two possible circuit configurations.

            # Config 1: R1, R2, R3 in Parallel
            # After R2 fails (open), V3 = I_source * Req_fail = 26
            # Req_fail = (R1 * R3) / (R1 + R3)
            # I_source = V3_fail / Req_fail = 26 * (R1 + R3) / (R1 * R3)
            # Before failure, V_intact = I_source * Req_intact
            # Req_intact = 1 / (1/R1 + 1/R2 + 1/R3)
            # I3_parallel = V_intact / R3
            I_source_p = V3_fail * (R1 + R3) / (R1 * R3)
            Req_intact_p = 1 / (1/R1 + 1/R2 + 1/R3)
            V_intact_p = I_source_p * Req_intact_p
            I3_parallel = V_intact_p / R3
            
            if I3_parallel > max_I3:
                max_I3 = I3_parallel
                best_params = {'R1': R1, 'R2': R2, 'R3': R3, 'config': 'Parallel', 'I3': I3_parallel}

            # Config 2: R1 in series with (R2 || R3)
            # After R2 fails (open), V3 = I_source * R3 = 26
            # I_source = V3_fail / R3
            # Before failure, I3 is found by the current divider rule:
            # I3_series_parallel = I_source * R2 / (R2 + R3)
            I_source_sp = V3_fail / R3
            I3_series_parallel = I_source_sp * R2 / (R2 + R3)

            if I3_series_parallel > max_I3:
                max_I3 = I3_series_parallel
                best_params = {'R1': R1, 'R2': R2, 'R3': R3, 'config': 'Series-Parallel', 'I3': I3_series_parallel}

    # Print the results for the optimal case found
    print("The maximum possible current is achieved under the following conditions:")
    print(f"Resistor values: R1 = {best_params['R1']} Ω, R2 = {best_params['R2']} Ω, R3 = {best_params['R3']} Ω")
    print(f"Circuit configuration: {best_params['config']}")
    print("\nCalculation of the current through R3 (I3) before failure:")

    R1 = best_params['R1']
    R2 = best_params['R2']
    R3 = best_params['R3']
    
    if best_params['config'] == 'Parallel':
        # Step 1: Find the source current I_source from the failure condition
        # V3_fail = I_source * Req_fail => I_source = V3_fail / Req_fail
        Req_fail = (R1 * R3) / (R1 + R3)
        I_source = V3_fail / Req_fail
        print(f"1. When R2 fails, the equivalent resistance is (R1*R3)/(R1+R3) = ({R1}*{R3})/({R1}+{R3}) = {Req_fail:.4f} Ω.")
        print(f"   The source current is I_source = V3_fail / Req_fail = {V3_fail} / {Req_fail:.4f} = {I_source:.4f} A.")

        # Step 2: Calculate the voltage across the intact circuit
        # V_intact = I_source * Req_intact
        Req_intact = 1 / (1/R1 + 1/R2 + 1/R3)
        V_intact = I_source * Req_intact
        print(f"2. With all resistors intact, the equivalent resistance is 1/(1/R1 + 1/R2 + 1/R3) = 1/(1/{R1} + 1/{R2} + 1/{R3}) = {Req_intact:.4f} Ω.")
        print(f"   The voltage across the intact parallel circuit is V_intact = I_source * Req_intact = {I_source:.4f} * {Req_intact:.4f} = {V_intact:.4f} V.")

        # Step 3: Calculate the current through R3
        # I3 = V_intact / R3
        I3 = V_intact / R3
        print(f"3. The current through R3 is I3 = V_intact / R3 = {V_intact:.4f} / {R3} = {I3:.4f} A.")
        print("\nFinal equation for I3:")
        print(f"I3 = (({V3_fail} * ({R1} + {R3})) / ({R1} * {R3})) * (1 / (1/{R1} + 1/{R2} + 1/{R3})) / {R3}")
        
    else: # Series-Parallel
        # Step 1: Find the source current I_source from the failure condition
        I_source = V3_fail / R3
        print(f"1. When R2 fails, the current flows through R1 and R3 in series. The source current is I_source = V3_fail / R3 = {V3_fail} / {R3} = {I_source:.4f} A.")

        # Step 2: Use current divider to find I3 in the intact circuit
        I3 = I_source * (R2 / (R2 + R3))
        print(f"2. With all resistors intact, I_source splits between R2 and R3. The current through R3 is I3 = I_source * (R2 / (R2 + R3))")
        print(f"   I3 = {I_source:.4f} * ({R2} / ({R2} + {R3})) = {I3:.4f} A.")

    print(f"\nThe maximum possible current is {best_params['I3']:.4f} A.")

solve_circuit_problem()
<<<4.5659>>>