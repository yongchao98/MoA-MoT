import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def find_max_current():
    """
    Finds the maximum possible current through R3 by searching for valid
    resistor values and calculating the current for possible circuit scenarios.
    """
    R2 = 6
    max_current = 0.0
    best_params = {}

    # Search for primes R3 up to a reasonable limit. The current functions
    # decrease with resistance, so the max will be found with small R values.
    for R3 in range(2, 200):
        if not is_prime(R3):
            continue

        # Search for R1 values that satisfy the R3 - R1 > 2 constraint.
        for R1 in range(1, R3 - 2):
            # Constraint: R1, R2, R3 must be distinct.
            if R1 == R2 or R1 == R3 or R2 == R3:
                continue
            
            # Constraint: z({R1, R2, R3}) = 6
            # This means 6 is the element furthest from the mean.
            # This can only be true if 6 is the min or max value in the set.
            # Let's check this and the distance from the mean.
            
            resistors = [R1, R2, R3]
            mean = sum(resistors) / 3
            
            dist_R1 = abs(R1 - mean)
            dist_R2 = abs(R2 - mean)
            dist_R3 = abs(R3 - mean)
            
            # Check if R2 is the element furthest from the mean
            if dist_R2 >= dist_R1 and dist_R2 >= dist_R3:
                # This set of (R1, R3) is valid. Now calculate I3 for both scenarios.

                # Scenario 1: R1, R2, R3 in parallel
                # I3 = 26 * (6*(R1+R3)) / (6*R1 + R1*R3 + 6*R3)
                # Numerator = 26 * 6 * (R1 + R3)
                # Denominator = R1*R2 + R1*R3 + R2*R3, with R2=6
                
                # Equation for I_source from post-failure condition: 
                # I_source * (R1 * R3 / (R1 + R3)) = 26
                # Equation for I3 with intact circuit:
                # I3 = V_parallel / R3 = (I_source * R_eq) / R3
                # After substitution and simplification:
                current_scen1 = (26 * R2 * (R1 + R3)) / (R1 * R2 + R1 * R3 + R2 * R3)

                if current_scen1 > max_current:
                    max_current = current_scen1
                    best_params = {
                        "scenario": "Parallel",
                        "R1": R1, "R2": R2, "R3": R3,
                        "current": current_scen1
                    }

                # Scenario 2: (R2 || R3) in series with R1
                # Equation for I_source from post-failure condition: I_source * R3 = 26
                # Equation for I3 with intact circuit (current divider):
                # I3 = I_source * (R2 / (R2 + R3))
                # After substitution and simplification:
                current_scen2 = (26 / R3) * (R2 / (R2 + R3))

                if current_scen2 > max_current:
                    max_current = current_scen2
                    best_params = {
                        "scenario": "Combination (R2||R3 in series with R1)",
                        "R1": R1, "R2": R2, "R3": R3,
                        "current": current_scen2
                    }

    print("Found optimal resistor values that satisfy all conditions.")
    R1 = best_params['R1']
    R2 = best_params['R2']
    R3 = best_params['R3']
    print(f"Resistor values: R1 = {R1} Ohms, R2 = {R2} Ohms, R3 = {R3} Ohms.")
    print(f"The maximum current occurs in the '{best_params['scenario']}' circuit configuration.")

    if best_params['scenario'] == 'Parallel':
        numerator = 26 * R2 * (R1 + R3)
        denominator = R1*R2 + R1*R3 + R2*R3
        print(f"\nThe calculation for the current I_3 is:")
        print(f"I_3 = (26 * R2 * (R1 + R3)) / (R1*R2 + R1*R3 + R2*R3)")
        print(f"I_3 = (26 * {R2} * ({R1} + {R3})) / ({R1}*{R2} + {R1}*{R3} + {R2}*{R3})")
        print(f"I_3 = {numerator} / {denominator}")
    else: # Combination scenario
        numerator = 26 * R2
        denominator = R3 * (R2 + R3)
        print(f"\nThe calculation for the current I_3 is:")
        print(f"I_3 = (26 / R3) * (R2 / (R2 + R3))")
        print(f"I_3 = (26 / {R3}) * ({R2} / ({R2} + {R3}))")
        print(f"I_3 = {numerator} / {denominator}")

    print(f"\nMaximum possible current through R3: {max_current:.4f} A")

find_max_current()
<<<13.3881>>>