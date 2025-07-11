import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def find_resistors_and_current():
    """
    Finds the resistor values that satisfy the problem conditions
    and calculates the maximum possible current through R3.
    """
    max_current = 0
    best_resistors = None

    # We search for R3 as a prime number, starting from 7 (since R3 > 6).
    # The current decreases as resistance values increase, so the first solution
    # found by iterating R3 and R1 upwards should give the maximum current.
    # We will limit the search range for practical purposes.
    for r3 in range(7, 100):
        if not is_prime(r3):
            continue

        # Condition: R3 - R1 > 2 => R1 < R3 - 2
        # Condition: R1 > 6 as 6 is the minimum resistance
        for r1 in range(7, r3 - 2):
            if r1 == r3:  # Must be distinct
                continue

            # The resistor set is {r1, 6, r3}.
            # It's already ordered since r1 > 6 and r3 > r1.
            # Ordered set: {r_a, r_b, r_c} = {6, r1, r3}
            # Condition for z(C) = 6 (the minimum value): 2*r_b > r_a + r_c
            if 2 * r1 > 6 + r3:
                # We found a valid set of resistors
                best_resistors = {'R1': r1, 'R2': 6, 'R3': r3}
                
                # Calculate the current for this set
                numerator = 156 * (r1 + r3)
                denominator = r3 * (6 * r1 + 6 * r3 + r1 * r3)
                max_current = numerator / denominator

                # Print the results and break the loops
                print(f"Found a valid set of resistor values satisfying all conditions:")
                print(f"R1 = {best_resistors['R1']} ohms, R2 = {best_resistors['R2']} ohms, R3 = {best_resistors['R3']} ohms (prime)")
                print("\nVerification of conditions:")
                print(f"Distinct integers: Yes, {best_resistors['R1']}, {best_resistors['R2']}, {best_resistors['R3']} are distinct.")
                print(f"R3 - R1 > 2: {best_resistors['R3']} - {best_resistors['R1']} = {best_resistors['R3']-best_resistors['R1']} > 2. Yes.")
                
                # Verify z(C)=6 condition
                mean = (r1 + 6 + r3) / 3
                dist_r1 = abs(r1 - mean)
                dist_r2 = abs(6 - mean)
                dist_r3 = abs(r3 - mean)
                print(f"Mean of {{6, {r1}, {r3}}} is {mean:.2f}. Distances: |6-mean|={dist_r2:.2f}, |{r1}-mean|={dist_r1:.2f}, |{r3}-mean|={dist_r3:.2f}.")
                print(f"The element furthest from the mean is 6, so z(C)=6. Yes.")

                print("\nCalculating the current through R3 (I_3) when the circuit is intact:")
                print("Formula: I_3 = (26 * 6 * (R1 + R3)) / (R3 * (6*R1 + 6*R3 + R1*R3))")
                num_calc_1 = 26 * 6
                num_calc_2 = r1 + r3
                print(f"I_3 = ({num_calc_1} * ({r1} + {r3})) / ({r3} * (6*{r1} + 6*{r3} + {r1}*{r3}))")
                den_calc_1 = 6 * r1
                den_calc_2 = 6 * r3
                den_calc_3 = r1 * r3
                print(f"I_3 = ({num_calc_1} * {num_calc_2}) / ({r3} * ({den_calc_1} + {den_calc_2} + {den_calc_3}))")
                final_num = num_calc_1 * num_calc_2
                final_den = r3 * (den_calc_1 + den_calc_2 + den_calc_3)
                print(f"I_3 = {final_num} / {final_den}")

                # Simplify fraction for clarity
                common_divisor = math.gcd(final_num, final_den)
                print(f"I_3 = {final_num // common_divisor} / {final_den // common_divisor} Amperes")
                print(f"\nThe maximum possible current through R3 is approximately {max_current:.5f} A.")

                return max_current
    
    # In case no solution is found
    if not best_resistors:
        print("No solution found within the search range.")
        return None

final_current = find_resistors_and_current()
# The required output format is just the answer value in <<<...>>>
# As the code already prints the detailed explanation, we will just output the final answer string here.
if final_current is not None:
    print(f"\n<<<{final_current:.5f}>>>")
