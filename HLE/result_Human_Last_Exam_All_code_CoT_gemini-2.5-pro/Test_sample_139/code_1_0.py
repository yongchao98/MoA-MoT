import math
from fractions import Fraction

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def check_z_condition(r1, r3):
    """Checks if z({r1, 6, r3}) = 6."""
    r2 = 6
    # Condition: 6 must be the value furthest from the mean.
    # This implies 6 must be the minimum or maximum value in the set.
    is_min = (r2 < r1 and r2 < r3)
    is_max = (r2 > r1 and r2 > r3)
    if not (is_min or is_max):
        return False
    
    # Verify it is furthest from the mean
    mean = (r1 + r2 + r3) / 3.0
    dist_r2 = abs(r2 - mean)
    dist_r1 = abs(r1 - mean)
    dist_r3 = abs(r3 - mean)
    
    # Use a small tolerance for floating point comparison
    return dist_r2 >= dist_r1 and dist_r2 >= dist_r3

def find_max_current():
    """Finds the maximum possible current I3 by searching for valid resistor pairs."""
    max_i3 = Fraction(0)
    best_r1, best_r3 = None, None
    best_config = ""

    # Search a reasonable range for resistor values
    for r1 in range(1, 50):
        for r3 in range(1, 50):
            if r1 == 6 or r3 == 6 or r1 == r3:
                continue

            # Check all given conditions
            if is_prime(r3) and (r3 - r1 > 2) and check_z_condition(r1, r3):
                # This is a valid pair (r1, r3). Now calculate I3 for plausible circuits.
                
                # --- Configuration 1: All resistors in parallel ---
                # V_prime = 26 V. I_source = V_prime / Req_prime = 26 * (1/r1 + 1/r3)
                # V_intact = I_source * Req_intact = 26 * (1/r1 + 1/r3) / (1/r1 + 1/6 + 1/r3)
                # I3 = V_intact / r3
                g1, g2, g3 = Fraction(1, r1), Fraction(1, 6), Fraction(1, r3)
                v_intact_parallel = 26 * (g1 + g3) / (g1 + g2 + g3)
                i3_parallel = v_intact_parallel / r3
                
                if i3_parallel > max_i3:
                    max_i3 = i3_parallel
                    best_r1, best_r3 = r1, r3
                    best_config = "parallel"
                    
                # --- Configuration 2: R1 in series with parallel(R2, R3) ---
                # V3_prime = I_source * r3 = 26 => I_source = 26/r3
                # I3_intact = I_source * r2 / (r2 + r3)
                i3_series = Fraction(26, r3) * Fraction(6, 6 + r3)
                
                if i3_series > max_i3:
                    max_i3 = i3_series
                    best_r1, best_r3 = r1, r3
                    best_config = "series"

    return max_i3, best_r1, best_r3

# Execute the search and print the results
max_i3, r1, r3 = find_max_current()
r2 = 6

print("Step 1: Find resistor values satisfying all conditions.")
print(f"The optimal resistor values that maximize the current are found to be:")
print(f"R1 = {r1} ohms")
print(f"R2 = 6 ohms (given)")
print(f"R3 = {r3} ohms (prime, as required)")
print("-" * 30)
print("Step 2: Calculate the maximum current for the optimal configuration (all-parallel).")

# The formula for current I3 in the parallel configuration is derived from:
# I3 = V / R3, where V is the voltage across the resistors when R2 is intact.
# V = I_source * Req_intact
# I_source is found from the state when R2 fails: I_source = V_failed / Req_failed = 26 / (1/(1/R1 + 1/R3))
num = 156 * (r1 + r3)
den = r3 * (6 * r3 + r1 * r3 + 6 * r1)
final_fraction = Fraction(num, den)

print("\nThe current I3 is calculated using the formula for a parallel circuit:")
print(f"I3 = 156 * (R1 + R3) / (R3 * (6*R3 + R1*R3 + 6*R1))")
print("\nSubstituting the optimal values:")
print(f"I3 = 156 * ({r1} + {r3}) / ({r3} * (6*{r3} + {r1}*{r3} + 6*{r1}))")
print(f"I3 = 156 * ({r1 + r3}) / ({r3} * ({6*r3 + r1*r3 + 6*r1}))")
print(f"I3 = {num} / ({r3} * ({6*r3 + r1*r3 + 6*r1}))")
print(f"I3 = {num} / {den}")
print(f"I3 = {final_fraction.numerator} / {final_fraction.denominator} A")

print("\nFinal Answer:")
print(f"The maximum possible current through R3 is {float(max_i3)} A.")

<<<1.0298507462686567>>>