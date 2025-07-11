import math
from fractions import Fraction

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
    Solves the resistor circuit problem by finding the resistor values
    and then calculating the maximum possible current through R3.
    """
    best_r1 = None
    best_r3 = None

    # From the problem analysis, we need to find the smallest integer R1 > 8
    # for which a prime R3 exists in the range (R1 + 2, 2*R1 - 6).
    # The current I3 is maximized for the smallest valid (R1, R3) pair.
    # We search for the first such pair, which will be the solution.
    for r1_candidate in range(9, 200):  # Search within a reasonable range for R1
        lower_bound_r3 = r1_candidate + 2
        upper_bound_r3 = 2 * r1_candidate - 6

        # If the range is invalid, continue to the next R1
        if lower_bound_r3 >= upper_bound_r3:
            continue

        for r3_candidate in range(lower_bound_r3 + 1, upper_bound_r3):
            if is_prime(r3_candidate):
                # Found the smallest valid pair (R1, R3)
                best_r1 = r1_candidate
                best_r3 = r3_candidate
                break
        
        if best_r1 is not None:
            break

    if best_r1 is None:
        print("No solution found.")
        return

    r1 = best_r1
    r3 = best_r3
    r2 = 6

    print(f"Found the resistor values that satisfy all conditions:")
    print(f"R1 = {r1} ohms")
    print(f"R2 = {r2} ohms")
    print(f"R3 = {r3} ohms (which is a prime number)\n")

    print("Calculating the current I_3 through resistor R3 with all resistors intact.")
    print("The formula for the current is: I_3 = (26 * R2 * (R1 + R3)) / (R3 * (R1*R2 + R1*R3 + R2*R3))\n")
    
    print("Substituting the values into the equation:")
    
    # Numerator calculation
    num_term1 = r1 + r3
    numerator = 26 * r2 * num_term1
    
    # Denominator calculation
    den_term1 = r1 * r2
    den_term2 = r1 * r3
    den_term3 = r2 * r3
    den_sum = den_term1 + den_term2 + den_term3
    denominator = r3 * den_sum

    print(f"I_3 = (26 * {r2} * ({r1} + {r3})) / ({r3} * ({r1}*{r2} + {r1}*{r3} + {r2}*{r3}))")
    print(f"I_3 = (156 * {num_term1}) / ({r3} * ({den_term1} + {den_term2} + {den_term3}))")
    print(f"I_3 = {numerator} / ({r3} * {den_sum})")
    print(f"I_3 = {numerator} / {denominator}\n")

    # Simplify the fraction for a more elegant result
    result_fraction = Fraction(numerator, denominator)
    print(f"As a simplified fraction, the exact current is: I_3 = {result_fraction.numerator} / {result_fraction.denominator} A.")

    # Calculate the final numerical answer
    max_current = numerator / denominator
    print(f"The maximum possible current through R3 is {max_current:.5f} A.")
    
    # Final answer in the required format
    print(f"<<<{max_current}>>>")

solve_circuit_problem()