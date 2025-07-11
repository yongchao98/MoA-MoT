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

def find_max_current():
    """
    Finds the resistor pair (R1, R3) that maximizes the current I3
    and prints the step-by-step calculation.
    """
    R2 = 6
    V_fail = 26
    max_current = 0
    best_resistors = None

    # Analysis shows current I3 is maximized for the smallest valid R1 and R3.
    # We search for the smallest prime R3, and for each, the smallest valid R1.
    # The first valid pair found will be our answer.
    # The search range for R3 is arbitrary; we expect a solution quickly.
    for R3_val in range(7, 100):
        if not is_prime(R3_val):
            continue

        # Derived constraints from the problem statement:
        # R1 >= (R3 + 6) / 2
        # R1 < R3 - 2
        # R1 must be an integer > 6
        
        # Calculate integer bounds for the R1 search loop
        r1_lower_bound = math.ceil((R3_val + R2) / 2.0)
        r1_upper_bound = R3_val - 2
        
        for R1_val in range(r1_lower_bound, r1_upper_bound):
            # This is a valid pair satisfying all conditions.
            # Since we are iterating R3 and R1 from low to high, this first
            # pair is the one that minimizes the resistances.
            best_resistors = (R1_val, R3_val)
            break # Found the smallest R1 for this R3
        
        if best_resistors:
            break # Found the smallest R3 with a valid R1

    if not best_resistors:
        print("No solution found within the search range.")
        return

    R1_final, R3_final = best_resistors
    
    # Use Fraction for precise calculations
    R1_f = Fraction(R1_final)
    R3_f = Fraction(R3_final)
    R2_f = Fraction(R2)
    V_fail_f = Fraction(V_fail)

    # Calculate source current Is from the failed-R2 state
    Is_final = V_fail_f * (R1_f + R3_f) / (R1_f * R3_f)
    
    # Calculate conductances for the intact circuit
    G1_final = Fraction(1, R1_f)
    G2_final = Fraction(1, R2_f)
    G3_final = Fraction(1, R3_f)
    total_G_final = G1_final + G2_final + G3_final
    
    # Calculate I3 using the current divider rule
    I3_final = Is_final * G3_final / total_G_final

    # Print the detailed calculation as requested
    print(f"Based on the problem's constraints, the maximum current occurs for the smallest valid resistance values.")
    print(f"The optimal values are found to be R1 = {R1_final} Ω and R3 = {R3_final} Ω (with R2 = {R2} Ω).\n")
    print("--- Calculation for the Maximum Current ---")

    print(f"\n1. First, calculate the source current, I_s:")
    print(f"   When R2 fails, the voltage across R1 and R3 is {V_fail} V.")
    print(f"   V_fail = I_s * (R1 || R3) = I_s * (R1 * R3) / (R1 + R3)")
    print(f"   {V_fail} = I_s * ({R1_final} * {R3_final}) / ({R1_final} + {R3_final})")
    print(f"   {V_fail} = I_s * {R1_f * R3_f}")
    print(f"   I_s = {V_fail} * ({R1_final} + {R3_final}) / ({R1_final * R3_final}) = {V_fail * (R1_final + R3_final) / (R1_final*R3_final):.1f} A")

    print(f"\n2. Next, calculate the current through R3 (I_3) when the circuit is intact:")
    print(f"   Using the current divider rule: I_3 = I_s * (1/R3) / (1/R1 + 1/R2 + 1/R3)")
    print(f"   The total conductance is 1/{R1_final} + 1/{R2} + 1/{R3_final} = {G1_final.numerator}/{G1_final.denominator} + {G2_final.numerator}/{G2_final.denominator} + {G3_final.numerator}/{G3_final.denominator} = {total_G_final.numerator}/{total_G_final.denominator} S")
    print(f"   I_3 = {float(Is_final):.1f} A * (1/{R3_final}) / ({total_G_final.numerator}/{total_G_final.denominator})")
    
    # Final exact value as fraction
    print(f"\n   The exact current is I_3 = {Is_final.numerator}/{Is_final.denominator} * ({G3_final.numerator}/{G3_final.denominator}) / ({total_G_final.numerator}/{total_G_final.denominator}) = {I3_final.numerator}/{I3_final.denominator} A.")
    
    # Final value as float for the answer
    print(f"\nMaximum possible current through R3: {float(I3_final):.4f} A")

    global final_answer
    final_answer = float(I3_final)


if __name__ == '__main__':
    final_answer = 0.0
    find_max_current()
    # The final answer is wrapped in <<<>>> as requested.
    print(f"\n<<<{final_answer:.4f}>>>")
