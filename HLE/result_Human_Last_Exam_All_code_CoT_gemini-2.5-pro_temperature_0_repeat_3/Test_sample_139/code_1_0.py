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

def solve_and_print_max_current():
    """
    Finds the resistor pair (R1, R3) that satisfies all problem constraints
    and maximizes the current through R3, then prints the detailed calculation.
    """
    R2 = 6
    V_fail = 26
    valid_solutions = []

    # Analysis shows that for a solution to exist, R1 > 8.
    # We search a reasonable range for R1 to find valid pairs.
    # The condition z(C)=6 for the case R2=6 being the minimum and R3 > R1
    # simplifies to R1 + 2 < R3 < 2*R1 - 6.
    for R1 in range(9, 100):
        # Iterate through possible R3 values based on the derived inequality
        for R3 in range(R1 + 3, 2 * R1 - 6):
            # Check if R3 is prime
            if is_prime(R3):
                # All conditions are met by this construction.
                # Calculate the current I3 for this valid pair.
                R1_f, R2_f, R3_f = Fraction(R1), Fraction(R2), Fraction(R3)
                V_fail_f = Fraction(V_fail)

                I_total_f = V_fail_f / R1_f + V_fail_f / R3_f
                R_eq_intact_f = 1 / (1/R1_f + 1/R2_f + 1/R3_f)
                V_intact_f = I_total_f * R_eq_intact_f
                I3_f = V_intact_f / R3_f
                
                valid_solutions.append({'R1': R1, 'R3': R3, 'I3': I3_f})

    if not valid_solutions:
        print("No solution found.")
        return

    # Find the solution with the maximum current
    best_solution = max(valid_solutions, key=lambda x: x['I3'])
    
    R1 = best_solution['R1']
    R3 = best_solution['R3']
    max_I3 = best_solution['I3']

    # Present the final calculation for the optimal case
    print(f"The optimal resistor values that satisfy all conditions and maximize the current are R1 = {R1} ohms and R3 = {R3} ohms.")
    print("\nThe calculation for the maximum current is as follows:")

    R1_f, R2_f, R3_f = Fraction(R1), Fraction(R2), Fraction(R3)
    V_fail_f = Fraction(V_fail)
    
    I_total_f = V_fail_f / R1_f + V_fail_f / R3_f
    R_eq_intact_f = 1 / (1/R1_f + 1/R2_f + 1/R3_f)
    V_intact_f = I_total_f * R_eq_intact_f
    I3_f = V_intact_f / R3_f

    print(f"\n1. Total source current (I_total) based on the failure condition (V = {V_fail}V):")
    print(f"   I_total = V / R1 + V / R3 = {V_fail_f} / {R1_f} + {V_fail_f} / {R3_f} = {I_total_f.numerator}/{I_total_f.denominator} A")

    print(f"\n2. Voltage (V_intact) when R2 = {R2} ohms is intact:")
    print(f"   R_eq = 1 / (1/{R1_f} + 1/{R2_f} + 1/{R3_f}) = {R_eq_intact_f.numerator}/{R_eq_intact_f.denominator} Ohms")
    print(f"   V_intact = I_total * R_eq = ({I_total_f.numerator}/{I_total_f.denominator}) * ({R_eq_intact_f.numerator}/{R_eq_intact_f.denominator}) = {V_intact_f.numerator}/{V_intact_f.denominator} V")

    print(f"\n3. Current through R3 (I3) when intact:")
    print(f"   I3 = V_intact / R3 = ({V_intact_f.numerator}/{V_intact_f.denominator}) / {R3_f} = {I3_f.numerator}/{I3_f.denominator} A")
    
    print(f"\nThe maximum possible current through R3 is {float(max_I3):.4f} A.")

solve_and_print_max_current()