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

def solve_circuit():
    """
    Finds the maximum possible current through R3 by searching for valid
    resistor values R1 and R3 that satisfy all problem constraints.
    """
    R2 = 6
    max_I3 = 0.0
    best_R1 = None
    best_R3 = None

    # Search for R3 (prime number). As I3 generally decreases with higher R,
    # a limited search range is sufficient.
    for r3_candidate in range(2, 100):
        if not is_prime(r3_candidate):
            continue
        R3 = r3_candidate

        # Search for R1 based on the constraint R3 - R1 > 2 => R1 < R3 - 2
        for R1 in range(1, R3 - 2):
            
            # Constraint: R1, R2, R3 must be distinct integers
            # R1 and R3 are not equal to R2=6
            if R1 == R2 or R3 == R2:
                continue
            
            # Constraint: R1 and R3 are distinct, guaranteed by R1 < R3 - 2

            # Constraint: z(R1, R2, R3) = 6 means 6 is the min or max value.
            is_z_ok = (R2 < R1 and R2 < R3) or (R2 > R1 and R2 > R3)
            if not is_z_ok:
                continue

            # If all constraints are met, calculate I3
            # V = 156 * (R1 + R3) / (R1*R3 + 6*(R1 + R3))
            # I3 = V / R3
            numerator_v = 156 * (R1 + R3)
            denominator_v = R1 * R3 + R2 * (R1 + R3)
            voltage = numerator_v / denominator_v
            
            current_I3 = voltage / R3
            
            # Check if this is the maximum current found so far
            if current_I3 > max_I3:
                max_I3 = current_I3
                best_R1 = R1
                best_R3 = R3

    print("To find the maximum possible current through R3, we evaluate the conditions.")
    print(f"The optimal combination of resistors found is R1 = {best_R1} ohms, R2 = {R2} ohms, and R3 = {best_R3} ohms.")
    print("These values satisfy all conditions:")
    print(f"- R3 ({best_R3}) is a prime number.")
    print(f"- R3 - R1 = {best_R3} - {best_R1} = {best_R3 - best_R1}, which is greater than 2.")
    print(f"- R2=6 is the maximum value in the set {{{best_R1}, {R2}, {best_R3}}}, satisfying the z(C)=6 condition.")
    print("\nNow, we calculate the current through R3 (I3) for this configuration.")
    print("The formula for the voltage V in the intact circuit is:")
    print("V = 156 * (R1 + R3) / (R1*R3 + 6*(R1 + R3))")

    num_v = 156 * (best_R1 + best_R3)
    den_v = best_R1 * best_R3 + R2 * (best_R1 + best_R3)
    
    print("\nStep-by-step calculation:")
    print(f"V = 156 * ({best_R1} + {best_R3}) / ({best_R1}*{best_R3} + 6*({best_R1} + {best_R3}))")
    print(f"V = 156 * ({best_R1 + best_R3}) / ({best_R1 * best_R3} + 6*({best_R1 + best_R3}))")
    print(f"V = {num_v} / ({best_R1 * best_R3} + {6 * (best_R1 + best_R3)})")
    print(f"V = {num_v} / {den_v}")
    
    voltage = num_v / den_v
    print(f"V = {voltage} volts")
    
    print("\nThe current I3 is V / R3:")
    print(f"I3 = {voltage} / {best_R3}")
    print(f"I3 = {max_I3} Amperes")

solve_circuit()
print(f'<<<{936/205}>>>')