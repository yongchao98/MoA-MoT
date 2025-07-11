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

def solve_circuit_problem():
    """
    Finds resistor values and calculates the maximum possible current through R3.
    """
    R1, R3 = 0, 0
    R2 = 6
    found = False

    # Find the smallest integer R1 > 8 for which a prime R3 exists
    # in the range (R1 + 2, 2*R1 - 6). This will maximize the current.
    for r1_candidate in range(9, 200):  # Search R1 starting from 9
        lower_bound_r3 = r1_candidate + 2
        upper_bound_r3 = 2 * r1_candidate - 6
        for r3_candidate in range(lower_bound_r3 + 1, upper_bound_r3):
            if is_prime(r3_candidate):
                R1 = r1_candidate
                R3 = r3_candidate
                found = True
                break
        if found:
            break

    if not found:
        print("Could not find suitable resistor values.")
        return

    # Scenario 1: All resistors in parallel.
    # When R2 fails, V_parallel = 26V. The source current is I = 26 * (1/R1 + 1/R3).
    # With R2 intact, V_intact = I / (1/R1 + 1/R2 + 1/R3).
    # Current through R3 is I3 = V_intact / R3.
    I_source_1 = 26 * (1/R1 + 1/R3)
    R_eq_inv_1 = (1/R1 + 1/R2 + 1/R3)
    V_intact_1 = I_source_1 / R_eq_inv_1
    I3_scenario1 = V_intact_1 / R3

    # Scenario 2: R1 in series with a parallel combination of R2 and R3.
    # When R2 fails, V_R3 = I * R3 = 26V, so source current I = 26 / R3.
    # With R2 intact, current through R3 is I3 = I * (R2 / (R2 + R3)).
    I_source_2 = 26 / R3
    I3_scenario2 = I_source_2 * (R2 / (R2 + R3))

    # Determine the maximum possible current.
    max_I3 = max(I3_scenario1, I3_scenario2)
    
    print(f"The optimal resistor values that satisfy all conditions are:")
    print(f"R1 = {R1} Ω, R2 = {R2} Ω, R3 = {R3} Ω (prime)")

    print("\nOf the two possible circuit configurations, the parallel circuit yields the maximum current.")
    print("The calculation for this maximum current is:")
    
    # Printing the equation with the numbers plugged in.
    print(f"I_3 = ( (26 * (1/{R1} + 1/{R3})) / (1/{R1} + 1/{R2} + 1/{R3}) ) / {R3}")
    
    # Calculate and print intermediate fraction for clarity
    num = 26 * (R1 + R3) * R2
    den = R3 * (R1*R2 + R1*R3 + R2*R3)
    common_divisor = math.gcd(num, den)
    
    print(f"I_3 = {num//common_divisor} / {den//common_divisor} A")
    print(f"\nThe maximum possible current through R3 is {max_I3:.4f} Amperes.")

solve_circuit_problem()
<<<1.0299>>>