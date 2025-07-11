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

def calculate_i3(r1, r3):
    """Calculates the current I3 through resistor R3."""
    r2 = 6
    v_prime = 26
    numerator = v_prime * r2 * (r1 + r3)
    denominator = (r1 * r3) + (r1 * r2) + (r2 * r3)
    if denominator == 0:
        return 0
    return numerator / denominator

def solve_resistor_problem():
    """Finds the maximum possible current through R3."""
    max_i3 = 0
    best_r1 = None
    best_r3 = None
    r2 = 6
    
    # Store valid configurations found
    solutions = []

    # Case A: 6 is the maximum value (R1 < 6, R3 < 6)
    for r3_candidate in range(1, r2):
        if not is_prime(r3_candidate):
            continue
        for r1_candidate in range(1, r2):
            if r1_candidate == r3_candidate or r1_candidate == r2 or r3_candidate == r2:
                continue
            
            # Check the final condition
            if r3_candidate - r1_candidate > 2:
                solutions.append((r1_candidate, r3_candidate))

    # Case B: 6 is the minimum value (R1 > 6, R3 > 6)
    # We search a reasonable range; the function for I3 decreases as R values increase,
    # so the max will be found with small R values.
    for r3_candidate in range(r2 + 1, 50): # Search for prime R3 up to 50
        if not is_prime(r3_candidate):
            continue
        for r1_candidate in range(r2 + 1, 50): # Search for R1 up to 50
            if r1_candidate == r3_candidate:
                continue
            
            # Check the final condition
            if r3_candidate - r1_candidate > 2:
                 solutions.append((r1_candidate, r3_candidate))
    
    if not solutions:
        print("No valid resistor combinations found.")
        return

    # Calculate I3 for all valid solutions and find the maximum
    for r1, r3 in solutions:
        current_i3 = calculate_i3(r1, r3)
        if current_i3 > max_i3:
            max_i3 = current_i3
            best_r1 = r1
            best_r3 = r3

    # Output the result
    print(f"The maximum current is achieved with the following resistor values:")
    print(f"R1 = {best_r1} ohms")
    print(f"R2 = {r2} ohms")
    print(f"R3 = {best_r3} ohms (which is a prime number)")
    print("\nThe formula for the current I3 is:")
    print("I3 = (156 * (R1 + R3)) / (R1*R3 + 6*R1 + 6*R3)")
    print("\nPlugging in the values:")
    num_calc = 156 * (best_r1 + best_r3)
    den_calc = (best_r1 * best_r3) + 6 * best_r1 + 6 * best_r3
    print(f"I3 = (156 * ({best_r1} + {best_r3})) / ({best_r1}*{best_r3} + 6*{best_r1} + 6*{best_r3})")
    print(f"I3 = {num_calc} / {den_calc}")
    print(f"\nMaximum possible current through R3: {max_i3} A")


solve_resistor_problem()
<<<22.829268292682927>>>