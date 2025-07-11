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

def find_resistors_and_calculate_current():
    """
    Finds the resistor values based on the problem constraints and calculates the maximum current.
    """
    found_pair = None
    
    # From our analysis, R3 must be a prime number > 10. We check primes in increasing order.
    r3_candidate = 11
    while not found_pair:
        if is_prime(r3_candidate):
            # The constraints lead to the inequality: (R3 + 6)/2 < R1 < R3 - 2
            r1_lower_bound = (r3_candidate + 6) / 2
            r1_upper_bound = r3_candidate - 2
            
            # To maximize current, we need the smallest valid R1.
            # We check integers for R1 in its valid range.
            for r1_candidate in range(math.ceil(r1_lower_bound), r1_upper_bound):
                # The first pair we find will have the smallest R3 and R1, maximizing I_3.
                if r1_candidate > r1_lower_bound: # Precise check for the bound
                    found_pair = (r1_candidate, r3_candidate)
                    break
        if found_pair:
            break
        r3_candidate += 1
        
    R1, R3 = found_pair
    R2 = 6

    print("Step 1: Determine the resistor values.")
    print("The constraints on the resistors are:")
    print(f"  - R1, R2, R3 are distinct integers.")
    print(f"  - R2 = {R2} ohms.")
    print(f"  - R3 is a prime number.")
    print(f"  - z(R1, {R2}, R3) = {R2}, which implies R1 and R3 are both greater than {R2}.")
    print(f"  - R3 - R1 > 2.")
    print("\nSolving these constraints reveals that the optimal values to maximize the current are the smallest possible integers that satisfy them.")
    print(f"The unique solution is: R1 = {R1} ohms, R3 = {R3} ohms.\n")

    print("--------------------------------------------------\n")
    
    print("Step 2: Calculate the maximum possible current through R3.")
    print("The formula for the current I_3 through resistor R3 when the circuit is intact is derived from the circuit properties:")
    print("I_3 = 156 * (R1 + R3) / (R3 * (6*R1 + R1*R3 + 6*R3))\n")

    print("Substituting the determined resistor values into the equation:")
    
    term1 = 6 * R1
    term2 = R1 * R3
    term3 = 6 * R3
    inner_sum = term1 + term2 + term3
    numerator = 156 * (R1 + R3)
    denominator = R3 * inner_sum
    
    # Simplify fraction for cleaner output
    common_divisor = math.gcd(numerator, denominator)
    simple_num = numerator // common_divisor
    simple_den = denominator // common_divisor
    
    print(f"I_3 = 156 * ({R1} + {R3}) / ({R3} * (6*{R1} + {R1}*{R3} + 6*{R3}))")
    print(f"I_3 = {156 * (R1 + R3)} / ({R3} * ({term1} + {term2} + {term3}))")
    print(f"I_3 = {numerator} / ({R3} * {inner_sum})")
    print(f"I_3 = {numerator} / {denominator}")
    print(f"I_3 = {simple_num}/{simple_den} Amperes\n")

    final_current = numerator / denominator
    print(f"The maximum possible current through R3 is approximately: {final_current:.5f} A")


find_resistors_and_calculate_current()
<<<1.02985>>>