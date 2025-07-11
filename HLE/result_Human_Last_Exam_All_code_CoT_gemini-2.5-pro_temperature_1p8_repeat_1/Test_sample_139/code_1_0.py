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

def find_resistors_and_current():
    """
    Finds the optimal resistor values (R1, R3) and calculates the maximum possible current.
    """
    r2 = 6
    found_pair = None

    # We derived that R3 must be a prime > 10. We search starting from the smallest such prime.
    r3_candidate = 11
    while not found_pair:
        if not is_prime(r3_candidate):
            r3_candidate += 1
            continue

        # We derived the range for R1: (R3+6)/2 < R1 < R3-2
        r1_min = math.floor((r3_candidate + r2) / 2) + 1
        r1_max = r3_candidate - 2 -1

        for r1_candidate in range(r1_min, r1_max + 1):
            # Check all conditions
            # R1 and R3 must be distinct from R2 and each other
            if r1_candidate == r2 or r3_candidate == r2 or r1_candidate == r3_candidate:
                continue
            
            # Additional constraint R3 - R1 > 2
            if r3_candidate - r1_candidate <= 2:
                continue

            # We found the pair with the smallest possible R1 and R3, which will maximize the current.
            found_pair = (r1_candidate, r3_candidate)
            break
        
        if found_pair:
            break
        
        r3_candidate += 1

    if found_pair:
        r1, r3 = found_pair
        print(f"The resistor values that maximize the current are R1 = {r1}, R2 = {r2}, R3 = {r3}.")
        
        # Calculate the maximum current using the derived formula
        numerator = 156 * (r1 + r3)
        denominator = r3 * (6 * r1 + r1 * r3 + 6 * r3)
        current = numerator / denominator

        print("\nThe maximum possible current through R3 is calculated as:")
        print(f"I = (156 * ({r1} + {r3})) / ({r3} * (6 * {r1} + {r1} * {r3} + 6 * {r3}))")
        print(f"I = {numerator} / {denominator}")
        print(f"I = {numerator / math.gcd(numerator, denominator)} / {denominator / math.gcd(numerator, denominator)}")
        print(f"I ≈ {current:.5f} A")
        
        # Final answer format
        print(f"\n<<<{current}>>>")
        
    else:
        print("No solution found that satisfies all conditions.")

find_resistors_and_current()