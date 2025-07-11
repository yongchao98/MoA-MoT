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
    Finds the resistor values that satisfy the problem's conditions and
    calculates the maximum possible current through R3.
    """
    optimal_r1 = None
    optimal_r3 = None

    # We need to find the smallest R1 >= 10 that provides a valid prime R3.
    # This pair will maximize the current I3.
    # We search in a reasonable range for R1.
    for r1 in range(10, 100):
        # R3 must be between R1 + 3 and 2*R1 - 7 (inclusive)
        for r3 in range(r1 + 3, 2 * r1 - 6):
            if is_prime(r3):
                optimal_r1 = r1
                optimal_r3 = r3
                break
        if optimal_r1 is not None:
            break

    if optimal_r1 is None:
        print("No solution found in the given range.")
        return

    r2 = 6

    # Calculate the maximum current I3 using the derived formula
    # I_3 = (26 * 6 * (R1+R3)) / (R3 * (6*R1+6*R3+R1*R3))
    numerator = 156 * (optimal_r1 + optimal_r3)
    denominator = optimal_r3 * (6 * optimal_r1 + 6 * optimal_r3 + optimal_r1 * optimal_r3)
    max_current = numerator / denominator

    print(f"The conditions are met for the resistor values:")
    print(f"R1 = {optimal_r1} ohms")
    print(f"R2 = {r2} ohms")
    print(f"R3 = {optimal_r3} ohms")
    print("\nCalculating the maximum possible current through R3:")
    print(f"I_3 = (156 * ({optimal_r1} + {optimal_r3})) / ({optimal_r3} * (6*{optimal_r1} + 6*{optimal_r3} + {optimal_r1}*{optimal_r3}))")
    print(f"I_3 = {numerator} / {denominator}")
    print(f"I_3 = {numerator/math.gcd(numerator, denominator):.0f} / {denominator/math.gcd(numerator, denominator):.0f} Amps")
    print(f"\nThe maximum possible current is: {max_current} A")

find_resistors_and_current()