import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def find_resistors_and_current():
    """
    Finds the values of R1 and R3 that satisfy the problem's conditions and
    calculates the maximum possible current through R3.
    """
    optimal_pair = None
    
    # We are given R3 - R1 > 2, which means R3 > R1.
    # The condition z({R1, 6, R3}) = 6 implies 6 is either the minimum or maximum value.
    # Case 1: R1, R3 < 6. For 6 to be furthest from the mean, we need 6 + R1 >= 2*R3.
    # Given R3 must be a prime < 6 and R1 < R3-2, the only possibilities are (R1=1, R3=5)
    # and (R1=2, R3=5). Neither satisfies 6 + R1 >= 2*R3. Thus, this case has no solution.

    # Case 2: R1, R3 > 6. For 6 to be furthest from the mean, we need R3 <= 2*R1 - 6.
    # To maximize the current I_3, we should seek the smallest valid R1 and R3 values.
    # Let's search for the first valid pair (R1, R3).
    # We'll search a reasonable range for R1, starting from the minimum possible value.
    for r1 in range(7, 100):
        # Per conditions, R3 must be in the range [r1 + 3, 2*r1 - 6]
        # We add 1 to the end of range to make it inclusive.
        for r3 in range(r1 + 3, 2 * r1 - 6 + 1):
            if is_prime(r3):
                # We've found the first valid pair. This pair will maximize the current.
                optimal_pair = (r1, r3)
                break
        if optimal_pair:
            break

    if optimal_pair:
        r1, r3 = optimal_pair
        print(f"Found resistor values that satisfy all conditions:")
        print(f"R1 = {r1} ohms")
        print(f"R3 = {r3} ohms")
        print("-" * 20)
        
        # Calculate the maximum current I_3 using the derived formula:
        # I_3 = 156 * (R1 + R3) / (R3 * (6 * (R1 + R3) + R1 * R3))
        
        numerator = 156 * (r1 + r3)
        denominator = r3 * (6 * (r1 + r3) + r1 * r3)
        max_current = numerator / denominator

        print("Calculating the current I_3 through R3:")
        # The final answer requires outputting the equation with numbers
        print(f"I_3 = 156 * ({r1} + {r3}) / ({r3} * (6 * ({r1} + {r3}) + {r1} * {r3}))")
        print(f"I_3 = {numerator} / {denominator}")
        print(f"I_3 = {max_current:.5f} Amperes")
        
        return max_current
    else:
        print("No solution found within the search range.")
        return None

# Execute the function to find the answer
max_i3_value = find_resistors_and_current()
# The final answer in the required format
# We use f-string formatting to control the number of decimal places for a clean answer.
print(f'<<<{max_i3_value:.5f}>>>')
