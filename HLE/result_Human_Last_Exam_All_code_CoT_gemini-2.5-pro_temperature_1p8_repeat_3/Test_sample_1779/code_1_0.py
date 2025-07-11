import math

def check_language_membership(input_string):
    """
    This function acts as a decider for the language L.
    A string `w` is in L if its length is a perfect square.
    This algorithm always halts, proving that L is a decidable language.

    Args:
        input_string (str): The string to check for membership in L.
    """
    n = len(input_string)

    print(f"Analyzing string: \"{input_string}\"")
    print(f"The length of the string is: {n}")

    # A non-negative integer n is a perfect square if its integer 
    # square root, when squared, equals n itself.
    # We use math.isqrt for an efficient integer square root (Python 3.8+).
    # For older versions, int(n**0.5) can be used.
    is_perfect_square = False
    if n >= 0:
        sqrt_n = math.isqrt(n)
        if sqrt_n * sqrt_n == n:
            is_perfect_square = True
    
    if is_perfect_square:
        # As requested, outputting the equation for the perfect square.
        print(f"The equation is: {sqrt_n} * {sqrt_n} = {n}")
        print("Result: The length is a perfect square, so the string IS in the decidable language L.\n")
    else:
        print(f"Result: The length {n} is not a perfect square, so the string IS NOT in the decidable language L.\n")

# --- Demonstration with different strings ---

# A string whose length (9) is a perfect square.
check_language_membership("123456789")

# A string whose length (10) is not a perfect square.
check_language_membership("0123456789")

# An empty string, whose length (0) is a perfect square.
check_language_membership("")

# A string whose length (16) is a perfect square.
check_language_membership("sixteen_is_4x4_!")