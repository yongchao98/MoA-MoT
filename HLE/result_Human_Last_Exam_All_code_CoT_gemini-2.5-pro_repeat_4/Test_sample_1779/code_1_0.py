import math

def decide_language_membership(input_string):
    """
    This function acts as a decider for the language L where a string is in L
    if and only if its length is a perfect square.
    This algorithm is guaranteed to halt for any input.
    """
    n = len(input_string)
    print(f"Analyzing input string: \"{input_string}\"")
    print(f"Length of the string is n = {n}")

    is_perfect_square = False
    # A negative number cannot be a perfect square of an integer.
    # Length is always non-negative.
    if n >= 0:
        # Calculate the integer part of the square root
        root = int(math.sqrt(n))
        # Check if squaring the root gives back the original number
        if root * root == n:
            is_perfect_square = True
            # Print the equation as requested
            print(f"Checking if n is a perfect square... Yes, it is.")
            print(f"The equation is: {root} * {root} = {n}")

    if is_perfect_square:
        print(f"Result: The string's length is a perfect square. The string is in language L.")
    else:
        print(f"Checking if n is a perfect square... No, it is not.")
        print(f"Result: The string's length is not a perfect square. The string is not in language L.")
    print("-" * 30)

# --- Main execution with example cases ---
# Test cases include lengths that are and are not perfect squares.
test_strings = [
    "123456789",  # length 9 (3*3)
    "abcd",       # length 4 (2*2)
    "a",          # length 1 (1*1)
    "",           # length 0 (0*0)
    "hello",      # length 5 (not a perfect square)
    "undecidable" # length 11 (not a perfect square)
]

for test_str in test_strings:
    decide_language_membership(test_str)
