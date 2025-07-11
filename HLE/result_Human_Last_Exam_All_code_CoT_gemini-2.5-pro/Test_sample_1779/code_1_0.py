import math

def decide_language_L(w: str):
    """
    This function decides if a string 'w' belongs to the language L,
    where L = {w | len(w) is a perfect square}.
    It simulates the decider algorithm described.
    """
    print(f"--- Checking input string w = '{w}' ---")
    n = len(w)
    print(f"Step 1: The length of w is n = {n}.")

    if n < 0: # Should not happen for string length, but good practice
        print(f"Step 2: Length {n} is negative, so it cannot be a perfect square.")
        print(f"Result: The string is NOT in language L.")
        return

    print(f"Step 2: Checking if n = {n} is a perfect square.")
    
    # An efficient way to check is to take the integer part of the square root
    # and see if squaring it gives back the original number.
    # This is equivalent to the iterative checking described in the plan.
    
    if n == 0:
        print("Step 3: 0 is a perfect square (0 * 0 = 0).")
        print(f"Result: The string IS in language L.")
        return

    # Calculate integer part of the square root
    sqrt_n = int(math.sqrt(n))
    
    # Check if squaring this integer gives n
    square_of_sqrt_n = sqrt_n * sqrt_n

    print(f"Step 3: The integer part of the square root of {n} is {sqrt_n}.")
    print(f"Step 4: Checking if {sqrt_n} * {sqrt_n} = {n}.")
    print(f"   Calculation: {sqrt_n} * {sqrt_n} = {square_of_sqrt_n}")
    
    if square_of_sqrt_n == n:
        print(f"Result: Since the equation holds, {n} is a perfect square.")
        print(f"The string IS in language L.")
    else:
        print(f"Result: Since {square_of_sqrt_n} != {n}, {n} is not a perfect square.")
        print(f"The string IS NOT in language L.")

# Let's test the decider with a few example strings
test_strings = ["", "####", "#########", "hello", "decidable"]
for s in test_strings:
    decide_language_L(s)
    print()
