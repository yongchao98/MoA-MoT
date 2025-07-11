def represent_in_base(n, c):
    """
    Implements the IMAG(n, c) algorithm to find the representation of a number n
    in a given base c.

    Args:
        n (int): The number to convert (n >= 0).
        c (int): The base to convert to (c >= 2).
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 0 or c < 2:
        print("Error: Input must be integers where n >= 0 and c >= 2.")
        return

    original_n = n
    
    # Handle the edge case n = 0
    if n == 0:
        print("The base", c, "representation of 0 is: (0)")
        print("Check: 0 = 0 * 2^0")
        return

    digits = []
    
    # Step 1: Initialization
    x = n
    q = x // c
    a0 = x - q * c  # This is equivalent to x % c
    digits.append(a0)

    # Step 2: Loop while q > 0
    while q > 0:
        # Step 2.1: Update variables
        x = q
        q = x // c
        ai = x - q * c # This is equivalent to x % c
        digits.append(ai)

    # The algorithm produces digits from least significant to most significant.
    # We reverse the list to get the correct representation (a_i, a_i-1, ..., a_0).
    representation = digits[::-1]
    
    # Step 3: Return/Output the result
    print(f"The base {c} representation of {original_n} is: {tuple(representation)}")
    
    # As requested, output the numbers in the final equation to verify the result
    equation_parts = []
    for power, digit in enumerate(reversed(digits)):
        equation_parts.append(f"{digit} * {c}^{power}")
    
    equation_str = " + ".join(equation_parts)
    print(f"Check: {original_n} = {equation_str}")


# Example: Find the representation of n=25 in base c=3.
# The algorithm should produce (2, 2, 1) because 25 = 2*3^2 + 2*3^1 + 1*3^0 = 18 + 6 + 1.
represent_in_base(25, 3)

print("\n" + "="*20 + "\n")

# Another example: Represent n=100 in base c=2 (binary).
represent_in_base(100, 2)