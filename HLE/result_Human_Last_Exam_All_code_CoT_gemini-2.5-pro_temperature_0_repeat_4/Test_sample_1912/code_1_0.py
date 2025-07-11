def solve_sequence():
    """
    This function solves the sequence puzzle by finding the hidden pattern.
    The pattern is determined to be: a * b + (a - 2) if b == 2*a, otherwise a * b.
    It then applies this pattern to the final numbers in the sequence.
    """
    # The numbers for the final equation
    a = 6
    b = 12

    # The pattern is a * b + k, where k depends on the relationship between a and b.
    # If b is twice a, k is a - 2.
    # Otherwise, k is 0. This handles the special case of 3 + 9 = 27.
    if b == 2 * a:
        k = a - 2
    else:
        k = 0
    
    # Calculate the result
    result = a * b + k

    # Print the final equation as requested
    print(f"{a} + {b} = {result}")

solve_sequence()