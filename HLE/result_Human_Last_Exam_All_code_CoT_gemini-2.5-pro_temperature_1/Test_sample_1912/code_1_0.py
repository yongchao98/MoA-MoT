def solve_sequence():
    """
    This function solves the sequence puzzle by applying the discovered pattern.
    The pattern for 'a + b' is a * b + a - (b / a).
    """
    # The numbers from the final line of the sequence
    a = 6
    b = 12

    # Apply the discovered pattern to find the result
    # We use integer division // as the pattern implies it.
    result = a * b + a - (b // a)

    # Print the final equation with the calculated result
    print(f"{a} + {b} = {result}")

solve_sequence()