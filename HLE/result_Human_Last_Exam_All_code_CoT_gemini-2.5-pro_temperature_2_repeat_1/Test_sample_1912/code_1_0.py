def solve_puzzle():
    """
    Calculates the result for the final equation based on the discovered pattern.
    The pattern for 'a + b' is a * b + a - (b / a).
    """
    # The numbers for the final equation
    a = 6
    b = 12

    # Apply the discovered pattern
    # a * b + a - (b / a)
    result = a * b + a - (b / a)

    # The result should be an integer
    result = int(result)

    # Print the final equation with the calculated result
    print(f"{a} + {b} = {result}")

solve_puzzle()