def solve_puzzle():
    """
    This function solves the pattern puzzle and prints the final equation.
    The pattern is defined as:
    - If b is the square of a, the result is a * b.
    - Otherwise, the result is a * b + (a - 2).
    """
    a = 6
    b = 12

    # Check if b is the square of a
    if b == a * a:
        result = a * b
    else:
        # Otherwise, apply the second rule
        result = a * b + (a - 2)

    # Print the final equation with each number
    print(f"{a} + {b} = {result}")

solve_puzzle()