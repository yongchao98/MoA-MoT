def solve_sequence():
    """
    Solves for the next number in the sequence based on the hypothesis
    that the last two distinct numbers form a simple equation.
    """
    # The last part of the sequence suggests an equation.
    # The two main numbers appear to be 36 and 124.
    # Let's assume the equation is a + x = b.
    a = 36
    b = 124

    # Solve for the unknown next number, x.
    x = b - a

    # The prompt asks to output each number in the final equation.
    print(f"{a} + {x} = {b}")

solve_sequence()