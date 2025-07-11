def solve_bit_plus_plus():
    """
    Solves the Bit++ programming problem.
    The initial value of x is 0. The program reads a series of statements
    and calculates the final value of x.
    """
    try:
        # Read the number of statements.
        n = int(input())
    except (ValueError, EOFError):
        # Handle cases with no input or invalid input for n.
        print(0)
        return

    x = 0
    for _ in range(n):
        try:
            # Read the statement for the current line.
            statement = input()
            # Check for the operation type and update x accordingly.
            if '+' in statement:
                x += 1
            else:
                # The problem guarantees the statement will have either '++' or '--'.
                x -= 1
        except EOFError:
            # If input ends unexpectedly, break the loop.
            break
            
    # The problem asks to output the final value of x.
    # The instruction "output each number in the final equation" does not apply
    # to this problem's requirements, which ask for a single integer output.
    print(x)

solve_bit_plus_plus()