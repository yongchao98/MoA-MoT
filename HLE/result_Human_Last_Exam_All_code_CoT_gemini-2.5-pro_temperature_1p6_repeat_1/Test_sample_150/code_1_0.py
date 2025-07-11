import sys

def solve():
    """
    Solves the Bit++ problem.
    """
    # Read the number of statements from the first line of input.
    try:
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle cases with no input or invalid input
        print(0)
        return

    # Initialize the variable x to 0.
    x = 0

    # Process each of the n statements.
    for _ in range(n):
        # Read the statement from a new line.
        statement = sys.stdin.readline().strip()

        # Check if the statement contains an increment operation.
        # The '+' character is unique to increment operations (++X or X++).
        if '+' in statement:
            x += 1
        # Otherwise, it must contain a decrement operation.
        else:
            x -= 1

    # Print the final value of x.
    print(x)

solve()