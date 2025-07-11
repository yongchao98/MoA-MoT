import sys

def solve():
    """
    Solves the Bit++ problem by reading statements from standard input
    and calculating the final value of x.
    """
    try:
        # Read the number of statements.
        num_statements = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # If input is empty or not a number, assume 0 statements.
        num_statements = 0

    x = 0
    for _ in range(num_statements):
        statement = sys.stdin.readline().strip()

        # The operation is uniquely determined by the character at index 1.
        # If it's '+', it's an increment operation.
        # If it's '-', it's a decrement operation.
        if '+' in statement:
            x += 1
        else:
            x -= 1

    # Print the final result.
    print(x)

solve()