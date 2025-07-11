import sys

def solve():
    """
    Solves the Bit++ problem by reading statements and calculating the final value of x.
    """
    try:
        # Read the number of statements from the first line.
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle cases where the first line is empty or not a valid integer.
        n = 0

    x = 0
    # Loop for each statement.
    for _ in range(n):
        # Read the statement from a new line.
        statement = sys.stdin.readline()
        
        # The most efficient logic is to check if the increment operator '+' is present.
        # If '+' is in the statement, it's an increment operation.
        # Otherwise, it must be a decrement operation.
        if '+' in statement:
            x += 1
        else:
            x -= 1
            
    # Print the final calculated value of x.
    print(x)

solve()