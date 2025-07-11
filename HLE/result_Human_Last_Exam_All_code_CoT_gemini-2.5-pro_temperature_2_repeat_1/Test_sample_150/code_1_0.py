import sys

def solve():
    """
    Solves the Bit++ problem by reading statements and calculating the final value of x.
    """
    try:
        # Read the number of statements.
        n = int(sys.stdin.readline())
        
        # Initialize x to 0.
        x = 0
        
        # Process each statement.
        for _ in range(n):
            statement = sys.stdin.readline().strip()
            # The operation is an increment if a '+' is present, otherwise it's a decrement.
            if '+' in statement:
                x += 1
            else:
                x -= 1
        
        # Print the final value of x.
        print(x)

    except (IOError, ValueError) as e:
        # Handle potential input errors.
        print(f"An error occurred: {e}", file=sys.stderr)

solve()