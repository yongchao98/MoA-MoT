import sys

def solve():
    """
    Solves the Bit++ problem by reading statements and calculating the final value of x.
    """
    # Read the number of statements from input.
    # The first line contains a single integer n.
    try:
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle cases where input is empty or not a valid integer
        print(0)
        return

    # The initial value of x is 0.
    x = 0

    # Process each of the n statements.
    for _ in range(n):
        # Read the statement from the next line.
        statement = sys.stdin.readline().strip()
        
        # A statement contains exactly one operation (++ or --).
        # We can check for the '+' character to determine the operation.
        if '+' in statement:
            x += 1
        else:
            x -= 1
            
    # Print the final value of x.
    print(x)

solve()