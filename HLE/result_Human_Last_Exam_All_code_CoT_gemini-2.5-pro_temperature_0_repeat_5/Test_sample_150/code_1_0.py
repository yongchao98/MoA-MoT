import sys

def solve():
    """
    This function solves the Bit++ problem by reading statements from standard input
    and calculating the final value of the variable x.
    """
    try:
        # Read the number of statements. The input might have trailing whitespace.
        n_str = sys.stdin.readline()
        if not n_str:
            return
        n = int(n_str.strip())
        
        # Initialize the variable x to 0.
        x = 0
        
        # Loop n times to process each statement.
        for _ in range(n):
            # Read one statement from input.
            statement = sys.stdin.readline().strip()
            
            # A statement is an increment if it contains a '+'.
            # Otherwise, it must be a decrement as it will contain a '-'.
            # This correctly handles all four cases: "++X", "X++", "--X", "X--".
            if '+' in statement:
                x += 1
            else:
                x -= 1
                
        # Print the final value of x.
        print(x)

    except (IOError, ValueError) as e:
        # Handle potential input errors gracefully.
        # For this problem's constraints, this is unlikely to be triggered.
        pass

solve()