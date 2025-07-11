import sys

def solve():
    """
    This function solves the Bit++ problem by reading from standard input.
    """
    try:
        # Read the number of statements
        n = int(sys.stdin.readline())
        
        # Initialize variables
        initial_value = 0
        x = initial_value
        increments = 0
        decrements = 0
        
        # Process each statement
        for _ in range(n):
            statement = sys.stdin.readline().strip()
            if '+' in statement:
                x += 1
                increments += 1
            else:
                x -= 1
                decrements -= 1 # Stored as a negative number for the equation
                
        # The problem asks for the final value of x.
        # The instructions for me, the AI, ask to print the equation.
        # The final equation is: initial_value + increments - decrements = final_value
        # Note: decrements are stored as a negative, so we print it with a '+'
        # e.g., 0 + 3 + (-2) = 1 is the same as 0 + 3 - 2 = 1
        
        # To match the format "0 + {increments} - {decrements}" we use abs()
        print(f"{initial_value} + {increments} - {abs(decrements)} = {x}")

    except (IOError, ValueError) as e:
        print(f"An error occurred: {e}", file=sys.stderr)

solve()