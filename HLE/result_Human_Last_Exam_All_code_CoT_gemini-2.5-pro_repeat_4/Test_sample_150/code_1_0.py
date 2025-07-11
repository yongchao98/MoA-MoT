import sys

def solve_bit_plus_plus():
    """
    Solves the Bit++ problem by reading statements and calculating the final value of x.
    """
    try:
        # Read the number of statements from the first line of input.
        num_statements = int(sys.stdin.readline())
    except (IOError, ValueError):
        # Handle cases with invalid or no input for the number of statements.
        num_statements = 0

    # The initial value of the variable x is 0.
    x = 0

    # Process each statement.
    for _ in range(num_statements):
        # Read the statement from a new line and remove any surrounding whitespace.
        statement = sys.stdin.readline().strip()

        # An efficient way to check the operation is to look at the middle character.
        # For "X++", "++X", "X--", and "--X", the character at index 1 is always
        # '+' for an increment and '-' for a decrement.
        if statement and statement[1] == '+':
            x += 1
        elif statement:
            x -= 1
    
    # Print the final value of x.
    print(x)

# Execute the solution function.
solve_bit_plus_plus()