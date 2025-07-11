import sys

def solve_bit_plus_plus():
    """
    This function solves the Bit++ problem by reading a series of statements
    and calculating the final value of a variable 'x' initialized to 0.
    """
    try:
        # Read the number of statements from the first line of input.
        num_statements = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # If input is empty or invalid, assume 0 statements.
        num_statements = 0

    final_value_x = 0
    # Process each statement.
    for _ in range(num_statements):
        statement = sys.stdin.readline()
        
        # An increment operation must contain the '+' character.
        if '+' in statement:
            final_value_x += 1
        # A decrement operation must contain the '-' character.
        else:
            final_value_x -= 1
            
    # Print the final calculated value of x.
    print(final_value_x)

solve_bit_plus_plus()