import sys

def solve_bit_plus_plus():
    """
    This function solves the Bit++ problem by reading statements from standard input,
    calculating the final value of a variable 'x', and printing the result.
    
    The plan is as follows:
    1. Initialize a variable, x, to 0. This represents the single variable in Bit++.
    2. Read the first line of input to get the number of statements, n.
    3. Loop n times to process each statement.
    4. In each iteration, read a statement line.
    5. Check if the statement string contains a '+' character. This is a simple and
       effective way to identify an increment operation ('++X' or 'X++').
    6. If '+' is present, increment x by 1.
    7. If '+' is not present, the operation must be a decrement ('--X' or 'X--'),
       so decrement x by 1.
    8. After the loop completes, print the final integer value of x.
    """
    try:
        # Read the number of statements from the first line of input.
        num_statements = int(sys.stdin.readline())
        
        # Initialize the variable x to 0 as per the problem description.
        x = 0
        
        # Loop through each of the statements.
        for _ in range(num_statements):
            # Read the statement and remove any leading/trailing whitespace.
            statement = sys.stdin.readline().strip()
            
            # If a '+' is in the statement, it's an increment operation.
            if '+' in statement:
                x += 1
            # Otherwise, it must be a decrement operation.
            else:
                x -= 1
        
        # Print the final value of x.
        print(x)

    except (ValueError, IndexError):
        # Handle cases with malformed input.
        pass

solve_bit_plus_plus()