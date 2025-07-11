import sys

def solve():
    """
    This function implements a correct interpreter for the X++ language.
    It reads the number of statements, then iterates through them,
    adjusting a counter variable 'x' accordingly.
    """
    try:
        # Read the total number of statements.
        # The input() function reads the entire line, so we strip whitespace and convert to int.
        num_statements_str = sys.stdin.readline()
        if not num_statements_str.strip():
            # Handle empty input for number of statements if necessary
            print(0)
            return
            
        n = int(num_statements_str)
        
        # Initialize the variable X to 0.
        x = 0
        
        # Process each statement
        for _ in range(n):
            statement = sys.stdin.readline().strip()
            
            # Check if the statement is an increment or decrement
            if "++" in statement:
                x += 1
            elif "--" in statement:
                x -= 1
        
        # Print the final value of X
        print(x)

    except (IOError, ValueError) as e:
        # Handle potential errors like bad input format
        sys.stderr.write(f"An error occurred: {e}\n")
        # In a real scenario, you might want to exit with an error code
        # For this problem, we can assume input is always well-formed.
        print(0)

solve()