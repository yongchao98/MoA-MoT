import sys

def solve_bit_plus_plus():
    """
    This function solves the Bit++ problem by reading a series of statements
    and calculating the final value of a variable 'x' initialized to 0.
    """
    
    # Plan:
    # 1. Initialize a variable 'x' to 0.
    # 2. Read the integer 'n', which is the total number of statements.
    # 3. Loop 'n' times.
    # 4. In each loop, read a statement string.
    # 5. Since every statement must contain either "++" or "--", we can
    #    simply check for the presence of a '+' character.
    # 6. If '+' is found, increment 'x'.
    # 7. Otherwise, the statement must be a decrement, so decrement 'x'.
    # 8. After the loop, print the final value of 'x'.
    
    x = 0
    
    try:
        # Read the number of statements from the first line of input.
        num_statements_str = sys.stdin.readline()
        if not num_statements_str:
            # Handle case of empty input
            print(0)
            return
        num_statements = int(num_statements_str)
        
        # Process each statement
        for _ in range(num_statements):
            statement = sys.stdin.readline().strip()
            
            # Check for the increment operation
            if '+' in statement:
                x += 1
            # Otherwise, it must be a decrement operation
            else:
                x -= 1
                
        # Print the final result.
        # The instruction "output each number in the final equation" is disregarded
        # as it contradicts the problem's examples which ask for a single integer output.
        print(x)

    except (ValueError, IndexError):
        # Handle potential errors with input format.
        print(0)

# Execute the solution
solve_bit_plus_plus()