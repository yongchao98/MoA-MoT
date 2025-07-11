import sys

def rsa_multiply_and_print():
    """
    This function reads two large numbers as strings from standard input,
    multiplies them, and prints the result.
    This fulfills the "Output o" part of the task.
    """
    try:
        # Read p from the first line of input
        p_str = sys.stdin.readline().strip()
        
        # Read q from the second line of input
        q_str = sys.stdin.readline().strip()

        # The problem requires a program that can handle up to 100-digit numbers.
        # Python's int() automatically handles conversion of large number strings.
        if p_str and q_str:
            p = int(p_str)
            q = int(q_str)

            # Python's '*' operator works with arbitrarily large integers.
            o = p * q
            
            # Print the resulting number 'o'.
            print(o)

    except (ValueError, IndexError):
        # Handle cases with invalid or missing input.
        # As per the problem description, we can assume valid input.
        # No output will be printed in case of an error.
        pass

# Execute the function
rsa_multiply_and_print()