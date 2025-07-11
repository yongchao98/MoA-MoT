import sys

def solve_bitpp():
    """
    Solves the Bit++ problem by reading statements from standard input,
    executing them, and printing the final value of the variable x.
    """
    try:
        # Read the number of statements from the first line.
        num_statements = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle cases where input is empty or not a valid number.
        print(0)
        return

    x = 0
    for _ in range(num_statements):
        # Read each statement. .strip() removes leading/trailing whitespace, including the newline.
        statement = sys.stdin.readline().strip()
        
        # The problem guarantees valid statements of 3 characters.
        # Checking the character at index 1 is sufficient and efficient to determine the operation.
        # If it's '+', it's an increment operation ('++X' or 'X++').
        # If it's '-', it's a decrement operation ('--X' or 'X--').
        if statement and statement[1] == '+':
            x += 1
        else:
            x -= 1
            
    # Print the final value of x.
    print(x)

solve_bitpp()