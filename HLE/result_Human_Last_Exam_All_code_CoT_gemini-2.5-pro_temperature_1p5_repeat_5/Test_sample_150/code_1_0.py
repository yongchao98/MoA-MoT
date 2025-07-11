def solve():
    """
    Solves the Bit++ problem.
    """
    # Read the number of statements
    try:
        n = int(input())
    except (ValueError, EOFError):
        # Handle empty input or non-integer input
        print(0)
        return

    # Initialize the variable x
    x = 0

    # Process each statement
    for _ in range(n):
        try:
            statement = input()
            # Check if the statement contains '++'
            if '+' in statement:
                x += 1
            # Check if the statement contains '--'
            elif '-' in statement:
                x -= 1
        except EOFError:
            # Handle unexpected end of input
            break
            
    # Print the final value of x
    print(f"{x}")

solve()