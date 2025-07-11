def solve_puzzle():
    """
    This function solves the puzzle based on the observed mathematical sequence.
    The pattern is deduced from the first three items:
    1st item resolves to 8.
    2nd item resolves to 5.
    3rd item resolves to 3.
    This implies an equation: 8 - 5 = 3.

    The pattern continues with the next logical equation in the sequence.
    """
    # The numbers from the sequence
    num1 = 5
    num2 = 3
    
    # The next step in the sequence
    result = num1 - num2
    
    # The final decoded word corresponds to the result '2', which is 'TWO'.
    # The prompt asks for each number in the final equation to be printed.
    
    print(f"{num1} - {num2} = {result}")

solve_puzzle()