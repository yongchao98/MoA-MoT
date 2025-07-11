def solve_grid_puzzle():
    """
    Calculates the sum of the numbers found in the grid.
    The decrypted message is interpreted as a directive to sum these numbers.
    """
    # The numbers present in the grid provided in the puzzle
    numbers_in_grid = [75, 45, 86, 46, 20, 88, 49, 85, 22, 30]

    # Calculate the sum
    total_sum = sum(numbers_in_grid)

    # Create the equation string with each number
    equation = " + ".join(map(str, numbers_in_grid))

    # Print the full equation and the final answer
    print(f"The calculation is: {equation} = {total_sum}")
    print(f"The final sum is: {total_sum}")

solve_grid_puzzle()