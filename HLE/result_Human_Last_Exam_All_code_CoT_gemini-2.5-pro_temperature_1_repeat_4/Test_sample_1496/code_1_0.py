def solve_grid_puzzle():
    """
    This function solves the final part of the puzzle by summing the numbers in the grid.
    """
    grid = [
        ['▣', 75, '◧', '◩', '◫', 45, 86],
        ['◨', '◨', '◪', '◨', '▨', '◪', '◫'],
        ['▤', '◫', '▦', 46, '◩', 20, '▣'],
        ['▥', '▧', '◨', 88, '▤', '▦', '◩'],
        ['◧', '◫', '◪', '◪', '▨', '◧', '▦'],
        ['▥', '▤', '▨', '▥', 49, '◨', 85],
        ['▩', 22, '▣', '◧', 30, '▥', '▩']
    ]

    numbers_to_sum = []
    # Iterate through each row in the grid
    for row in grid:
        # Iterate through each item in the row
        for item in row:
            # Check if the item is an integer
            if isinstance(item, int):
                numbers_to_sum.append(item)

    # Calculate the sum
    total_sum = sum(numbers_to_sum)

    # Create the equation string
    equation_str = " + ".join(map(str, numbers_to_sum))

    # Print the final equation and the result
    print(f"{equation_str} = {total_sum}")

solve_grid_puzzle()