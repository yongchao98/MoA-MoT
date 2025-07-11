def solve_grid_puzzle():
    """
    Parses the grid to find all integer values, calculates their sum,
    and prints the equation and the final result.
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

    total_sum = 0
    numbers_to_sum = []

    # Iterate through each row in the grid
    for row in grid:
        # Iterate through each item in the row
        for item in row:
            # Check if the item is an integer
            if isinstance(item, int):
                numbers_to_sum.append(item)
                total_sum += item

    # Build the equation string
    equation = " + ".join(map(str, numbers_to_sum))
    
    # Print the full equation and the final sum
    print(f"The hidden question in the deciphered text likely asks for the sum of the numbers in the grid.")
    print(f"The equation is: {equation} = {total_sum}")
    print(f"The final sum is: {total_sum}")

solve_grid_puzzle()