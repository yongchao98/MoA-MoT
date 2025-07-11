import itertools

def solve_grid_puzzle():
    """
    Extracts all numbers from the provided grid, calculates their sum,
    and prints the full equation as requested.
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

    # Flatten the grid into a single list and filter for integers
    numbers = [item for item in itertools.chain(*grid) if isinstance(item, int)]
    
    # Calculate the sum
    total = sum(numbers)

    # Create the equation string as requested
    # The instruction says "output each number in the final equation"
    equation = " + ".join(map(str, numbers))
    
    # Print the equation and the final answer
    print(f"Equation: {equation} = {total}")
    print(f"The final answer is: {total}")

solve_grid_puzzle()