import sys
import io

# Capture original stdout to restore it later
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve_puzzle():
    """
    Solves the multi-step puzzle by extracting and summing numbers from the grid
    based on a decrypted instruction.
    """
    
    # The grid provided in the puzzle
    grid = [
        ['▣', 75, '◧', '◩', '◫', 45, 86],
        ['◨', '◨', '◪', '◨', '▨', '◪', '◫'],
        ['▤', '◫', '▦', 46, '◩', 20, '▣'],
        ['▥', '▧', '◨', 88, '▤', '▦', '◩'],
        ['◧', '◫', '◪', '◪', '▨', '◧', '▦'],
        ['▥', '▤', '▨', '▥', 49, '◨', 85],
        ['▩', 22, '▣', '◧', 30, '▥', '▩']
    ]

    # As per the decrypted instructions, we need to sum all numbers in the grid.
    # The condition "indicated by one of the following unicode symbols" is met
    # by every number in the grid when adjacency (including diagonals) is considered.
    
    numbers_to_sum = []
    for row in grid:
        for cell in row:
            if isinstance(cell, int):
                numbers_to_sum.append(cell)
    
    # Calculate the total sum
    total_sum = sum(numbers_to_sum)
    
    # Format the equation string
    equation_str = " + ".join(map(str, numbers_to_sum))
    
    # Print the final equation and the result
    print(f"{equation_str} = {total_sum}")

# Execute the function
solve_puzzle()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)