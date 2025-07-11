import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_word_search():
    """
    Finds common metals in a grid, sorts them alphabetically, and prints the first 12.
    """
    grid = [
        ['N', 'T', 'I', 'T', 'A', 'N', 'I', 'U', 'M', 'M', 'I', 'T'],
        ['E', 'C', 'D', 'M', 'C', 'R', 'A', 'S', 'G', 'A', 'T', 'Z'],
        ['T', 'O', 'M', 'T', 'I', 'I', 'E', 'T', 'O', 'N', 'A', 'I'],
        ['S', 'B', 'I', 'C', 'T', 'L', 'S', 'I', 'L', 'G', 'I', 'N'],
        ['G', 'A', 'M', 'A', 'N', 'E', 'Z', 'L', 'D', 'A', 'R', 'C'],
        ['N', 'L', 'N', 'I', 'S', 'A', 'C', 'I', 'U', 'N', 'I', 'G'],
        ['U', 'T', 'M', 'I', 'I', 'I', 'D', 'I', 'A', 'E', 'D', 'M'],
        ['T', 'L', 'R', 'E', 'V', 'L', 'I', 'S', 'C', 'S', 'I', 'C'],
        ['S', 'I', 'C', 'T', 'A', 'I', 'R', 'O', 'N', 'E', 'U', 'O'],
        ['M', 'P', 'M', 'U', 'I', 'S', 'E', 'N', 'G', 'A', 'M', 'P'],
        ['C', 'M', 'U', 'N', 'I', 'T', 'A', 'L', 'P', 'N', 'P', 'P'],
        ['R', 'C', 'C', 'G', 'M', 'N', 'I', 'C', 'K', 'E', 'L', 'E'],
        ['L', 'R', 'N', 'M', 'C', 'A', 'D', 'M', 'I', 'U', 'M', 'M'],
        ['O', 'Y', 'R', 'U', 'C', 'R', 'E', 'M', 'I', 'M', 'I', 'L']
    ]

    metal_list = [
        'ALUMINUM', 'BISMUTH', 'CADMIUM', 'CHROMIUM', 'COBALT', 'COPPER',
        'GOLD', 'IRON', 'LEAD', 'LITHIUM', 'MAGNESIUM', 'MANGANESE',
        'MERCURY', 'NICKEL', 'PLATINUM', 'SILVER', 'SODIUM', 'TIN',
        'TITANIUM', 'TUNGSTEN', 'URANIUM', 'ZINC'
    ]

    num_rows = len(grid)
    num_cols = len(grid[0])

    # Generate horizontal and vertical lines from the grid
    lines = []
    # Horizontal lines
    for r in range(num_rows):
        lines.append("".join(grid[r]))
    # Vertical lines
    for c in range(num_cols):
        lines.append("".join(grid[r][c] for r in range(num_rows)))

    # Create a full search space including reversed lines
    search_space = lines + [line[::-1] for line in lines]

    # Find all metals in the search space
    found_metals = set()
    for metal in metal_list:
        for line in search_space:
            if metal in line:
                found_metals.add(metal)

    # Sort the found metals and get the first 12
    sorted_metals = sorted(list(found_metals))
    result = sorted_metals[:12]

    # Print the result in the required format
    print(", ".join(result))

solve_word_search()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the final output and the special formatted answer
print(output)
final_answer = output.strip()
print(f"<<<{final_answer}>>>")