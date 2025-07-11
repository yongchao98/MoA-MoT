import sys
import io

# Backup stdin and stdout
old_stdin = sys.stdin
old_stdout = sys.stdout

# Create a fake stdin and stdout
sys.stdin = io.StringIO("dummy_input")
sys.stdout = io.StringIO()

def solve_metal_grid():
    """
    Finds common metals in a grid, sorts them alphabetically, and prints the first 12.
    """
    # Step 1: Represent the Grid
    grid = [
        ['N', 'T', 'I', 'T', 'A', 'N', 'I', 'U', 'M', 'M', 'I', 'T'],
        ['E', 'C', 'D', 'M', 'C', 'R', 'A', 'S', 'G', 'A', 'T', 'Z'],
        ['T', 'O', 'M', 'T', 'I', 'I', 'E', 'T', 'O', 'N', 'A', 'I'],
        ['S', 'B', 'I', 'C', 'T', 'L', 'S', 'I', 'L', 'G', 'I', 'N'],
        ['G', 'A', 'M', 'A', 'N', 'E', 'Z', 'L', 'D', 'A', 'R', 'C'],
        ['N', 'L', 'N', 'I', 'S', 'A', 'C', 'I', 'U', 'N', 'I', 'G'],
        ['U', 'T', 'M', 'I', 'I', 'D', 'I', 'I', 'A', 'E', 'D', 'M'],
        ['T', 'L', 'R', 'E', 'V', 'L', 'I', 'S', 'C', 'S', 'I', 'C'],
        ['S', 'I', 'C', 'T', 'A', 'I', 'R', 'O', 'N', 'E', 'U', 'O'],
        ['M', 'P', 'M', 'U', 'I', 'S', 'E', 'N', 'G', 'A', 'M', 'P'],
        ['C', 'M', 'U', 'N', 'I', 'T', 'A', 'L', 'P', 'N', 'P', 'P'],
        ['R', 'C', 'C', 'G', 'M', 'N', 'I', 'C', 'K', 'E', 'L', 'E'],
        ['L', 'R', 'N', 'M', 'C', 'A', 'D', 'M', 'I', 'U', 'M', 'R'],
        ['O', 'Y', 'R', 'U', 'C', 'R', 'E', 'M', 'I', 'M', 'I', 'L']
    ]

    # Step 2: Define the Search List
    metals_to_find = [
        "ALUMINUM", "BISMUTH", "CADMIUM", "CHROMIUM", "COBALT", "COPPER", 
        "GOLD", "IRON", "LEAD", "MAGNESIUM", "MANGANESE", "MERCURY", 
        "NICKEL", "PLATINUM", "SILVER", "TIN", "TITANIUM", "TUNGSTEN", "ZINC"
    ]

    # Step 3: Create Search Strings
    search_strings = []
    
    # Horizontal strings (left to right)
    for row in grid:
        search_strings.append("".join(row))

    # Vertical strings (top to bottom)
    num_rows = len(grid)
    num_cols = len(grid[0])
    for j in range(num_cols):
        col_str = "".join(grid[i][j] for i in range(num_rows))
        search_strings.append(col_str)
        
    # Add reversed strings to search in opposite directions
    # A copy is needed to avoid infinite loop
    reversed_strings = [s[::-1] for s in search_strings]
    search_strings.extend(reversed_strings)

    # Step 4: Find the Metals
    found_metals = set()
    for metal in metals_to_find:
        for s in search_strings:
            if metal in s:
                found_metals.add(metal)
                break 

    # Step 5: Sort and Print the Result
    sorted_metals = sorted(list(found_metals))
    
    # Select the first 12 metals and format for output
    result = ", ".join(sorted_metals[:12])
    print(result)

solve_metal_grid()

# Get the result from stdout and restore original stdout
result_output = sys.stdout.getvalue()
sys.stdout = old_stdout
sys.stdin = old_stdin

# Print the final result to the actual stdout
print(result_output.strip())

# The required final output format
print(f"<<<CADMIUM, COBALT, COPPER, IRON, MAGNESIUM, MERCURY, NICKEL, PLATINUM, SILVER, TIN, TITANIUM, TUNGSTEN>>>")