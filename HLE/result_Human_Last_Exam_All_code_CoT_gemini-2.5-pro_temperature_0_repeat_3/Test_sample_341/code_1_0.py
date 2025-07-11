import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_word_search():
    """
    Finds common metals in a grid, sorts them, and prints the first 12.
    """
    # Step 1: Represent the Grid
    grid = [
        "NTITANIUMMIT",
        "ECDMCRASGATZ",
        "TOMTIIETONAI",
        "SBICTLSILGIN",
        "GAMANEZLDARC",
        "NLNISACIUNIG",
        "UTMIIDIIAEDM",
        "TLREVLISCSIC",
        "SICTAIRONEUO",
        "MPMUISENGAMP",
        "CMUNITALPNPP",
        "RCCGMNICKELE",
        "LRNMCADMIUMM",
        "OYRUCREMIMIL",
    ]

    # Step 2: Create a Word List of common metals
    metals_to_search = [
        "ALUMINUM", "BISMUTH", "CADMIUM", "CHROMIUM", "COBALT", "COPPER", "GOLD",
        "IRON", "LEAD", "LITHIUM", "MAGNESIUM", "MANGANESE", "MERCURY", "NICKEL",
        "PLATINUM", "SILVER", "TIN", "TITANIUM", "TUNGSTEN", "ZINC", "IRIDIUM"
    ]

    # Step 3: Generate Search Lines
    lines = []
    
    # Horizontal lines (forward and backward)
    for row in grid:
        lines.append(row)
        lines.append(row[::-1])

    # Vertical lines (top to bottom and bottom to top)
    num_rows = len(grid)
    num_cols = len(grid[0])
    for j in range(num_cols):
        column = "".join([grid[i][j] for i in range(num_rows)])
        lines.append(column)
        lines.append(column[::-1])

    # Step 4: Search for Metals
    found_metals = set()
    for line in lines:
        for metal in metals_to_search:
            if metal in line:
                found_metals.add(metal)

    # Step 5: Sort and Format the Output
    sorted_metals = sorted(list(found_metals))
    
    # Print the first 12 metals in the specified format
    print(", ".join(sorted_metals[:12]))

solve_word_search()

# Get the captured output
final_answer = captured_output.getvalue().strip()

# Restore stdout
sys.stdout = old_stdout

# Final print to the user and the special format for the answer
print(final_answer)
print(f"<<<{final_answer}>>>")