import sys

def find_metals_in_grid():
    """
    Finds common metals in a grid of letters, sorts them, and prints the first 12.
    """
    # Step 1: Represent the Grid
    # The grid has 14 rows and 12 columns. Row 4 appears to be missing a character.
    # By aligning with the last column, which spells 'T-Z-I-N...', the missing character is 'N'.
    grid = [
        "NTITANIUMMIT",
        "ECDMCRASGATZ",
        "TOMTIIETONAI",
        "SBICTLSLGINN", # Original has 11 chars, 'N' is deduced from vertical alignment.
        "GAMANEZLDARC",
        "NLNISACIUNIG",
        "UTMIIIDIIAEM",
        "TLREVLISCSIC",
        "SICTAIRONEUO",
        "MPMUISENGAMP",
        "CMUNITALPNPP",
        "RCCGMNICKELE",
        "LRNMCADMIUMM",
        "OYRUCREMIMIL"
    ]

    # Step 2: Create a Word List
    metals = [
        "ALUMINIUM", "BISMUTH", "CADMIUM", "CHROMIUM",
        "COBALT", "COPPER", "GOLD", "IRON", "LEAD", "MAGNESIUM",
        "MANGANESE", "MERCURY", "NICKEL", "PLATINUM", "SILVER",
        "TIN", "TITANIUM", "TUNGSTEN", "ZINC"
    ]

    # Step 3: Generate Search Strings
    all_lines = []
    
    # Horizontal lines (forward and backward)
    horizontal_lines = grid
    reversed_horizontal_lines = [line[::-1] for line in horizontal_lines]
    all_lines.extend(horizontal_lines)
    all_lines.extend(reversed_horizontal_lines)

    # Vertical lines (forward and backward)
    rows = len(grid)
    cols = len(grid[0])
    for j in range(cols):
        column = "".join([grid[i][j] for i in range(rows)])
        all_lines.append(column)
        all_lines.append(column[::-1])

    # Step 4: Search for Metals
    found_metals = set()
    for metal in metals:
        for line in all_lines:
            if metal in line:
                found_metals.add(metal)

    # Step 5: Sort and Format
    sorted_metals = sorted(list(found_metals))
    
    # Select the first 12 metals
    result_list = sorted_metals[:12]
    
    # Print the result in the required format
    print(", ".join(result_list))

find_metals_in_grid()