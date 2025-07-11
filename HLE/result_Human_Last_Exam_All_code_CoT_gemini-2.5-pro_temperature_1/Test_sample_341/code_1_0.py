def solve_word_search():
    """
    Finds common metals in a grid, sorts them alphabetically,
    and prints the first 12.
    """
    # 1. Represent the Grid
    grid = [
        "NTITANIUMMIT",
        "ECDMCRASGATZ",
        "TOMTIIETONAI",
        "SBICLSILGIN",
        "GAMANEZLDARC",
        "NLNISACIUNIG",
        "UTMIIDIIAEDM",
        "TLREVLISCSIC",
        "SICTAIRONEUO",
        "MPMUISENGAMP",
        "CMUNITALPNPP",
        "RCCGMNICKELE",
        "LRNMCADMIUME",
        "OYRUCREMIMIL",
    ]

    # 2. Create a Word List of common metals
    metals_to_search = [
        "ALUMINUM", "BISMUTH", "CADMIUM", "CHROMIUM", "COBALT", "COPPER",
        "GOLD", "IRON", "LEAD", "MAGNESIUM", "MANGANESE", "MERCURY",
        "NICKEL", "PALLADIUM", "PLATINUM", "SILVER", "TIN", "TITANIUM",
        "TUNGSTEN", "ZINC"
    ]

    # 3. Build Search Space
    # Horizontal lines are just the grid itself
    horizontal_lines = grid

    # Vertical lines
    num_rows = len(grid)
    num_cols = len(grid[0])
    vertical_lines = []
    for j in range(num_cols):
        column = "".join([grid[i][j] for i in range(num_rows)])
        vertical_lines.append(column)

    all_lines = horizontal_lines + vertical_lines
    
    # 4. Search for Words
    found_metals = set()
    for metal in metals_to_search:
        # Search forwards
        for line in all_lines:
            if metal in line:
                found_metals.add(metal)
        
        # Search backwards
        reversed_metal = metal[::-1]
        for line in all_lines:
            if reversed_metal in line:
                found_metals.add(metal)

    # 5. Collect and Sort
    sorted_metals = sorted(list(found_metals))

    # 6. Format the Output
    # Get the first 12 metals
    result_list = sorted_metals[:12]
    
    print(", ".join(result_list))

solve_word_search()