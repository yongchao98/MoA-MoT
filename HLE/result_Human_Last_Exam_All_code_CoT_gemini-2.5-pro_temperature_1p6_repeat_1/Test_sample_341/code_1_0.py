def solve_word_search():
    """
    Finds common metals in a grid, sorts them, and prints the first 12.
    """
    # 1. Represent the Grid
    grid = [
        "NTITANIUMMIT", "ECDMCRASGATZ", "TOMTIIETONAI", "SBICTLSILGIN",
        "GAMANEZLDARC", "NLNISACIUNIG", "UTMIIDIIAEDM", "TLREVLISCSIC",
        "SICTAIRONEUO", "MPMUISENGAMP", "CMUNITALPNPP", "RCCGMNICKELE",
        "LRNMCADMIUMR", "OYRUCREMIMIL"
    ]

    # 2. Create a Word List of common metals
    # Metals with 'H' (BISMUTH, CHROMIUM) or 'LL' (PALLADIUM) are excluded
    # as those letters/combinations are not in the grid.
    metals_to_find = [
        "ALUMINUM", "CADMIUM", "COBALT", "COPPER", "GOLD", "IRON", "LEAD",
        "MAGNESIUM", "MANGANESE", "MERCURY", "NICKEL", "PLATINUM", "SILVER",
        "TIN", "TITANIUM", "TUNGSTEN", "ZINC"
    ]

    found_metals = set()
    num_rows = len(grid)
    num_cols = len(grid[0])

    # 3. Implement Search Logic
    # Generate horizontal and vertical strings from the grid
    horizontal_lines = grid
    vertical_lines = [''.join([grid[r][c] for r in range(num_rows)]) for c in range(num_cols)]
    
    all_lines = horizontal_lines + vertical_lines

    for line in all_lines:
        reversed_line = line[::-1]
        for metal in metals_to_find:
            if metal in line or metal in reversed_line:
                found_metals.add(metal)

    # 4. Collect and Sort Results
    sorted_metals = sorted(list(found_metals))

    # 5. Format and Print the first 12 metals
    result = sorted_metals[:12]
    print(", ".join(result))

solve_word_search()