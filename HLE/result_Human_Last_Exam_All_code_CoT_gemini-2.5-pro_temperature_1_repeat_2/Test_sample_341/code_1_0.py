def solve_word_search():
    """
    Finds common metals in a grid, sorts them alphabetically,
    and prints the first 12.
    """
    # The grid of letters provided in the problem
    grid = [
        "NTITANIUMMIT", "ECDMCRASGATZ", "TOMTIIETONAI", "SBICTLSILGIN",
        "GAMANEZLDARC", "NLNISACIUNIG", "UTMIIDIIAEDM", "TLREVLISCSIC",
        "SICTAIRONEUO", "MPMUISENGAMP", "CMUNITALPNPP", "RCCGMNICKELE",
        "LRNMCADMIUMM", "OYRUCREMIMIL",
    ]

    # A list of common metals to search for
    metals_list = [
        "ALUMINUM", "ANTIMONY", "BERYLLIUM", "BISMUTH", "CADMIUM", "CALCIUM", "CERIUM",
        "CHROMIUM", "COBALT", "COPPER", "GOLD", "INDIUM", "IRIDIUM", "IRON", "LEAD",
        "LITHIUM", "MAGNESium", "MANGANESE", "MERCURY", "MOLYBDENUM", "NICKEL", "NIOBIUM",
        "OSMIUM", "PALLADIUM", "PLATINUM", "PLUTONIUM", "POTASSIUM", "RHODIUM", "SILVER",
        "SODIUM", "STRONTIUM", "TANTALUM", "TECHNETIUM", "THALLIUM", "THORIUM", "TIN",
        "TITANIUM", "TUNGSTEN", "URANIUM", "VANADIUM", "ZINC", "ZIRCONIUM"
    ]

    # Use a set to store found metals to avoid duplicates
    found_metals = set()

    # Get grid dimensions
    rows = len(grid)
    cols = len(grid[0])

    # Create a list of all strings to search in (horizontal, vertical, and their reverses)
    search_space = []
    
    # Add horizontal lines and their reverses
    for row_str in grid:
        search_space.append(row_str)
        search_space.append(row_str[::-1])
        
    # Add vertical lines and their reverses
    for j in range(cols):
        col_str = "".join([grid[i][j] for i in range(rows)])
        search_space.append(col_str)
        search_space.append(col_str[::-1])

    # Search for each metal in the generated search space
    for metal in metals_list:
        for line in search_space:
            if metal.upper() in line:
                found_metals.add(metal.upper())
                break  # Move to the next metal once found

    # Sort the found metals alphabetically
    sorted_metals = sorted(list(found_metals))

    # Format the first 12 metals for the output
    result = ", ".join(sorted_metals[:12])
    print(result)

solve_word_search()