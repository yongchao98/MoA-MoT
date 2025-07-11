def find_metals_in_grid():
    """
    Finds common metals in a predefined grid of letters.
    The words can be horizontal or vertical, in either direction.
    """
    
    # The grid of letters from the problem
    grid = [
        "NTITANIUMMIT", "ECDMCRASGATZ", "TOMTIIETONAI", "SBICTLSIILGIN",
        "GAMANEZLDARC", "NLNISACIUNIG", "UTMIIDIIAEDM", "TLREVLISCSIC",
        "SICTAIRONEUO", "MPMUISENGAMP", "CMUNITAPLPNPP", "RCCGMNICKELE",
        "LRNMCADMIUMM", "OYRUCREMIMIL"
    ]

    # A list of common metals to search for
    metals_to_search = [
        "ALUMINUM", "BISMUTH", "BRASS", "BRONZE", "CADMIUM", "CHROMIUM",
        "COBALT", "COPPER", "GOLD", "IRON", "LEAD", "MAGNESIUM",
        "MANGANESE", "MERCURY", "NICKEL", "PLATINUM", "SILVER", "STEEL",
        "TIN", "TITANIUM", "TUNGSTEN", "ZINC"
    ]

    # The rows are the grid itself
    rows = grid
    # Transpose the grid to get the columns as strings
    num_cols = len(grid[0])
    cols = ["".join([row[i] for row in grid]) for i in range(num_cols)]

    # Combine rows and columns into a single search space
    search_space = rows + cols
    
    found_metals = set()

    # Iterate through each metal and search for it in the grid
    for metal in metals_to_search:
        for line in search_space:
            # Check for the metal name forwards
            if metal in line:
                found_metals.add(metal.capitalize())
            # Check for the metal name backwards
            if metal in line[::-1]:
                found_metals.add(metal.capitalize())

    # Sort the found metals alphabetically
    sorted_metals = sorted(list(found_metals))

    # Select the first 12 metals from the sorted list
    result_list = sorted_metals[:12]

    # Join the list into a comma-separated string for the final output
    final_answer = ", ".join(result_list)
    
    print(final_answer)

find_metals_in_grid()
<<<Cadmium, Cobalt, Gold, Iron, Lead, Manganese, Mercury, Nickel, Platinum, Silver, Tin, Titanium>>>