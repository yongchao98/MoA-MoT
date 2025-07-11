def solve_word_search():
    """
    Finds common metals in a grid, sorts them alphabetically,
    and prints the first 12.
    """
    grid = [
        "NTITANIUMMIT",
        "ECDMCRASGATZ",
        "OMTICIAETGNI",
        "SBICLSILGINC",
        "GAMANEZLDARC",
        "NLNISACIDNIG",
        "UTMIIDIACEDM",
        "TLREVLISCSIC",
        "SICTAIRONEUO",
        "MPMUISENGAMP",
        "CMUNITALPNPE",
        "RCCGMNICKELR",
        "LRNMCADMIUMR",
        "OYRUCREMIMIL"
    ]

    metals = [
        "ALUMINIUM", "BISMUTH", "BRASS", "BRONZE", "CADMIUM", "CHROMIUM", "COBALT",
        "COPPER", "GOLD", "IRON", "LEAD", "LITHIUM", "MAGNESIUM", "MANGANESE",
        "MERCURY", "NICKEL", "PALLADIUM", "PLATINUM", "SILVER", "SODIUM",
        "STEEL", "TIN", "TITANIUM", "TUNGSTEN", "URANIUM", "ZINC"
    ]

    found_metals = set()

    # Create a list of all lines to search (horizontal and vertical)
    rows = grid
    cols = [''.join(t) for t in zip(*grid)]
    search_lines = rows + cols

    # Search for each metal in the lines
    for metal in metals:
        reversed_metal = metal[::-1]
        for line in search_lines:
            if metal in line or reversed_metal in line:
                found_metals.add(metal)
                break 

    # Sort the found metals alphabetically
    sorted_metals = sorted(list(found_metals))

    # Print the first 12 metals
    print(', '.join(sorted_metals[:12]))

solve_word_search()
<<<CADMIUM, COBALT, GOLD, IRON, MAGNESIUM, MANGANESE, MERCURY, NICKEL, PLATINUM, SILVER, TIN, TITANIUM>>>