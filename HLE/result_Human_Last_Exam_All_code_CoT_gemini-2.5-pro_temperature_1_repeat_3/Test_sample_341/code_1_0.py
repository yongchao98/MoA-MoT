import collections

def solve():
    """
    Finds all common metals in a grid, sorts them alphabetically, and prints the first 12.
    """
    # The grid of letters, with empty columns removed.
    grid = [
        "NTITANIUMMIT",
        "ECDMCRASGATZ",
        "TOMTIIETONAI",
        "SBICTLSLGIN",
        "GAMANEZLDARC",
        "NLNISACIUNIG",
        "UTMIIDIIAEDM",
        "TLREVLISCSIC",
        "SICTAIRONEUO",
        "MPMUISENGAMP",
        "CMUNITALPNPP",
        "RCCGMNICKELE",
        "LRNMCADMIUMM",
        "OYRUCREMIMIL"
    ]

    # A list of common metals to search for.
    metals = [
        "ALUMINUM", "BISMUTH", "CADMIUM", "CHROMIUM", "COBALT", "COPPER",
        "GOLD", "IRON", "LEAD", "MAGNESIUM", "MANGANESE", "MERCURY",
        "NICKEL", "PLATINUM", "SILVER", "TIN", "TITANIUM", "TUNGSTEN", "ZINC"
    ]

    # Create strings for rows (horizontal).
    rows = grid

    # Create strings for columns (vertical).
    num_cols = len(grid[0])
    cols = []
    for j in range(num_cols):
        col_str = "".join([row[j] for row in grid])
        cols.append(col_str)

    # Combine all searchable strings.
    search_strings = rows + cols

    found_metals = set()

    # Search for each metal in all directions.
    for metal in metals:
        # Also check for the reversed word.
        metal_rev = metal[::-1]
        for s in search_strings:
            if metal in s or metal_rev in s:
                found_metals.add(metal)
                # Found it, move to the next metal.
                break

    # Sort the found metals alphabetically.
    sorted_metals = sorted(list(found_metals))
    
    # Get the first 12 metals.
    result = sorted_metals[:12]

    # Print the result as a comma-separated string.
    print(", ".join(result))

solve()