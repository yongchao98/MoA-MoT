import collections

def solve():
    """
    Finds all common metals in the grid, sorts them, and prints the first 12.
    """
    grid = [
        "NTITANIUMMIT",
        "ECDMCRASGATZ",
        "TOMTIIETONAI",
        "SBICTLSILGIN",
        "GAMANEZLDARC",
        "NLNISACIUNIG",
        "UTMIIDIAEDM",
        "TLREVLISCSIC",
        "SICTAIRONEUO",
        "MPMUISENGAMP",
        "CMUNITALPNPP",
        "RCCGMNICKELE",
        "LRNMCADMIUMR",
        "OYRUCREMIMIL",
    ]

    # A list of common metals to search for
    metals_to_find = [
        "ALUMINUM", "CADMIUM", "COBALT", "COPPER", "GOLD", "IRON",
        "LEAD", "MAGNESIUM", "MERCURY", "NICKEL", "PLATINUM", "SILVER",
        "TIN", "TITANIUM", "TUNGSTEN", "ZINC"
    ]

    found_metals = set()
    num_rows = len(grid)
    num_cols = len(grid[0])

    # Create strings for horizontal and vertical search
    rows = grid
    cols = ["".join([grid[r][c] for r in range(num_rows)]) for c in range(num_cols)]
    
    search_lines = rows + cols

    for metal in metals_to_find:
        reversed_metal = metal[::-1]
        for line in search_lines:
            if metal in line or reversed_metal in line:
                found_metals.add(metal.capitalize())
                break
    
    # Sort the found metals alphabetically
    sorted_metals = sorted(list(found_metals))

    # Get the first 12 metals
    result_list = sorted_metals[:12]
    
    # Format and print the result
    print(", ".join(result_list))

solve()
<<<Cadmium, Cobalt, Copper, Iron, Magnesium, Mercury, Nickel, Platinum, Silver, Tin, Titanium, Tungsten>>>