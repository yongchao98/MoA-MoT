import collections

def solve():
    """
    Finds all common metals in a grid, sorts them alphabetically, and prints the first 12.
    """
    # 1. Define the grid from the problem description. Note that row 6 is longer.
    grid = [
        "NTITANIUMMIT",
        "ECDMCRASGATZ",
        "TOMTIIETONAI",
        "SBICTLSILGIN",
        "GAMANEZLDARC",
        "NLNISACIUNIG",
        "UTMIIIDIAEDM",
        "TLREVLISCSIC",
        "SICTAIRONEUO",
        "MPMUISENGAMP",
        "CMUNITALPNPP",
        "RCCGMNICKELE",
        "LRNMCADMIUMR",
        "OYRUCREMIMIL",
    ]

    # 2. Define the list of common metals to search for.
    metals = [
        "ALUMINUM", "BISMUTH", "BRASS", "BRONZE", "CADMIUM", "CHROMIUM", 
        "COBALT", "COPPER", "GOLD", "IRON", "LEAD", "LITHIUM", "MAGNESIUM", 
        "MANGANESE", "MERCURY", "NICKEL", "PALLADIUM", "PLATINUM", "SILVER", 
        "SODIUM", "STEEL", "TIN", "TITANIUM", "TUNGSTEN", "ZINC"
    ]

    found_metals = set()
    rows = len(grid)
    if rows == 0:
        return
    max_cols = max(len(r) for r in grid)

    # 3. Iterate through each metal and search for it in the grid.
    for metal in metals:
        word_len = len(metal)
        metal_rev = metal[::-1]

        # Search horizontally
        for r in range(rows):
            if word_len <= len(grid[r]):
                if metal in grid[r] or metal_rev in grid[r]:
                    found_metals.add(metal)
                    break # Found this metal, move to the next one
        if metal in found_metals:
            continue

        # Search vertically
        found_vertically = False
        for c in range(max_cols):
            if found_vertically:
                break
            for r in range(rows - word_len + 1):
                # Build forward vertical string
                v_word_fwd = ""
                # Build backward vertical string
                v_word_bwd = ""
                possible = True
                for i in range(word_len):
                    # Check forward
                    if c >= len(grid[r + i]):
                        possible = False
                        break
                    v_word_fwd += grid[r + i][c]
                    
                    # Check backward
                    if c >= len(grid[r + word_len - 1 - i]):
                        possible = False
                        break
                    v_word_bwd += grid[r + word_len - 1 - i][c]

                if possible and (v_word_fwd == metal or v_word_bwd == metal):
                    found_metals.add(metal)
                    found_vertically = True
                    break
    
    # 4. Sort the found metals and print the first 12.
    sorted_metals = sorted(list(found_metals))
    print(", ".join(sorted_metals[:12]))

solve()
<<<CADMIUM, COBALT, COPPER, GOLD, IRON, MAGNESIUM, MERCURY, NICKEL, PLATINUM, SILVER, TITANIUM, TUNGSTEN>>>