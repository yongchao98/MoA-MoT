def find_metals_in_grid():
    """
    Finds common metals in a grid, sorts them alphabetically,
    and prints the first 12.
    """
    # The grid of letters from the problem description
    grid = [
        ['N', 'T', 'I', 'T', 'A', 'N', 'I', 'U', 'M', 'M', 'I', 'T'],
        ['E', 'C', 'D', 'M', 'C', 'R', 'A', 'S', 'G', 'A', 'T', 'Z'],
        ['T', 'O', 'M', 'T', 'I', 'I', 'E', 'T', 'O', 'N', 'A', 'I'],
        ['S', 'B', 'I', 'C', 'T', 'L', 'S', 'I', 'L', 'G', 'I', 'N'],
        ['G', 'A', 'M', 'A', 'N', 'E', 'Z', 'L', 'D', 'A', 'R', 'C'],
        ['N', 'L', 'N', 'I', 'S', 'A', 'C', 'I', 'U', 'N', 'I', 'G'],
        ['U', 'T', 'M', 'I', 'I', 'I', 'D', 'I', 'A', 'E', 'D', 'M'],
        ['T', 'L', 'R', 'E', 'V', 'L', 'I', 'S', 'C', 'S', 'I', 'C'],
        ['S', 'I', 'C', 'T', 'A', 'I', 'R', 'O', 'N', 'E', 'U', 'O'],
        ['M', 'P', 'M', 'U', 'I', 'S', 'E', 'N', 'G', 'A', 'M', 'P'],
        ['C', 'M', 'U', 'N', 'I', 'T', 'A', 'L', 'P', 'N', 'P', 'P'],
        ['R', 'C', 'C', 'G', 'M', 'N', 'I', 'C', 'K', 'E', 'L', 'E'],
        ['L', 'R', 'N', 'M', 'C', 'A', 'D', 'M', 'I', 'U', 'M', 'M'],
        ['O', 'Y', 'R', 'U', 'C', 'R', 'E', 'M', 'I', 'M', 'I', 'L']
    ]

    # A list of common metals to search for
    metals_to_find = [
        "ALUMINUM", "ANTIMONY", "BERYLLIUM", "BISMUTH", "CADMIUM", "CALCIUM", 
        "CHROMIUM", "COBALT", "COPPER", "GOLD", "IRON", "LEAD", "LITHIUM", 
        "MAGNESIUM", "MANGANESE", "MERCURY", "MOLYBDENUM", "NICKEL", "NIOBIUM", 
        "PLATINUM", "POTASSIUM", "SILVER", "SODIUM", "STEEL", "TANTALUM", "TIN", 
        "TITANIUM", "TUNGSTEN", "VANADIUM", "ZINC", "ZIRCONIUM"
    ]

    found_metals = set()

    # Search horizontally (left-to-right and right-to-left)
    for row in grid:
        row_str = "".join(row)
        for metal in metals_to_find:
            if metal in row_str or metal in row_str[::-1]:
                found_metals.add(metal)

    # Search vertically (top-to-bottom and bottom-to-top)
    num_rows = len(grid)
    if num_rows == 0:
        return
    num_cols = len(grid[0])
    for c in range(num_cols):
        col_str = "".join([grid[r][c] for r in range(num_rows)])
        for metal in metals_to_find:
            if metal in col_str or metal in col_str[::-1]:
                found_metals.add(metal)

    # Sort the found metals alphabetically and take the first 12
    sorted_metals = sorted(list(found_metals))
    result = sorted_metals[:12]

    # Print the result in the required format
    print(", ".join(result))

find_metals_in_grid()