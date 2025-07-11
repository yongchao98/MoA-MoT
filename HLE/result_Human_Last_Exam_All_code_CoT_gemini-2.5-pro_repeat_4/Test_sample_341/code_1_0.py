def solve_word_search():
    """
    Finds all common metals in a grid, sorts them alphabetically,
    and prints the first 12.
    """
    # Step 1: Represent the Grid
    grid = [
        ['N', 'T', 'I', 'T', 'A', 'N', 'I', 'U', 'M', 'M', 'I', 'T'],
        ['E', 'C', 'D', 'M', 'C', 'R', 'A', 'S', 'G', 'A', 'T', 'Z'],
        ['T', 'O', 'M', 'T', 'I', 'I', 'E', 'T', 'O', 'N', 'A', 'I'],
        ['S', 'B', 'I', 'C', 'T', 'L', 'S', 'I', 'L', 'G', 'I', 'N'],
        ['G', 'A', 'M', 'A', 'N', 'E', 'Z', 'L', 'D', 'A', 'R', 'C'],
        ['N', 'L', 'N', 'I', 'S', 'A', 'C', 'I', 'U', 'N', 'I', 'G'],
        ['U', 'T', 'M', 'I', 'I', 'D', 'I', 'I', 'A', 'E', 'D', 'M'],
        ['T', 'L', 'R', 'E', 'V', 'L', 'I', 'S', 'C', 'S', 'I', 'C'],
        ['S', 'I', 'C', 'T', 'A', 'I', 'R', 'O', 'N', 'E', 'U', 'O'],
        ['M', 'P', 'M', 'U', 'I', 'S', 'E', 'N', 'G', 'A', 'M', 'P'],
        ['C', 'M', 'U', 'N', 'I', 'T', 'A', 'L', 'P', 'N', 'P', 'P'],
        ['R', 'C', 'C', 'G', 'M', 'N', 'I', 'C', 'K', 'E', 'L', 'E'],
        ['L', 'R', 'N', 'M', 'C', 'A', 'D', 'M', 'I', 'U', 'M', 'R'],
        ['O', 'Y', 'R', 'U', 'C', 'R', 'E', 'M', 'I', 'M', 'I', 'L']
    ]

    # Step 2: Define Metal List
    metals_to_find = [
        "ALUMINUM", "BISMUTH", "CADMIUM", "CHROMIUM", "COBALT", "COPPER",
        "GOLD", "IRON", "LEAD", "MAGNESIUM", "MANGANESE", "MERCURY",
        "NICKEL", "PLATINUM", "SILVER", "TIN", "TITANIUM", "TUNGSTEN", "ZINC"
    ]

    # Step 3: Generate Search Strings
    num_rows = len(grid)
    num_cols = len(grid[0])
    
    search_space = []

    # Horizontal (left to right)
    for r in range(num_rows):
        search_space.append("".join(grid[r]))

    # Vertical (top to bottom)
    for c in range(num_cols):
        search_space.append("".join([grid[r][c] for r in range(num_rows)]))

    # Add reversed strings for right-to-left and bottom-to-top searches
    search_space.extend([s[::-1] for s in search_space])

    # Step 4: Find the Metals
    found_metals = set()
    for metal in metals_to_find:
        for s in search_space:
            if metal in s:
                found_metals.add(metal)

    # Step 5: Sort and Format
    sorted_metals = sorted(list(found_metals))
    
    # Print the first 12 metals
    print(", ".join(sorted_metals[:12]))

solve_word_search()