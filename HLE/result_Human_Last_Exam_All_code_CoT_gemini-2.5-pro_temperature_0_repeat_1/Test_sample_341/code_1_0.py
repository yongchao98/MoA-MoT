def find_metals_in_grid():
    """
    Finds common metals in a grid, sorts them alphabetically, and prints the first 12.
    """
    # The grid of letters provided in the problem
    grid_str = """
    N T I T A N I U M M I T
    E C D M C R A S G A T Z
    T O M T I I E T O N A I
    S B I C T L S I L G I N
    G A M A N E Z L D A R C
    N L N I S A C I U N I G
    U T M I I D I I A E D M
    T L R E V L I S C S I C
    S I C T A I R O N E U O
    M P M U I S E N G A M P
    C M U N I T A L P N P P
    R C C G M N I C K E L E
    L R N M C A D M I U M M
    O Y R U C R E M I M I L
    """
    grid = [line.strip().split() for line in grid_str.strip().split('\n')]

    # A list of common metals to search for
    metals = [
        "ALUMINUM", "BISMUTH", "CADMIUM", "CHROMIUM", "COBALT", "COPPER",
        "GOLD", "IRON", "LEAD", "LITHIUM", "MAGNESIUM", "MANGANESE",
        "MERCURY", "NICKEL", "PLATINUM", "SILVER", "SODIUM", "TIN",
        "TITANIUM", "TUNGSTEN", "URANIUM", "ZINC"
    ]

    found_metals = set()
    rows = len(grid)
    cols = len(grid[0])

    # Generate all horizontal and vertical strings from the grid
    lines_to_search = []
    # Horizontal lines (left to right)
    for r in range(rows):
        lines_to_search.append("".join(grid[r]))
    
    # Vertical lines (top to bottom)
    for c in range(cols):
        lines_to_search.append("".join(grid[r][c] for r in range(rows)))

    # Search for each metal and its reverse in the generated lines
    for metal in metals:
        reversed_metal = metal[::-1]
        for line in lines_to_search:
            if metal in line or reversed_metal in line:
                found_metals.add(metal)

    # Sort the found metals alphabetically
    sorted_metals = sorted(list(found_metals))
    
    # Select the first 12 metals from the sorted list
    result_to_print = sorted_metals[:12]
    
    # Print the result in the specified format
    print(", ".join(result_to_print))

find_metals_in_grid()
<<<CADMIUM, COBALT, GOLD, IRON, LEAD, MAGNESIUM, MERCURY, NICKEL, PLATINUM, SILVER, TIN, TITANIUM>>>