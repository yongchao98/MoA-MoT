def solve_word_search():
    """
    Finds common metals in a grid, sorts them alphabetically, and prints the first 12.
    """
    # 1. Represent the grid from the provided text.
    grid_text = """
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
    L R N M C A D M I U M R
    O Y R U C R E M I M I L
    """
    lines = grid_text.strip().split('\n')
    grid = [line.strip().split(' ') for line in lines]

    # 2. Define the list of common metals to search for.
    metals_to_find = [
        "ALUMINIUM", "BISMUTH", "CADMIUM", "CHROMIUM", "COBALT", "COPPER",
        "GOLD", "IRON", "LEAD", "MAGNESIUM", "MANGANESE", "MERCURY",
        "NICKEL", "PLATINUM", "SILVER", "SODIUM", "TIN", "TITANIUM",
        "TUNGSTEN", "URANIUM", "ZINC"
    ]

    # 3. Generate all possible strings from the grid (rows, columns, and their reverses).
    search_space = []
    
    # Horizontal strings (left-to-right and right-to-left)
    for row in grid:
        line = "".join(row)
        search_space.append(line)
        search_space.append(line[::-1])

    # Vertical strings (top-to-bottom and bottom-to-top)
    if grid:
        num_rows = len(grid)
        num_cols = len(grid[0])
        for j in range(num_cols):
            col_str = "".join([grid[i][j] for i in range(num_rows)])
            search_space.append(col_str)
            search_space.append(col_str[::-1])

    # 4. Find all metals present in the search space.
    found_metals = set()
    for metal in metals_to_find:
        for s in search_space:
            if metal in s:
                found_metals.add(metal)

    # 5. Sort the found metals, take the first 12, and format for printing.
    sorted_metals = sorted(list(found_metals))
    result = ", ".join(sorted_metals[:12])
    print(result)

solve_word_search()