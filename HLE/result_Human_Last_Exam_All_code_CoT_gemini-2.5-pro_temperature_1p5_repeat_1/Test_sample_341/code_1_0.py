import collections

def solve_word_search():
    """
    Finds common metal names in a grid, sorts them alphabetically,
    and prints the first 12.
    """
    # 1. Represent the grid from the problem description.
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
    L R N M C A D M I U M R
    O Y R U C R E M I M I L
    """
    grid = [line.strip().split() for line in grid_str.strip().split('\n')]

    # 2. Define the list of common metals to search for.
    metals_to_find = [
        "ALUMINUM", "BISMUTH", "BRASS", "BRONZE", "CADMIUM", "CHROMIUM",
        "COBALT", "COPPER", "GOLD", "IRIDIUM", "IRON", "LEAD", "LITHIUM",
        "MAGNESIUM", "MANGANESE", "MERCURY", "NICKEL", "PALLADIUM",
        "PLATINUM", "SILVER", "SODIUM", "STEEL", "TIN", "TITANIUM",
        "TUNGSTEN", "URANIUM", "ZINC"
    ]

    # 3. Generate all possible search strings from the grid.
    search_space = []
    
    # Horizontal strings (left to right)
    rows_lr = [''.join(row) for row in grid]
    search_space.extend(rows_lr)
    
    # Horizontal strings (right to left)
    rows_rl = [s[::-1] for s in rows_lr]
    search_space.extend(rows_rl)

    # Vertical strings (top to bottom)
    cols_tb = [''.join(col) for col in zip(*grid)]
    search_space.extend(cols_tb)

    # Vertical strings (bottom to top)
    cols_bt = [s[::-1] for s in cols_tb]
    search_space.extend(cols_bt)

    # 4. Search for the metals in the generated strings.
    found_metals = set()
    for metal in metals_to_find:
        for line in search_space:
            if metal in line:
                found_metals.add(metal)
                # Once found, no need to check other orientations for the same metal
                break

    # 5. Sort the found metals alphabetically and get the first 12.
    sorted_found_metals = sorted(list(found_metals))
    result = sorted_found_metals[:12]

    # Print the final result in the specified format.
    print(", ".join(result))

solve_word_search()
<<<CADMIUM, COBALT, COPPER, GOLD, IRON, MAGNESIUM, MERCURY, NICKEL, PLATINUM, SILVER, TIN, TITANIUM>>>