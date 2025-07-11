import re

def find_metals_in_grid():
    """
    Finds all common metals in a grid, sorts them alphabetically, and prints the first 12.
    Words can be horizontal or vertical, in both forward and reverse directions.
    """
    # Step 1: Define and parse the grid from the problem description.
    # The grid is represented as a multiline string.
    grid_str = """
    N  T  I  T  A  N  I  U  M  M  I  T
    E  C  D  M  C  R  A  S  G  A  T  Z
    T  O  M  T  I  I  E  T  O  N  A  I
    S  B  I  C  T  L  S  I  L  G  I  N
    G  A  M  A  N  E  Z  L  D  A  R  C
    N  L  N  I  S  A  C  I  U  N  I  G
    U  T  M  I  I  D  I  I  A  E  D  M
    T  L  R  E  V  L  I  S  C  S  I  C
    S  I  C  T  A  I  R  O  N  E  U  O
    M  P  M  U  I  S  E  N  G  A  M  P
    C  M  U  N  I  T  A  L  P  N  P  P
    R  C  C  G  M  N  I  C  K  E  L  E
    L  R  N  M  C  A  D  M  I  U  M  M
    O  Y  R  U  C  R  E  M  I  M  I  L
    """
    # Process the string into a 2D list of characters.
    grid = [re.split(r'\s+', line.strip()) for line in grid_str.strip().split('\n')]

    # Step 2: Create a list of all possible strings to search in.
    # This includes rows, columns, and their reversed versions.
    search_space = []
    
    # Horizontal strings (rows)
    rows = ["".join(row) for row in grid]
    for row in rows:
        search_space.append(row)
        search_space.append(row[::-1])

    # Vertical strings (columns)
    num_rows = len(grid)
    num_cols = len(grid[0])
    cols = []
    for j in range(num_cols):
        col_str = "".join([grid[i][j] for i in range(num_rows)])
        cols.append(col_str)
    
    for col in cols:
        search_space.append(col)
        search_space.append(col[::-1])

    # Step 3: Define a comprehensive list of common metals.
    metals = [
        "ALUMINIUM", "ANTIMONY", "BARIUM", "BERYLLIUM", "BISMUTH", "CADMIUM",
        "CALCIUM", "CHROMIUM", "COBALT", "COPPER", "GOLD", "IRIDIUM", "IRON", 
        "LEAD", "LITHIUM", "MAGNESIUM", "MANGANESE", "MERCURY", "MOLYBDENUM", 
        "NICKEL", "PALLADIUM", "PLATINUM", "POTASSIUM", "SILVER", "SODIUM", "STEEL",
        "TIN", "TITANIUM", "TUNGSTEN", "URANIUM", "VANADIUM", "ZINC", "BRASS", "BRONZE"
    ]

    # Step 4: Find all metals from the list that are present in the grid.
    found_metals = set()
    for metal in metals:
        for s in search_space:
            if metal in s:
                found_metals.add(metal)
                break
    
    # Step 5: Sort the found metals alphabetically.
    sorted_metals = sorted(list(found_metals))

    # Step 6: Print the first 12 metals from the sorted list.
    print(", ".join(sorted_metals[:12]))

if __name__ == "__main__":
    find_metals_in_grid()