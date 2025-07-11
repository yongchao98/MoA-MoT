import itertools

def solve_braid_index():
    """
    Calculates the braid index of a knot from its grid diagram representation.
    """
    # The grid number
    n = 7

    # Define the positions of the 'o' and 'x' markers as sets for efficient lookup.
    # Positions are (column, row) from (1,1) to (7,7).
    o_pos = frozenset([
        (1, 1), (2, 7), (3, 4), (4, 5), (5, 3), (6, 6), (7, 2)
    ])
    x_pos = frozenset([
        (1, 2), (2, 6), (3, 3), (4, 1), (5, 7), (6, 5), (7, 4)
    ])
    all_pos = o_pos.union(x_pos)

    def is_empty(i1, j1, i2, j2):
        """
        Checks if the interior of the rectangle defined by (i1, j1) and (i2, j2)
        is empty of any markers.
        """
        for i in range(i1 + 1, i2):
            for j in range(j1 + 1, j2):
                if (i, j) in all_pos:
                    return False
        return True

    # --- Calculate N_plus ---
    # N_plus is the count of empty rectangles with 'o's on the main diagonal
    # and 'x's on the anti-diagonal.
    n_plus = 0
    # Iterate through all unique pairs of 'o' markers
    for o1, o2 in itertools.combinations(o_pos, 2):
        i1, i2 = sorted([o1[0], o2[0]])
        j1, j2 = sorted([o1[1], o2[1]])

        # Check if the corners form a rectangle with the correct markers
        # Bottom-left and top-right must be 'o's
        if (i1, j1) in o_pos and (i2, j2) in o_pos:
            # Top-left and bottom-right must be 'x's
            if (i1, j2) in x_pos and (i2, j1) in x_pos:
                if is_empty(i1, j1, i2, j2):
                    n_plus += 1

    # --- Calculate N_minus ---
    # N_minus is the count of empty rectangles with 'x's on the main diagonal
    # and 'o's on the anti-diagonal.
    n_minus = 0
    # Iterate through all unique pairs of 'x' markers
    for x1, x2 in itertools.combinations(x_pos, 2):
        i1, i2 = sorted([x1[0], x2[0]])
        j1, j2 = sorted([x1[1], x2[1]])

        # Check if the corners form a rectangle with the correct markers
        # Bottom-left and top-right must be 'x's
        if (i1, j1) in x_pos and (i2, j2) in x_pos:
            # Top-left and bottom-right must be 'o's
            if (i1, j2) in o_pos and (i2, j1) in o_pos:
                if is_empty(i1, j1, i2, j2):
                    n_minus += 1
    
    # Apply Dynnikov's formula for the braid index
    braid_index = n - max(n_plus, n_minus)

    print(f"The grid number is n = {n}")
    print(f"Number of positive empty rectangles found: N_plus = {n_plus}")
    print(f"Number of negative empty rectangles found: N_minus = {n_minus}")
    print("\nThe braid index is calculated using the formula: n - max(N_plus, N_minus)")
    print(f"Braid Index = {n} - max({n_plus}, {n_minus}) = {braid_index}")

solve_braid_index()