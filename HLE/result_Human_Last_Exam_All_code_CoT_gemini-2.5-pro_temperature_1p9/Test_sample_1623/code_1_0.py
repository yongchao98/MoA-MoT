import collections

def solve():
    """
    Calculates the Thurston-Bennequin number for the given 5x5 grid diagram.
    """
    n = 5
    # O's at (1,1), (2,2), (3,3), (4,4), (5,5)
    # (col, row)
    o_coords = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]

    # X's at (1,4), (2,5), (3,1), (4,2), (5,3)
    # (col, row)
    x_coords = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]

    # Map row index j to the column of O and X in that row
    o_by_row = {r: c for c, r in o_coords}
    x_by_row = {r: c for c, r in x_coords}

    # Map col index i to the row of O and X in that column
    o_by_col = {c: r for c, r in o_coords}
    x_by_col = {c: r for c, r in x_coords}

    # Step 3: Determine Segment Orientations
    # Horizontal orientations (by row j=1..5)
    horiz_orientations = []
    print("Horizontal Orientations (by row):")
    for j in range(1, n + 1):
        o_col = o_by_row[j]
        x_col = x_by_row[j]
        # West-to-East is positive
        orientation = 1 if x_col > o_col else -1
        horiz_orientations.append(orientation)
        orient_str = "W->E (+)" if orientation == 1 else "E->W (-)"
        print(f"Row {j}: O at column {o_col}, X at column {x_col}. Orientation: {orient_str}")

    num_horiz_pos = horiz_orientations.count(1)
    num_horiz_neg = horiz_orientations.count(-1)
    print(f"\nNumber of W->E (+) rows: {num_horiz_pos}")
    print(f"Number of E->W (-) rows: {num_horiz_neg}\n")

    # Vertical orientations (by column i=1..5)
    vert_orientations = []
    print("Vertical Orientations (by column):")
    for i in range(1, n + 1):
        o_row = o_by_col[i]
        x_row = x_by_col[i]
        # South-to-North is positive
        orientation = 1 if x_row > o_row else -1
        vert_orientations.append(orientation)
        orient_str = "S->N (+)" if orientation == 1 else "N->S (-)"
        print(f"Column {i}: O at row {o_row}, X at row {x_row}. Orientation: {orient_str}")

    num_vert_pos = vert_orientations.count(1)
    num_vert_neg = vert_orientations.count(-1)
    print(f"\nNumber of S->N (+) columns: {num_vert_pos}")
    print(f"Number of N->S (-) columns: {num_vert_neg}\n")

    # Step 4: Count Corners
    # A square (i,j) is a NE corner if row j is W->E (+) and col i is S->N (+)
    num_ne = num_horiz_pos * num_vert_pos
    # A square (i,j) is a NW corner if row j is E->W (-) and col i is S->N (+)
    num_nw = num_horiz_neg * num_vert_pos

    print("Step 5: Calculate the Thurston-Bennequin Number")
    print(f"The number of North-East (NE) corners is the product of the number of W->E rows and S->N columns.")
    print(f"#(NE) = {num_horiz_pos} * {num_vert_pos} = {num_ne}")
    
    print(f"\nThe number of North-West (NW) corners is the product of the number of E->W rows and S->N columns.")
    print(f"#(NW) = {num_horiz_neg} * {num_vert_pos} = {num_nw}")

    # Final calculation
    tb = num_ne - num_nw
    print("\nThe Thurston-Bennequin number is calculated as tb = #(NE) - #(NW).")
    print(f"tb = {num_ne} - {num_nw} = {tb}")

    # Output final answer in the requested format
    print(f"\n<<< {tb} >>>")

solve()