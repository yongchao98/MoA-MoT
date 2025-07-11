import math

def solve_thurston_bennequin():
    """
    Calculates the maximal Thurston-Bennequin number for a given 5x5 grid diagram.
    """
    # 1. Represent the Grid Diagram
    n = 5
    # o_pos[column] = row
    o_pos = {1:1, 2:2, 3:3, 4:4, 5:5}
    # x_pos[column] = row
    x_pos = {1:4, 2:5, 3:1, 4:2, 5:3}

    # Inverse mappings: inv[row] = column
    o_inv = {v: k for k, v in o_pos.items()}
    x_inv = {v: k for k, v in x_pos.items()}

    # 2. Calculate the Writhe (w)
    writhe = 0
    crossing_details = []

    for i in range(1, n + 1):  # Iterate through columns
        for j in range(1, n + 1):  # Iterate through rows
            # Get the row span for the vertical line in column i
            o_row = o_pos[i]
            x_row = x_pos[i]

            # Get the column span for the horizontal line in row j
            o_col = o_inv[j]
            x_col = x_inv[j]

            # Check if a crossing exists at (i,j)
            # A vertical line passes through (i,j) if j is between the o and x in that column
            is_vert_pass = (j > min(o_row, x_row)) and (j < max(o_row, x_row))
            # A horizontal line passes through (i,j) if i is between the o and x in that row
            is_horiz_pass = (i > min(o_col, x_col)) and (i < max(o_col, x_col))

            if is_vert_pass and is_horiz_pass:
                # Determine vertical strand direction (o to x)
                # +1 for South-to-North (up), -1 for North-to-South (down)
                vert_dir = 1 if x_row > o_row else -1

                # Determine horizontal strand direction (o to x)
                # +1 for West-to-East (right), -1 for East-to-West (left)
                horiz_dir = 1 if x_col > o_col else -1

                # Calculate the sign of the crossing.
                # Convention: vertical strand crosses over horizontal.
                # Positive crossing: up-vert over left-moving horiz, or down-vert over right-moving horiz.
                # This corresponds to sign = - sgn(vert_dir) * sgn(horiz_dir)
                sign = -1 * vert_dir * horiz_dir
                writhe += sign
                crossing_details.append(f"  - Crossing at ({i},{j}) has sign {sign}")

    print(f"The grid size is n = {n}.")
    print("\nThe crossings found are:")
    for detail in crossing_details:
        print(detail)
    
    print(f"\nThe total writhe is the sum of signs: w = {writhe}")
    
    # 3. Calculate Thurston-Bennequin Numbers
    tb1 = writhe - n
    tb2 = -writhe - n

    print("\nThe Thurston-Bennequin number (tb) is calculated as w - n.")
    print(f"tb = {writhe} - {n} = {tb1}")

    print("\nFor the knot with reverse orientation, the calculation is -w - n.")
    print(f"tb_rev = -{writhe} - {n} = {tb2}")
    
    # 4. Find the Maximum
    max_tb = max(tb1, tb2)

    print("\nThe maximal Thurston-Bennequin number is the maximum of these two values.")
    print(f"Maximal tb = max({writhe} - {n}, -{writhe} - {n}) = max({tb1}, {tb2}) = {max_tb}")

    print(f"\n<<< {max_tb} >>>")

solve_thurston_bennequin()