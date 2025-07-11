import numpy as np

def solve_knot_tb():
    """
    Calculates the maximal Thurston-Bennequin number for a given 5x5 grid diagram.
    """
    n = 5
    # O positions are given as (column, row)
    # The problem statement says (i,j) means i-th column and j-th row.
    o_positions = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    x_positions = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]

    # Create mappings for easy lookup.
    # sigma_o[col] = row
    sigma_o = {pos[0]: pos[1] for pos in o_positions}
    sigma_x = {pos[0]: pos[1] for pos in x_positions}

    # o_x[row] = col
    o_x = {pos[1]: pos[0] for pos in o_positions}
    x_x = {pos[1]: pos[0] for pos in x_positions}
    
    writhe = 0
    crossings = []

    # Iterate through all potential crossing locations (i, j)
    # where i is the column of the vertical strand
    # and j is the row of the horizontal strand.
    for i in range(1, n + 1):
        # Vertical segment in column i connects (i, sigma_o[i]) and (i, sigma_x[i])
        yo = sigma_o[i]
        yx = sigma_x[i]
        
        for j in range(1, n + 1):
            # Horizontal segment in row j connects (o_x[j], j) and (x_x[j], j)
            xo = o_x[j]
            xx = x_x[j]

            # A crossing exists at (i, j) if the vertical segment spans row j
            # and the horizontal segment spans column i.
            # The intervals are open, meaning the crossing can't be at a marker.
            vertical_spans_j = min(yo, yx) < j < max(yo, yx)
            horizontal_spans_i = min(xo, xx) < i < max(xo, xx)

            if vertical_spans_j and horizontal_spans_i:
                # Determine orientation signs. +1 for increasing coord, -1 for decreasing.
                # The knot is oriented from O to X.
                vertical_dir_sign = np.sign(yx - yo)
                horizontal_dir_sign = np.sign(xx - xo)
                
                # Sign of crossing with Vertical over Horizontal is -det(v_vert, v_horiz)
                # where v_vert=(0, vert_dir), v_horiz=(horiz_dir, 0)
                # sign = - (0*0 - vert_dir * horiz_dir) = vert_dir * horiz_dir. This is a common convention.
                # Let's check the other convention sign = - sign(V_y) * sign(H_x)
                crossing_sign = -1 * vertical_dir_sign * horizontal_dir_sign
                writhe += crossing_sign
                crossings.append(f"Crossing at ({i},{j}) with sign {int(crossing_sign)}")

    print(f"Grid size n = {n}")
    print("Found the following crossings:")
    for c in crossings:
        print(c)
    
    tb = writhe - n
    
    print("\nThe writhe (w) is the sum of the signs of all crossings.")
    # Build the writhe equation string from the signs
    signs = [str(int(c.split(' ')[-1])) for c in crossings]
    writhe_eq = " + ".join(signs).replace('+ -', '- ')
    print(f"w = {writhe_eq} = {int(writhe)}")

    print("\nThe Thurston-Bennequin number (tb) is calculated as w - n.")
    print(f"tb = {int(writhe)} - {n} = {int(tb)}")

solve_knot_tb()