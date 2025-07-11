import math

def solve():
    """
    Calculates the maximal Thurston-Bennequin number for a given 5x5 grid diagram.
    """
    n = 5
    # O's at positions (1,1), (2,2), (3,3), (4,4), (5,5)
    # X's at positions (1,4), (2,5), (3,1), (4,2), (5,3)
    o_positions = [(i, i) for i in range(1, n + 1)]
    x_positions = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]

    # --- Step 1: Create helper maps for coordinates ---
    # These maps find the row/column of an O or X given the other coordinate.
    o_map_c = {c: r for c, r in o_positions}
    x_map_c = {c: r for c, r in x_positions}
    o_map_r = {r: c for c, r in o_positions}
    x_map_r = {r: c for c, r in x_positions}

    # --- Step 2: Trace the knot path to determine segment orientations ---
    # The knot path alternates between vertical and horizontal segments.
    # Vertical segments connect an X and an O in the same column.
    # Horizontal segments connect an O and an X in the same row.
    vert_dirs = {}
    horiz_dirs = {}

    start_pos = x_positions[0]
    curr_pos = start_pos
    
    # Follow the path for 2n segments to find all orientations
    for _ in range(n):
        # Move from X to O (vertically)
        x_col, x_row = curr_pos
        o_row = o_map_c[x_col]
        vert_dirs[x_col] = 1 if o_row > x_row else -1  # +1 for Up, -1 for Down
        o_pos_curr = (x_col, o_row)
        
        # Move from O to X (horizontally)
        o_col, o_row_curr = o_pos_curr
        x_col_next = x_map_r[o_row_curr]
        horiz_dirs[o_row_curr] = 1 if x_col_next > o_col else -1  # +1 for Right, -1 for Left
        curr_pos = (x_col_next, o_row_curr)


    # --- Step 3: Calculate the writhe (w) ---
    # Writhe is the sum of signs of crossings.
    # A crossing at (i,j) exists if vertical segment V_i crosses horizontal segment H_j.
    # Sign = sign(V_i direction) * sign(H_j direction)
    writhe = 0
    crossings = []
    for i in range(1, n + 1):  # column for vertical segment
        for j in range(1, n + 1):  # row for horizontal segment
            # Get vertical segment y-range
            y1_v = o_map_c[i]
            y2_v = x_map_c[i]
            y_min_v, y_max_v = min(y1_v, y2_v), max(y1_v, y2_v)

            # Get horizontal segment x-range
            x1_h = o_map_r[j]
            x2_h = x_map_r[j]
            x_min_h, x_max_h = min(x1_h, x2_h), max(x1_h, x2_h)

            # Check if the interiors of the segments cross
            if x_min_h < i < x_max_h and y_min_v < j < y_max_v:
                sign = vert_dirs[i] * horiz_dirs[j]
                writhe += sign
                crossings.append(((i,j), sign))
    
    # --- Step 4: Calculate the rotation number (ro) ---
    # ro is the sum of contributions from cusps (at O and X markers).
    # Up-cusp contributes +1, Down-cusp contributes -1.
    rotation_number = 0
    
    # At 'O' markers, path goes from Vertical to Horizontal
    # A turn from Up(S->N) is an up-cusp. A turn from Down(N->S) is a down-cusp.
    for col in range(1, n+1):
        rotation_number += vert_dirs[col]
        
    # At 'X' markers, path goes from Horizontal to Vertical
    # A turn to Up(S->N) is an up-cusp. A turn to Down(N->S) is a down-cusp.
    for col in range(1, n+1):
        rotation_number += vert_dirs[col]

    # --- Step 5: Calculate the Thurston-Bennequin number (tb) ---
    # The formula is tb = w - ro/2
    tb = writhe - rotation_number / 2
    
    print("To find the maximal Thurston-Bennequin number, we perform the following calculations:")
    print(f"1. The writhe (w) is the sum of signed crossings. The crossings are at {[c[0] for c in crossings]}.")
    print(f"   The calculated writhe is: w = {writhe}")
    print(f"2. The rotation number (ro) is calculated from the cusps in the diagram.")
    print(f"   The calculated rotation number is: ro = {int(rotation_number)}")
    print(f"3. The Thurston-Bennequin number is calculated using the formula: tb = w - ro / 2")
    print(f"   tb = {writhe} - {int(rotation_number)} / 2 = {int(tb)}")
    
    return int(tb)

final_answer = solve()
print(f"\nThe final computed maximal Thurston-Bennequin number is {final_answer}.")
print(f'<<<{final_answer}>>>')