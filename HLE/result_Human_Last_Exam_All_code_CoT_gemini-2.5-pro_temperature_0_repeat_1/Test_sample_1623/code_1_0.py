import numpy as np

def solve_knot_tb():
    """
    Calculates the maximal Thurston-Bennequin number for a knot
    defined by a 5x5 grid diagram.
    """
    n = 5
    # 1-based coordinates from the problem description
    o_coords_1based = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    x_coords_1based = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]

    # Convert to 0-based coordinates (col, row)
    o_pos = [(c - 1, r - 1) for c, r in o_coords_1based]
    x_pos = [(c - 1, r - 1) for c, r in x_coords_1based]

    # Create lookup dictionaries for easy access
    o_row_from_col = {c: r for c, r in o_pos}
    o_col_from_row = {r: c for c, r in o_pos}
    x_row_from_col = {c: r for c, r in x_pos}
    x_col_from_row = {r: c for c, r in x_pos}
    
    all_marker_pos = set(o_pos + x_pos)

    # --- 2. Calculate the Writhe (w) ---
    writhe = 0
    v_dir = {i: np.sign(o_row_from_col[i] - x_row_from_col[i]) for i in range(n)}
    h_dir = {j: np.sign(x_col_from_row[j] - o_col_from_row[j]) for j in range(n)}

    crossings = []
    for i in range(n):  # column of potential crossing
        for j in range(n):  # row of potential crossing
            # Check for vertical segment crossing row j
            y_o, y_x = o_row_from_col[i], x_row_from_col[i]
            if not (min(y_o, y_x) < j < max(y_o, y_x)):
                continue

            # Check for horizontal segment crossing column i
            x_o, x_x = o_col_from_row[j], x_col_from_row[j]
            if not (min(x_o, x_x) < i < max(x_o, x_x)):
                continue
            
            sign = v_dir[i] * h_dir[j]
            writhe += sign
            crossings.append(((i, j), sign))

    # --- 3. Calculate North-East Corners (c_NE) ---
    c_ne = 0
    ne_corners = []
    for i in range(n):
        for j in range(n):
            if (i, j) in all_marker_pos:
                continue
            # Check for X to the North (same column, higher row)
            has_x_north = x_row_from_col[i] > j
            # Check for O to the East (same row, greater column)
            has_o_east = o_col_from_row[j] > i
            
            if has_x_north and has_o_east:
                c_ne += 1
                ne_corners.append((i, j))

    # --- 4. Calculate South-West Corners (c_SW) ---
    c_sw = 0
    sw_corners = []
    for i in range(n):
        for j in range(n):
            if (i, j) in all_marker_pos:
                continue
            # Check for X to the South (same column, lower row)
            has_x_south = x_row_from_col[i] < j
            # Check for O to the West (same row, lesser column)
            has_o_west = o_col_from_row[j] < i

            if has_x_south and has_o_west:
                c_sw += 1
                sw_corners.append((i, j))

    # --- 5. Compute the Maximal Thurston-Bennequin Number ---
    tb = writhe - c_ne - c_sw

    # --- 6. Final Output ---
    print("Step 1: Calculating the writhe (w)")
    print(f"Found {len(crossings)} crossings:")
    for (pos, sign) in crossings:
        print(f"  - Crossing at (col={pos[0]+1}, row={pos[1]+1}) with sign {sign}")
    print(f"Total writhe w = {writhe}\n")
    
    print("Step 2: Calculating North-East corners (c_NE)")
    print(f"Found {c_ne} NE corners at (col, row) positions (1-based):")
    print(f"  {[(c+1, r+1) for c, r in ne_corners]}")
    print(f"Total NE corners c_NE = {c_ne}\n")

    print("Step 3: Calculating South-West corners (c_SW)")
    print(f"Found {c_sw} SW corners at (col, row) positions (1-based):")
    print(f"  {[(c+1, r+1) for c, r in sw_corners]}")
    print(f"Total SW corners c_SW = {c_sw}\n")

    print("Step 4: Calculating the maximal Thurston-Bennequin number (TB)")
    print("Formula: TB = w - c_NE - c_SW")
    print(f"TB = {writhe} - {c_ne} - {c_sw} = {tb}")

solve_knot_tb()
<<< -4 >>>