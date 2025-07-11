import math

def solve_knot_tb():
    """
    Calculates the maximal Thurston-Bennequin number for a knot defined by a grid diagram.
    """
    n = 5
    o_coords = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    x_coords = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]

    # P_O maps column i to row of O in that column
    # P_X maps column i to row of X in that column
    P_O = {c[0]: c[1] for c in o_coords}
    P_X = {c[0]: c[1] for c in x_coords}

    # Inverse maps row j to column of O/X in that row
    P_O_inv = {r[1]: r[0] for r in o_coords}
    P_X_inv = {r[1]: r[0] for r in x_coords}

    # --- 1. Calculate the Writhe (w) ---
    writhe = 0
    
    # sv_i is sign of vertical orientation in column i (O -> X)
    sv = {i: 1 if P_X[i] > P_O[i] else -1 for i in range(1, n + 1)}
    
    # sh_j is sign of horizontal orientation in row j (X -> O)
    sh = {j: 1 if P_O_inv[j] > P_X_inv[j] else -1 for j in range(1, n + 1)}

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            # Vertical segment in column i covers y-range (min(O_y, X_y), max(O_y, X_y))
            y_o_i, y_x_i = P_O[i], P_X[i]
            if not (min(y_o_i, y_x_i) < j < max(y_o_i, y_x_i)):
                continue

            # Horizontal segment in row j covers x-range (min(O_x, X_x), max(O_x, X_x))
            x_o_j, x_x_j = P_O_inv[j], P_X_inv[j]
            if not (min(x_o_j, x_x_j) < i < max(x_o_j, x_x_j)):
                continue
            
            # If we are here, a crossing exists at (i, j)
            writhe += sv[i] * sh[j]

    # --- 2. Count North-East (NE) corners (c_NE) ---
    c_ne = 0
    
    # Check corners at O markings (In from West, Out to North)
    for i in range(1, n + 1):
        row_o = P_O[i]
        col_x_in = P_X_inv[row_o]
        row_x_out = P_X[i]
        if i > col_x_in and row_x_out > row_o: # Arriving West->East (i>col_x_in), Departing South->North
            c_ne += 1

    # Check corners at X markings (In from South, Out to East)
    for i in range(1, n + 1):
        row_x = P_X[i]
        row_o_in = P_O[i]
        col_o_out = P_O_inv[row_x]
        if row_x > row_o_in and col_o_out > i: # Arriving South->North, Departing West->East
            c_ne += 1

    # --- 3. Calculate Thurston-Bennequin number (tb) ---
    # Formula: tb = w - c_NE / 2
    tb = writhe - c_ne / 2

    print(f"The writhe (w) is: {writhe}")
    print(f"The number of NE corners (c_NE) is: {c_ne}")
    print(f"The maximal Thurston-Bennequin number (tb) is calculated by the equation w - c_NE / 2:")
    print(f"{writhe} - {c_ne} / 2 = {int(tb)}")


solve_knot_tb()