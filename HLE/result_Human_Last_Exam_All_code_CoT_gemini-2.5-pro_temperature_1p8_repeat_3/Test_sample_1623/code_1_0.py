import math

def calculate_tb():
    """
    Calculates the Thurston-Bennequin number for the given 5x5 grid diagram.
    """
    n = 5
    
    # O locations: (col, row)
    o_coords = {(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)}
    # X locations: (col, row)
    x_coords = {(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)}

    # Map column to row for O and X markers
    o_col_to_row = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
    x_col_to_row = {1: 4, 2: 5, 3: 1, 4: 2, 5: 3}

    # Map row to column for O and X markers
    o_row_to_col = {v: k for k, v in o_col_to_row.items()}
    x_row_to_col = {v: k for k, v in x_col_to_row.items()}
    
    writhe = 0
    crossings = []

    # Iterate through all possible grid cells (i,j) for crossings
    for i in range(1, n + 1):  # i = column
        for j in range(1, n + 1):  # j = row
            # A cell must be empty to be a crossing
            if (i, j) in o_coords or (i, j) in x_coords:
                continue

            # Check for vertical span
            yo, yx = o_col_to_row[i], x_col_to_row[i]
            if not (min(yo, yx) < j < max(yo, yx)):
                continue

            # Check for horizontal span
            xo, xx = o_row_to_col[j], x_row_to_col[j]
            if not (min(xo, xx) < i < max(xo, xx)):
                continue
                
            # If we are here, (i, j) is a crossing
            crossings.append((i,j))
            
            # Determine the sign of the crossing
            # Vertical strand O -> X
            vertical_dir = math.copysign(1, yx - yo)
            # Horizontal strand X -> O
            horizontal_dir = math.copysign(1, xo - xx)
            
            # The sign of the crossing is the product of the directions
            # with standard {h, v} basis sign convention.
            sign = vertical_dir * horizontal_dir
            writhe += sign

    tb = writhe - n
    
    print("Thurston-Bennequin Number Calculation")
    print(f"Grid size n: {n}")
    print(f"Calculated Writhe w(D): {int(writhe)}")
    print(f"Formula: tb = w(D) - n")
    print(f"The equation is: {int(writhe)} - {n} = {int(tb)}")
    print(f"The maximal Thurston-Bennequin number is {int(tb)}.")

calculate_tb()