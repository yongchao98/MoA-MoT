import numpy as np

def solve():
    """
    Calculates the maximal Thurston-Bennequin number for a given 5x5 grid diagram.
    """
    n = 5
    
    # 'o' positions: (1,1), (2,2), (3,3), (4,4), (5,5)
    # 'x' positions: (1,4), (2,5), (3,1), (4,2), (5,3)
    # Using 1-based indexing for clarity, will convert to 0-based for list access.
    
    # Map column index to row index for 'o's and 'x's
    o_col_to_row = {i: i for i in range(1, n + 1)}
    x_col_to_row = {1: 4, 2: 5, 3: 1, 4: 2, 5: 3}
    
    # Map row index to column index for 'o's and 'x's
    o_row_to_col = {i: i for i in range(1, n + 1)}
    x_row_to_col = {4: 1, 5: 2, 1: 3, 2: 4, 3: 5}

    # --- 1. Calculate the Writhe (w) ---
    writhe = 0
    # Iterate through each potential crossing point (i, j)
    for i in range(1, n + 1):  # i is the column index
        for j in range(1, n + 1):  # j is the row index
            
            # Vertical segment in column i is between o_row and x_row
            o_row = o_col_to_row[i]
            x_row = x_col_to_row[i]
            
            # Horizontal segment in row j is between o_col and x_col
            o_col = o_row_to_col[j]
            x_col = x_row_to_col[j]
            
            # Check if a crossing exists at (i, j)
            # j must be strictly between the rows of the vertical segment
            # i must be strictly between the columns of the horizontal segment
            is_crossing = (min(o_row, x_row) < j < max(o_row, x_row) and
                           min(o_col, x_col) < i < max(o_col, x_col))
                           
            if is_crossing:
                # Orientation: O -> H -> X -> V -> O
                # Vertical sign: sgn(O_row - X_row)
                # Horizontal sign: sgn(X_col - O_col)
                v_sign = np.sign(o_row - x_row)
                h_sign = np.sign(x_col - o_col)
                writhe += v_sign * h_sign

    # --- 2. Calculate the Rotation Number (rot) ---
    # rot = 0.5 * ( sum(sgn(x'_j - o'_j)) - sum(sgn(x_i - o_i)) )
    # Term 1: Sum over rows
    sum_h_signs = 0
    for j in range(1, n + 1):
        sum_h_signs += np.sign(x_row_to_col[j] - o_row_to_col[j])
    
    # Term 2: Sum over columns. Note the formula is sgn(x_i-o_i), not o_i-x_i
    sum_v_signs = 0
    for i in range(1, n + 1):
        sum_v_signs += np.sign(x_col_to_row[i] - o_col_to_row[i])
        
    rotation_number = 0.5 * (sum_h_signs - sum_v_signs)
    
    # --- 3. Calculate Thurston-Bennequin Number (tb) ---
    # tb = w - rot
    tb_number = writhe - rotation_number

    # Print the final equation
    print(f"Writhe (w) = {int(writhe)}")
    print(f"Rotation Number (rot) = {int(rotation_number)}")
    print("Thurston-Bennequin Number (tb) = w - rot")
    print(f"{int(writhe)} - {int(rotation_number)} = {int(tb_number)}")
    print(f"\nThe maximal Thurston-Bennequin number is: {int(tb_number)}")

solve()