import sys

def solve():
    """
    Calculates the maximal Thurston-Bennequin number for a knot from a given grid diagram.
    """
    # Step 1: Define the grid parameters from the problem description.
    # Grid size
    n = 5
    # 1-based coordinates of 'o's and 'x's
    o_coords_1based = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    x_coords_1based = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]

    # We use 0-indexed lists for easier processing in Python.
    # pi maps a column index to a row index.
    # sigma maps a row index to a column index.
    pi_O = [0] * n
    pi_X = [0] * n
    sigma_O = [0] * n
    sigma_X = [0] * n

    # Populate the permutation maps from the 1-based coordinates.
    for c, r in o_coords_1based:
        # c-1 and r-1 convert to 0-based indices
        pi_O[c - 1] = r - 1
        sigma_O[r - 1] = c - 1

    for c, r in x_coords_1based:
        pi_X[c - 1] = r - 1
        sigma_X[r - 1] = c - 1

    # Step 2: Calculate the writhe of the knot diagram.
    # Writhe is the number of crossings. A crossing exists at (i, j) if the
    # vertical strand at column i intersects the horizontal strand at row j.
    writhe = 0
    for i in range(n):  # i is the column index
        for j in range(n):  # j is the row index
            # The vertical line in column 'i' spans between the rows of 'o' and 'x'.
            v_min = min(pi_O[i], pi_X[i])
            v_max = max(pi_O[i], pi_X[i])
            
            # The horizontal line in row 'j' spans between the columns of 'o' and 'x'.
            h_min = min(sigma_O[j], sigma_X[j])
            h_max = max(sigma_O[j], sigma_X[j])
            
            # Check if the grid point (i, j) is an intersection.
            # This happens if j is within the vertical span and i is within the horizontal span.
            is_vertical_span = (v_min < j < v_max)
            is_horizontal_span = (h_min < i < h_max)

            if is_vertical_span and is_horizontal_span:
                writhe += 1

    # Step 3: Calculate the maximal Thurston-Bennequin number.
    # The formula is tb = writhe - n.
    tb = writhe - n
    
    # Step 4: Print the final calculation and result.
    print(f"The grid size is n = {n}.")
    print(f"The calculated writhe of the knot is {writhe}.")
    print("The formula for the maximal Thurston-Bennequin number is: writhe - n.")
    print(f"Maximal Thurston-Bennequin number = {writhe} - {n} = {tb}")

solve()
print("<<<-2>>>")