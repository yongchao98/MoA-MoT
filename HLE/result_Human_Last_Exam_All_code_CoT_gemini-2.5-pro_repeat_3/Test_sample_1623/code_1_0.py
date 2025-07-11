def solve_knot_tb():
    """
    Calculates the maximal Thurston-Bennequin number for a knot
    specified by a 5x5 grid diagram.
    """
    # Step 1: Define the grid size and the positions of the markers.
    # The grid is 5x5.
    n = 5
    # The O's are on the diagonal (1,1), (2,2), ..., (5,5), which makes this a normalized grid.
    # The X's positions are given as (column, row) with 1-based indexing.
    x_positions = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]

    print(f"The grid size is n = {n}.")
    print(f"The X markers are at positions: {x_positions}")
    print("-" * 30)

    # Step 2: Construct the permutation for the X markers.
    # The permutation pi_x maps a column index `i` to the row index of the X in that column.
    # We use a 0-indexed list, so we adjust the 1-indexed positions from the problem.
    pi_x = [0] * n
    for col, row in x_positions:
        pi_x[col - 1] = row
    
    print(f"The corresponding permutation for X markers is: {pi_x}")

    # Step 3: Calculate the writhe (w) of the diagram.
    # For a normalized grid diagram, the writhe is the number of inversions
    # in the permutation pi_x. An inversion is a pair of indices (i, j)
    # such that i < j and pi_x[i] > pi_x[j].
    writhe = 0
    for i in range(n):
        for j in range(i + 1, n):
            if pi_x[i] > pi_x[j]:
                writhe += 1

    print(f"The writhe (w) is the number of inversions in the permutation.")
    print(f"Calculated writhe: w = {writhe}")
    
    # Step 4: Calculate the Thurston-Bennequin number (tb).
    # The formula is tb = w - n.
    # As explained in the plan, this gives the maximal tb number for this knot.
    tb = writhe - n
    
    print("\nThe maximal Thurston-Bennequin number is calculated as tb = w - n.")
    print("The final equation is:")
    
    # Step 5: Print the final calculation and the result, showing each number.
    print(f"{writhe} - {n} = {tb}")

solve_knot_tb()
<<<1>>>