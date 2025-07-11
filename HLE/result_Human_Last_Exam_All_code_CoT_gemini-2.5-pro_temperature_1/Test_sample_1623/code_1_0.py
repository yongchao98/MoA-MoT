def solve_knot_problem():
    """
    This function identifies the knot from the grid diagram and calculates its 
    maximal Thurston-Bennequin number.
    """
    # The problem defines a 5x5 grid.
    # The knot type is determined by the permutation pi mapping the row of an 'O'
    # to the column of the 'X' in the same row.
    # O positions: (1,1), (2,2), (3,3), (4,4), (5,5)
    # X positions: (1,4), (2,5), (3,1), (4,2), (5,3)
    
    # We find the X column for each row.
    # The X positions are given as (column, row).
    x_coords = [(1,4), (2,5), (3,1), (4,2), (5,3)]
    
    # Create a map from row to X-column
    pi_map = {row: col for col, row in x_coords}
    
    # The permutation pi is defined by pi(row) = X_column
    # pi(1) = 3 (from X at (3,1))
    # pi(2) = 4 (from X at (4,2))
    # pi(3) = 5 (from X at (5,3))
    # pi(4) = 1 (from X at (1,4))
    # pi(5) = 2 (from X at (2,5))
    
    print("Step 1: Identify the knot.")
    print("The permutation is pi = (1->3, 2->4, 3->5, 4->1, 5->2).")
    print("This corresponds to the torus knot T(p,q) with p=5, q=2, also known as the 5_1 knot.")
    print("-" * 20)
    
    # Parameters for the T(p,q) knot
    p = 5
    q = 2
    
    print("Step 2: Calculate the maximal Thurston-Bennequin number.")
    print("The formula for the maximal tb of a positive torus knot T(p,q) is: pq - p - q.")
    
    # Calculate the maximal Thurston-Bennequin number
    tb_max = p * q - p - q
    
    print("-" * 20)
    print("Final Calculation:")
    print(f"{p} * {q} - {p} - {q} = {tb_max}")

solve_knot_problem()