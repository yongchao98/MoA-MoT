import numpy as np

def solve_knot_problem():
    """
    Solves the knot theory problem by calculating winding numbers and the specified sum.
    """
    # Step 1: Define the grid diagram for the left-hand trefoil knot.
    # O's on the anti-diagonal. This grid is the mirror image of the standard
    # right-hand trefoil grid.
    # Grid size
    n = 3
    # Marker positions (row, col), with (0,0) at top-left
    # 'X' = 1, 'O' = -1
    grid = np.zeros((n, n), dtype=int)
    grid[0, 1] = 1  # X
    grid[1, 0] = 1  # X
    grid[2, 2] = 1  # X
    grid[0, 2] = -1 # O
    grid[1, 1] = -1 # O
    grid[2, 0] = -1 # O

    # Step 2: Calculate the winding number matrix w.
    # The formulas calculate w for the positive knot (mirror image).
    # We will flip the signs at the end for the left-hand trefoil.
    w_pos = np.zeros((n + 1, n + 1), dtype=int)
    
    # Use the row-update rule: w(i+1, j) = w(i, j) + sum_{k<j} (X_{i,k} - O_{i,k})
    for i in range(n):
        for j in range(n + 1):
            # Contribution from markers in row i to the left of column j
            row_sum = 0
            if j > 0:
                for k in range(j):
                    if grid[i, k] == 1: # X
                        row_sum += 1
                    elif grid[i, k] == -1: # O
                        row_sum -= 1
            w_pos[i + 1, j] = w_pos[i, j] + row_sum

    # For the left-hand trefoil (a negative knot), we flip the signs.
    w = -w_pos

    # Step 3: Identify the corner groups by calculating the k-matrix.
    # k[i,j] = number of markers adjacent to lattice point (i,j).
    k_matrix = np.zeros((n + 1, n + 1), dtype=int)
    for i in range(n + 1):
        for j in range(n + 1):
            count = 0
            # Check the four adjacent cells: (i-1,j-1), (i-1,j), (i,j-1), (i,j)
            # Top-left cell
            if i > 0 and j > 0 and grid[i - 1, j - 1] != 0:
                count += 1
            # Top-right cell
            if i > 0 and j < n and grid[i - 1, j] != 0:
                count += 1
            # Bottom-left cell
            if i < n and j > 0 and grid[i, j - 1] != 0:
                count += 1
            # Bottom-right cell
            if i < n and j < n and grid[i, j] != 0:
                count += 1
            k_matrix[i, j] = count

    # Step 4: Calculate the final sum.
    total_sum = 0
    sum_terms = []
    
    print("The final sum is calculated as follows:")
    print("Sum = 1 * (Sum of w in mho_1) + 2 * (Sum of w in mho_2) + 3 * (Sum of w in mho_3) + 4 * (Sum of w in mho_4)")
    print("")

    for k in range(1, 5):
        w_in_mho_k = w[k_matrix == k]
        sum_w_k = np.sum(w_in_mho_k)
        term = k * sum_w_k
        total_sum += term
        
        # Format the numbers for the equation string
        w_list_str = [str(val) for val in w_in_mho_k]
        if not w_list_str:
            w_list_str.append("0")
            
        sum_str = " + ".join(f"({val})" if val < 0 else str(val) for val in w_in_mho_k)
        if not sum_str:
            sum_str = "0"

        print(f"For k = {k}:")
        print(f"  The set of winding numbers is mho_{k} = {{{', '.join(w_list_str)}}}")
        print(f"  The term is {k} * ({sum_str}) = {k} * ({sum_w_k}) = {term}")
        print("")

    print(f"The total value of the expression is: {total_sum}")

solve_knot_problem()
<<<-8>>>