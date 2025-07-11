import collections

def solve_knot_grid_sum():
    """
    Solves the problem of calculating the weighted sum of winding numbers for a specific trefoil knot grid diagram.
    """
    n = 3  # Grid size

    # Step 1: Define the grid diagram for the left-hand trefoil knot.
    # Cell indices (row, col) are 1-based, from bottom-left.
    # O's are on the diagonal from bottom-left to top-right.
    # The X's are placed to form the left-hand trefoil knot.
    # Value: 1 for 'O', -1 for 'X', 0 for empty.
    markers = {
        (1, 1): 1, (2, 2): 1, (3, 3): 1,  # The 'O' markers
        (2, 1): -1, (3, 2): -1, (1, 3): -1 # The 'X' markers
    }

    print("Step 1: Grid Diagram Configuration")
    grid_repr = [['.' for _ in range(n)] for _ in range(n)]
    for (r,c), val in markers.items():
        if val == 1:
            grid_repr[n-r][c-1] = 'O'
        elif val == -1:
            grid_repr[n-r][c-1] = 'X'
    print("The 3x3 grid has markers at these cell coordinates (row, col) from bottom-left:")
    print("O: (1,1), (2,2), (3,3)")
    print("X: (2,1), (3,2), (1,3)")
    print("Grid visual:")
    for row in grid_repr:
        print("  ".join(row))
    print("-" * 20)

    # Step 2: Calculate winding numbers w(r, c) for each vertex (0 <= r, c <= n).
    # Vertex indices (r, c) are 0-based from bottom-left.
    winding_numbers = collections.defaultdict(int)

    # The equation for cell (r,c) relates the vertices around it:
    # w(r,c) - w(r, c-1) + w(r-1, c-1) - w(r-1, c) = marker_value(r,c)
    # Solve for w(r,c):
    # w(r,c) = w(r, c-1) - w(r-1, c-1) + w(r-1, c) + marker_value(r,c)
    # We can iterate from r=1 to n and c=1 to n.
    for r in range(1, n + 1):
        for c in range(1, n + 1):
            epsilon = markers.get((r, c), 0)
            w_nw = winding_numbers[(r, c - 1)]
            w_sw = winding_numbers[(r - 1, c - 1)]
            w_se = winding_numbers[(r - 1, c)]
            winding_numbers[(r, c)] = w_nw - w_sw + w_se + epsilon

    print("Step 2: Calculated Winding Numbers w(r,c)")
    w_grid = [[winding_numbers[(r,c)] for c in range(n+1)] for r in range(n, -1, -1)]
    print("The winding numbers form a 4x4 grid (vertex indices r,c from 0 to 3):")
    for row in w_grid:
        print(" ".join(map(str, row)))
    print("-" * 20)
    
    # Step 3: Classify vertices by k and sum winding numbers.
    sums_by_k = {1: 0, 2: 0, 3: 0, 4: 0}
    mho = {1: [], 2: [], 3: [], 4: []}

    for r in range(n + 1):
        for c in range(n + 1):
            # Check the four cells sharing vertex (r,c).
            # These are cells (r,c), (r+1,c), (r,c+1), (r+1,c+1) in 0-based cell coordinates.
            # In our 1-based cell system, these are (r+1, c+1), (r+1,c), (r,c+1), (r,c).
            
            k = 0
            # Cells are C(i, j) where 1 <= i,j <= n.
            # Vertex (r,c) is adjacent to cells C(r,c), C(r+1,c), C(r,c+1), C(r+1,c+1)
            # using 1-based indexing for cells:
            adjacent_cells = [(r, c), (r, c + 1), (r + 1, c), (r + 1, c + 1)]
            for cell_r, cell_c in adjacent_cells:
                if (cell_r, cell_c) in markers:
                    k += 1
            
            if k in sums_by_k:
                w = winding_numbers[(r,c)]
                sums_by_k[k] += w
                mho[k].append(w)
    
    print(f"Step 3: Sum of winding numbers for each k")
    for k in sorted(sums_by_k.keys()):
        print(f"k={k}: The set of winding numbers is \u03D8_{k} = {mho[k]}. The sum is {sums_by_k[k]}.")
    print("-" * 20)

    # Step 4: Compute the final result.
    total_sum = 0
    sum_terms = []
    for k in sorted(sums_by_k.keys()):
        total_sum += k * sums_by_k[k]
        sum_terms.append(f"{k}*({sums_by_k[k]})")

    final_equation = " + ".join(sum_terms)
    
    print("Step 4: Final Calculation")
    print(f"The formula is: \u03A3 k * (\u03A3 w) for k from 1 to 4.")
    print(f"This expands to: {final_equation}")
    print(f"Final value = {total_sum}")
    
    return total_sum

if __name__ == "__main__":
    solve_knot_grid_sum()