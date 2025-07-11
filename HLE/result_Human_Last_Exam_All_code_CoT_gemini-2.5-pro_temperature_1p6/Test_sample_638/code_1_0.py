import collections

def solve_knot_problem():
    """
    Solves the problem by constructing the grid diagram, calculating winding numbers,
    and then computing the specified weighted sum.
    """
    n = 3
    # Using (row, col) coordinates, 0-indexed from top-left.
    # O markers are on the anti-diagonal.
    o_pos = {(0, 2), (1, 1), (2, 0)}
    # X markers are placed to form the left-hand trefoil with clockwise orientation.
    x_pos = {(0, 1), (1, 0), (2, 2)}
    all_markers = o_pos.union(x_pos)

    # --- Step 2: Calculate Winding Numbers ---
    # w[i][j] is the winding number at lattice point (i, j)
    # The grid of cells is n x n, so the grid of vertices is (n+1) x (n+1)
    w = [[0] * (n + 1) for _ in range(n + 1)]

    # We determine the change in w when crossing a strand.
    # Moving Right: +1 for UP strand, -1 for DOWN strand.
    # Moving Down: +1 for RIGHT strand, -1 for LEFT strand.

    # Vertical strands: path is O->X
    # V at x=0.5 (col 0): O(2,0) -> X(1,0) is UP. Crosses y=1.5. Affects crossing from y=1 to y=2. No, wait, crosses grid line at y=2.
    # Path is from cell (2,0) to (1,0). On x=0.5 line. Crosses y=1.5. So w changes along path y=2.
    # V at x=1.5 (col 1): O(1,1) -> X(0,1) is UP. Crosses y=0.5. w changes along path y=1.
    # V at x=2.5 (col 2): O(0,2) -> X(2,2) is DOWN. Crosses y=0.5, y=1.5. w changes along paths y=1, y=2.
    def get_dv(i, j_idx): # change crossing vertical line x=j_idx+0.5 at height i
        if j_idx == 0: # V at x=0.5 UP
            return 1 if i == 2 else 0
        if j_idx == 1: # V at x=1.5 UP
            return 1 if i == 1 else 0
        if j_idx == 2: # V at x=2.5 DOWN
            return -1 if i in [1, 2] else 0
        return 0

    # Horizontal strands: path is X->O
    # H at y=0.5 (row 0): X(0,1) -> O(0,2) is RIGHT. Crosses x=1.5. w changes along path x=2.
    # H at y=1.5 (row 1): X(1,0) -> O(1,1) is RIGHT. Crosses x=0.5. w changes along path x=1.
    # H at y=2.5 (row 2): X(2,2) -> O(2,0) is LEFT. Crosses x=1.5, x=0.5. w changes along paths x=1, x=2.
    def get_dh(i_idx, j): # change crossing horiz line y=i_idx+0.5 at position j
        if i_idx == 0: # H at y=0.5 RIGHT
            return 1 if j == 2 else 0
        if i_idx == 1: # H at y=1.5 RIGHT
            return 1 if j == 1 else 0
        if i_idx == 2: # H at y=2.5 LEFT
            return -1 if j in [1, 2] else 0
        return 0
    
    # We can fill the w matrix row by row.
    # w[i][j] can be computed from w[i-1][j] or w[i][j-1].
    # Let's fill using a wavefront approach starting from w[0][0]=0.
    for i in range(n + 1):
        for j in range(n + 1):
            if i == 0 and j == 0:
                continue
            if i > 0:
                w[i][j] = w[i - 1][j] + get_dh(i - 1, j)
            elif j > 0:
                w[i][j] = w[i][j - 1] + get_dv(i, j - 1)
    
    # --- Step 3: Calculate Corner Counts ---
    # k[i][j] is the number of markers for which vertex (i,j) is a corner
    k = [[0] * (n + 1) for _ in range(n + 1)]
    marked_cells = [[0] * n for _ in range(n)]
    for r, c in all_markers:
        marked_cells[r][c] = 1

    for i in range(n + 1):
        for j in range(n + 1):
            count = 0
            # Vertex (i,j) is a corner for cells (i-1,j-1), (i-1,j), (i,j-1), (i,j)
            if i > 0 and j > 0 and marked_cells[i - 1][j - 1]: count += 1
            if i > 0 and j < n and marked_cells[i - 1][j]: count += 1
            if i < n and j > 0 and marked_cells[i][j - 1]: count += 1
            if i < n and j < n and marked_cells[i][j]: count += 1
            k[i][j] = count
    
    # --- Step 4: Compute the Final Sum ---
    sums_by_k = collections.defaultdict(list)
    for i in range(n + 1):
        for j in range(n + 1):
            k_val = k[i][j]
            if k_val > 0:
                sums_by_k[k_val].append(w[i][j])
    
    print("This program calculates the sum based on the winding numbers of a Seifert surface on a grid diagram.")
    print("The grid diagram is the minimal one for the left-hand trefoil knot with clockwise orientation.")
    
    total_sum = 0
    equation_parts = []
    
    for k_val in range(1, 5):
        w_values = sums_by_k[k_val]
        sum_w = sum(w_values)
        term = k_val * sum_w
        total_sum += term
        
        print(f"\nFor k = {k_val}:")
        print(f"  The winding numbers are: {w_values}")
        print(f"  The sum of these winding numbers is: {sum_w}")
        print(f"  The contribution to the total sum is: {k_val} * {sum_w} = {term}")
        equation_parts.append(f"{k_val} * ({sum_w})")
    
    final_equation = " + ".join(equation_parts)
    print("\nThe final equation is:")
    print(f"{final_equation} = {total_sum}")

    print(f"\nThe value of the expression is {total_sum}.")
    return total_sum

if __name__ == "__main__":
    result = solve_knot_problem()
    print(f"\n<<<8>>>")