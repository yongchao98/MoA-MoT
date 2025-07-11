import collections

def solve_knot_problem():
    """
    Calculates the specified sum for the left-hand trefoil knot grid diagram.
    """
    n = 3  # Minimal grid number for the trefoil knot

    # o's on the main diagonal (bottom-left to top-right)
    o_pos = set((i, i) for i in range(n))
    # x's for the left-hand trefoil knot with o's at (i,i)
    x_pos = set((i, (i + 1) % n) for i in range(n))

    markers = o_pos.union(x_pos)

    # Step 1: Determine segment orientations
    # Path orientation: o -> V -> x, x -> H -> o gives a clockwise path
    # which is required.
    h_segments = {} # {row: orientation}
    v_segments = {} # {col: orientation}
    
    # Horizontal segments: connect x to o in the same row
    # Path is from x to o
    rows_with_markers = set(r for r,c in markers)
    for r in rows_with_markers:
        markers_in_row = [m for m in markers if m[0] == r]
        o_marker = [m for m in markers_in_row if m in o_pos][0]
        x_marker = [m for m in markers_in_row if m in x_pos][0]
        if x_marker[1] < o_marker[1]:
            h_segments[r] = "RIGHT"
        else:
            h_segments[r] = "LEFT"

    # Vertical segments: connect o to x in the same col
    # Path is from o to x
    cols_with_markers = set(c for r,c in markers)
    for c in cols_with_markers:
        markers_in_col = [m for m in markers if m[1] == c]
        o_marker = [m for m in markers_in_col if m in o_pos][0]
        x_marker = [m for m in markers_in_col if m in x_pos][0]
        if o_marker[0] < x_marker[0]:
            v_segments[c] = "UP"
        else:
            v_segments[c] = "DOWN"

    # Step 2: Calculate winding numbers w(i,j) for each lattice point (i,j)
    w_matrix = [[0] * (n + 1) for _ in range(n + 1)]
    
    # Contribution from horizontal segments
    w_r = [0] * (n + 1)
    for i in range(1, n + 1):
        # Contribution from segment in row i-1
        row = i - 1
        w_r[i] = w_r[i-1]
        if row in h_segments:
            orientation = h_segments[row]
            # Crossing from bottom to top
            if orientation == "LEFT": # is crossing from left to right of segment
                w_r[i] -= 1
            elif orientation == "RIGHT": # is crossing from right to left
                w_r[i] += 1

    # Contribution from vertical segments
    w_c = [0] * (n + 1)
    for j in range(1, n + 1):
        # Contribution from segment in column j-1
        col = j - 1
        w_c[j] = w_c[j-1]
        if col in v_segments:
            orientation = v_segments[col]
            # Crossing from left to right
            if orientation == "UP": # is crossing from left to right
                w_c[j] -= 1
            elif orientation == "DOWN": # is crossing from right to left
                w_c[j] += 1

    # Total winding number is the sum of contributions
    for i in range(n + 1):
        for j in range(n + 1):
            w_matrix[i][j] = w_r[i] + w_c[j]

    # Step 3: Count adjacent markers k for each lattice point
    k_matrix = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        for j in range(n + 1):
            count = 0
            # A lattice point (i,j) is a corner for cells (i-1,j-1), (i-1,j), (i,j-1), (i,j)
            for r_offset in [-1, 0]:
                for c_offset in [-1, 0]:
                    cell_pos = (i + r_offset, j + c_offset)
                    if cell_pos in markers:
                        count += 1
            k_matrix[i][j] = count
            
    # Step 4: Group winding numbers into mho_k
    mho = collections.defaultdict(list)
    for i in range(n + 1):
        for j in range(n + 1):
            k = k_matrix[i][j]
            w = w_matrix[i][j]
            if k > 0:
                mho[k].append(w)
                
    # Step 5: Calculate sums and final result
    total_sum = 0
    sum_by_k = {}
    equation_parts = []
    
    for k in range(1, 5):
        sum_for_k = sum(mho[k])
        sum_by_k[k] = sum_for_k
        total_sum += k * sum_for_k
        equation_parts.append(f"{k}*({sum_for_k})")
        
    equation = " + ".join(equation_parts)
    print(f"{equation} = {total_sum}")

solve_knot_problem()