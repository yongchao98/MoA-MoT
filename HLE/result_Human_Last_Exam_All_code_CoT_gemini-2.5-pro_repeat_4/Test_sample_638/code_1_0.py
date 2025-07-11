import collections

def solve_knot_problem():
    """
    Solves the grid diagram winding number problem for the left-hand trefoil.
    """
    # 1. Establish the Grid Diagram (3x3 grid, 0-indexed)
    # o's on the anti-diagonal
    o_pos = {(0, 2), (1, 1), (2, 0)}
    # x's that satisfy the grid diagram rules
    x_pos = {(0, 0), (1, 2), (2, 1)}
    markers = o_pos.union(x_pos)

    # 3. Calculate the Disk Signs s(r,c)
    # Determine row and column orientations (x->o for rows, o->x for cols)
    # R0: x(0,0)->o(0,2) = right (+1)
    # R1: x(1,2)->o(1,1) = left (-1)
    # R2: x(2,1)->o(2,0) = left (-1)
    row_orient = {0: 1, 1: -1, 2: -1}
    # C0: o(2,0)->x(0,0) = up (+1)
    # C1: o(1,1)->x(2,1) = down (-1)
    # C2: o(0,2)->x(1,2) = down (-1)
    col_orient = {0: 1, 1: -1, 2: -1}

    # The left-hand trefoil signature is -2.
    # The sum of disk signs must be 1, since sigma = -2 * sum(s).
    # This corresponds to choosing the "white" squares (r+c is even) for disks.
    s = [[0] * 3 for _ in range(3)]
    sum_s = 0
    for r in range(3):
        for c in range(3):
            if (r + c) % 2 == 0:
                s[r][c] = row_orient[r] * col_orient[c]
                sum_s += s[r][c]

    # Verification: sum_s should be 1. It is.

    # 4. Compute the Winding Number Matrix w
    # w is a 4x4 matrix for vertices (0,0) to (3,3)
    w = [[0] * 4 for _ in range(4)]
    for i in range(1, 4):
        for j in range(1, 4):
            # The recursive formula for winding numbers
            w[i][j] = w[i - 1][j] + w[i][j - 1] - w[i - 1][j - 1] + s[i - 1][j - 1]

    # 5. Partition Winding Numbers by k
    # k is the number of marked squares adjacent to a vertex
    sums_by_k = collections.defaultdict(list)
    for i in range(4):
        for j in range(4):
            k = 0
            # A vertex (i,j) is a corner of up to 4 squares.
            # Check square (i-1, j-1)
            if i > 0 and j > 0 and (i - 1, j - 1) in markers:
                k += 1
            # Check square (i-1, j)
            if i > 0 and j < 3 and (i - 1, j) in markers:
                k += 1
            # Check square (i, j-1)
            if i < 3 and j > 0 and (i, j - 1) in markers:
                k += 1
            # Check square (i, j)
            if i < 3 and j < 3 and (i, j) in markers:
                k += 1
            if k > 0:
                sums_by_k[k].append(w[i][j])

    # 6. Calculate the Final Sum
    total_sum = 0
    calculation_str = []
    
    print("The winding number matrix w(i,j) is:")
    for row in w:
        print(" ".join(map(str, row)))
    print("-" * 20)
    
    for k in sorted(sums_by_k.keys()):
        sum_w_for_k = sum(sums_by_k[k])
        total_sum += k * sum_w_for_k
        
        # Format the numbers for the final equation string
        numbers_str = " + ".join(map(str, sums_by_k[k]))
        print(f"For k = {k}, the winding numbers are {sums_by_k[k]}.")
        print(f"Sum for k = {k}: {numbers_str} = {sum_w_for_k}")
        calculation_str.append(f"{k} * {sum_w_for_k}")
        
    print("-" * 20)
    print("The final calculation is:")
    final_expr_str = " + ".join(calculation_str)
    # Add k=4 term if it was missing
    if 4 not in sums_by_k:
         final_expr_str += " + 4 * 0"

    print(f"{final_expr_str} = {total_sum}")
    
    print("-" * 20)
    print(f"The value of the expression is {total_sum}")


solve_knot_problem()