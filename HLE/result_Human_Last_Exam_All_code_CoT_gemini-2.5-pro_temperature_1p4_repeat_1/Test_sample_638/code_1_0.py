import collections

def solve():
    """
    Solves the grid diagram problem for the left-hand trefoil knot.
    """
    # Step 1: Define the grid diagram for the left-hand trefoil knot
    # with o's on the main diagonal. n=3.
    # O-positions: (0,0), (1,1), (2,2) -> O_perm = [0, 1, 2]
    # X-positions: (0,1), (1,2), (2,0) -> X_perm = [1, 2, 0]
    n = 3
    O_perm = [0, 1, 2]
    X_perm = [1, 2, 0]
    print("Grid Diagram Setup:")
    print(f"O markers (col, row): {list(enumerate(O_perm))}")
    print(f"X markers (col, row): {list(enumerate(X_perm))}\n")
    
    # Step 2: Calculate winding numbers
    # Sign of a segment: sgn(end - start)
    # V_c sign: sgn(X_c - O_c)
    # H_r sign (direction): right is +1, left is -1
    # O_inv_perm and X_inv_perm map row to column
    O_inv_perm = [p[0] for p in sorted(enumerate(O_perm), key=lambda x: x[1])]
    X_inv_perm = [p[0] for p in sorted(enumerate(X_perm), key=lambda x: x[1])]

    def sgn(x):
        return 1 if x > 0 else -1 if x < 0 else 0

    V_signs = [sgn(X_perm[c] - O_perm[c]) for c in range(n)]
    H_dir_signs = [sgn(X_inv_perm[r] - O_inv_perm[r]) for r in range(n)]

    w = {}
    w[-1,-1] = 0 # Winding number for the unbounded region

    # Calculate winding numbers for all regions using the rules:
    # w_new = w_old + sign(V) when crossing a vertical segment from left
    # w_new = w_old - sign(H_dir) when crossing a horizontal segment from below
    for i in range(-1, n):
        for j in range(-1, n):
            if j > -1 and (i,j-1) in w:
                w[i,j] = w[i,j-1] - H_dir_signs[j]
            elif i > -1 and (i-1,j) in w:
                w[i,j] = w[i-1,j] + V_signs[i]
    
    # Extract the winding numbers for the four central regions
    w_values = {(i, j): w[i,j] for i in range(n-1) for j in range(n-1)}
    print("Calculated winding numbers w(i,j):")
    for (i,j), val in w_values.items():
        print(f"  w({i},{j}) = {val}")
    print()

    # Step 3: Classify winding numbers into sets Mho_k
    markers = set(enumerate(O_perm)) | set(enumerate(X_perm))
    mho = collections.defaultdict(list)

    print("Classifying winding numbers:")
    for (r_i, r_j), w_val in w_values.items():
        corners = [(r_i, r_j), (r_i + 1, r_j), (r_i, r_j + 1), (r_i + 1, r_j + 1)]
        k = sum(1 for c in corners if c in markers)
        mho[k].append(w_val)
        print(f"  Region ({r_i},{r_j}) with w={w_val} has k={k} markers at its corners.")
    print()

    # Step 4: Compute the final sum
    total_sum = 0
    sum_str_parts = []
    
    for k in range(1, 5):
        w_list = mho[k]
        if not w_list:
            sum_str_parts.append(f"({k} * 0)")
            continue
        
        current_sum = sum(w_list)
        total_sum += k * current_sum
        
        w_sum_str = " + ".join(map(str, w_list))
        sum_str_parts.append(f"({k} * ({w_sum_str}))")
        
    final_equation = " + ".join(sum_str_parts)
    print("The final calculation is:")
    print(f"{final_equation} = {total_sum}")

solve()