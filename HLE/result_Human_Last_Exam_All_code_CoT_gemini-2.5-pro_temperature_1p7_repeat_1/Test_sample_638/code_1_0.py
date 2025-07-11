import numpy as np

def solve_knot_problem():
    """
    Solves the problem based on the grid diagram of a left-hand trefoil knot.
    """
    n = 3
    # Step 1: Define the grid based on the problem description.
    # O's are on the anti-diagonal: (0,2), (1,1), (2,0)
    # X's for the left-hand trefoil are at (0,1), (1,0), (2,2)
    o_pos = [(0, 2), (1, 1), (2, 0)]
    x_pos = [(0, 1), (1, 0), (2, 2)]
    
    # Step 2: Calculate winding numbers w(i,j).
    # Create the 'c' matrix.
    c = np.zeros((n, n), dtype=int)
    for r, col in o_pos:
        c[r, col] = -1
    for r, col in x_pos:
        c[r, col] = 1

    # Compute the 'w' matrix of winding numbers.
    w = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            w_im1_j = w[i - 1, j] if i > 0 else 0
            w_i_jm1 = w[i, j - 1] if j > 0 else 0
            w_im1_jm1 = w[i - 1, j - 1] if i > 0 and j > 0 else 0
            w[i, j] = w_im1_j + w_i_jm1 - w_im1_jm1 + c[i, j]
            
    # Step 3: Classify winding numbers into sets mho_k.
    all_markers = set(o_pos + x_pos)
    mho = {k: [] for k in range(5)} # k=0 to 4
    
    for i in range(n):
        for j in range(n):
            corners = [(i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)]
            k = sum(1 for corner in corners if corner in all_markers)
            mho[k].append(w[i, j])

    # Step 4: Compute the final sum and display the calculation.
    total_sum = 0
    sum_strings = []
    
    print("The winding numbers w(i,j) are:")
    print(w)
    print("\nThe sets of winding numbers mho_k are:")
    for k in range(1, 5):
        print(f"mho_{k}: {mho[k]}")

    print("\nThe calculation for the final sum is:")
    
    for k in range(1, 5):
        sum_k = sum(mho[k])
        total_sum += k * sum_k
        
        # Format the sum components for readable output
        if not mho[k]:
            sum_str_k = "(0)"
        else:
            sum_str_k = "(" + "+".join(map(str, mho[k])).replace("+-", "-") + ")"
            
        term_str = f"{k} * {sum_str_k}"
        sum_strings.append(term_str)
        
    final_equation = " + ".join(sum_strings)
    print(f"{final_equation} = {total_sum}")
    
solve_knot_problem()