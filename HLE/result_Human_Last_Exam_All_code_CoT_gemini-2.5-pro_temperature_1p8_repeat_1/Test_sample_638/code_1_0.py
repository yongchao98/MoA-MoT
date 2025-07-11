import numpy as np

def solve_knot_problem():
    """
    Solves the problem by constructing the grid diagram, calculating Seifert surface winding numbers,
    and computing the specified weighted sum.
    """
    n = 3  # Grid size

    # 1. Define the marker matrix M for the left-hand trefoil knot.
    # O markers at (0,0), (1,1), (2,2) -> value -1
    # X markers at (1,0), (2,1), (0,2) -> value +1
    # Using [col][row] indexing
    M = np.zeros((n, n), dtype=int)
    M[0, 0] = -1  # O at (1,1)
    M[1, 1] = -1  # O at (2,2)
    M[2, 2] = -1  # O at (3,3)
    M[1, 0] = 1   # X at (2,1)
    M[2, 1] = 1   # X at (3,2)
    M[0, 2] = 1   # X at (1,3)

    # 2. Calculate the winding numbers w.
    # We use a padded grid for w to handle boundary conditions (w=0 outside).
    w_padded = np.zeros((n + 1, n + 1), dtype=int)
    for c in range(1, n + 1):
        for r in range(1, n + 1):
            w_padded[c, r] = (w_padded[c - 1, r] + 
                              w_padded[c, r - 1] - 
                              w_padded[c - 1, r - 1] + 
                              M[c - 1, r - 1])
    # Extract the 3x3 winding number matrix.
    w = w_padded[1:, 1:]
    
    # Let's map winding numbers to their respective k-sets.
    mho = {1: [], 2: [], 3: [], 4: []}
    
    # 3. Iterate through interior lattice points to populate mho_k.
    # Interior lattice points are (i,j) for i,j in {1, 2}.
    point_contributions = []
    for i in range(1, n):
        for j in range(1, n):
            # i,j are the 0-indexed interior lattice point coordinates
            
            # Count markers in the 4 surrounding cells.
            k = (abs(M[i - 1, j - 1]) + abs(M[i, j - 1]) + 
                 abs(M[i - 1, j]) + abs(M[i, j]))
            
            # Associate winding number of the bottom-left cell (i-1, j-1)
            # with this lattice point.
            w_val = w[i - 1, j - 1]
            
            if k > 0:
                mho[k].append(w_val)
                point_contributions.append({'k': k, 'w': w_val})

    # 4. Calculate and print the final sum.
    total_sum = 0
    sum_strings = []
    
    print("The final calculation is given by the formula:")
    print("S = sum_{k=1 to 4} (k * sum(w for w in Omega_k))")
    print("\nBased on the grid diagram for the left-hand trefoil knot, we find the following sets:")
    
    for k in range(1, 5):
        if not mho[k]:
            print(f"Omega_{k} = {{}}")
            continue
        
        sum_of_ws = sum(mho[k])
        total_sum += k * sum_of_ws
        
        w_strings = [str(val) for val in mho[k]]
        sum_string = f"{k} * ({' + '.join(w_strings)})"
        
        # Replace "+ -" with "- " for better readability
        sum_string = sum_string.replace('+ -', '- ')
        
        sum_strings.append(sum_string)
        
        print(f"Omega_{k} = {{{', '.join(map(str, mho[k]))}}}")

    final_equation_str = "S = " + " + ".join(sum_strings)
    final_equation_str = final_equation_str.replace('+ -', '- ')

    print("\nPlugging these into the formula, we get:")
    print(final_equation_str)
    
    final_calc_str = "S = "
    terms = []
    for k in range(1, 5):
        if mho[k]:
            term_val = k * sum(mho[k])
            terms.append(str(term_val))
    
    final_calc_str += " + ".join(terms).replace("+ -", "- ")
    print(final_calc_str)
    
    print(f"S = {total_sum}")
    print("\nThe final answer is:")
    print(f"<<<{total_sum}>>>")

solve_knot_problem()