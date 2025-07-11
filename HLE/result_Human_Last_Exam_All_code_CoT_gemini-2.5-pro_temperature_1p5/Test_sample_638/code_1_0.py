import collections

def solve_knot_sum():
    """
    Calculates the sum based on the Seifert surface of a left-hand trefoil knot.
    """
    # Grid markers: O at (0,0),(1,1),(2,2), X at (0,2),(1,0),(2,1)
    # Using (row, col) coordinates, bottom-left origin.
    markers = {
        (0, 0): 'O', (1, 1): 'O', (2, 2): 'O',
        (0, 2): 'X', (1, 0): 'X', (2, 1): 'X'
    }

    # Winding numbers w(i,j) for lattice points (col, row)
    # Based on standard Seifert surface construction for this knot diagram.
    windings = [[0] * 4 for _ in range(4)]
    windings[1][1] = 1
    windings[1][2] = 1
    windings[2][1] = 1
    windings[2][2] = 2

    w_for_k = collections.defaultdict(list)

    # Iterate over the 4x4 lattice points (c, r) = (col, row)
    for r in range(4):
        for c in range(4):
            # Calculate k, the number of adjacent markers
            k = 0
            # A point (c,r) is a corner for cells C(c-1,r-1), C(c,r-1), C(c-1,r), C(c,r).
            # Cell coordinates are from 0 to 2.
            potential_cells = [
                (r - 1, c - 1), (r - 1, c),
                (r, c - 1), (r, c)
            ]
            for cell_r, cell_c in potential_cells:
                if 0 <= cell_r <= 2 and 0 <= cell_c <= 2:
                    if (cell_r, cell_c) in markers:
                        k += 1
            
            w = windings[r][c]
            if w != 0:
                w_for_k[k].append(w)

    total_sum = 0
    equation_parts = []
    for k in sorted(w_for_k.keys()):
        if not w_for_k[k]:
             w_for_k[k].append(0) # For displaying "k * (0)"
        
    for k in range(1, 5):
        sum_w_k = sum(w_for_k[k])
        total_sum += k * sum_w_k
        
        # Format the sum part of the equation
        sum_str = "+".join(map(str, w_for_k[k])) if w_for_k[k] else "0"
        if len(w_for_k[k]) > 1:
            sum_str = f"({sum_str})"
            
        equation_parts.append(f"{k} * {sum_str}")

    final_equation = " + ".join(equation_parts)
    print("The final sum is calculated as:")
    print(final_equation)
    
    print("\nWhich evaluates to:")
    eval_parts = []
    for k in range(1,5):
        eval_parts.append(f"{k} * {sum(w_for_k[k])}")
    print(" + ".join(eval_parts) + f" = {total_sum}")
    
    # Return the final numeric answer as requested.
    # print(f"\nFinal Answer: {total_sum}")


solve_knot_sum()