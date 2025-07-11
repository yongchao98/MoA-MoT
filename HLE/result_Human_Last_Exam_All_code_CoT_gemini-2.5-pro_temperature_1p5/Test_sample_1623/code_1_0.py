import numpy as np

def calculate_writhe_and_crossings(p):
    """
    Calculates the writhe for a grid diagram given by permutation p.
    Orientation convention: horizontal o->x, vertical x->o.
    Crossing rule: vertical over horizontal.
    """
    n = len(p)
    # o is at (j,j)
    x_coords_by_col = p
    
    # x_coords by row for horizontal segments
    x_coords_by_row = {v: k for k, v in x_coords_by_col.items()}

    # Directions of segments
    # Vertical (x->o): sign(i - p[i])
    v_dir = {i: np.sign(i - x_coords_by_col[i]) for i in range(1, n + 1)}
    # Horizontal (o->x): sign(x_coords_by_row[j] - j)
    h_dir = {j: np.sign(x_coords_by_row.get(j, j) - j) for j in range(1, n + 1)}

    writhe = 0
    crossings = []

    # Iterate over all potential crossing points (i, j)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if i not in v_dir or j not in h_dir: continue
            
            # Check for crossing existence
            o_y, x_y = i, x_coords_by_col[i]
            o_x, x_x = j, x_coords_by_row[j]

            if (j - min(o_y, x_y)) * (j - max(o_y, x_y)) < 0 and \
               (i - min(o_x, x_x)) * (i - max(o_x, x_x)) < 0:
                
                # Crossing exists, calculate sign
                # v_dir_sign is -1 for down, +1 for up
                # h_dir_sign is -1 for left, +1 for right
                # Using CCW rule for sign: +1 if turn from vertical to horizontal is CCW
                # v_up(+1), h_right(+1) -> sign -1 (CW)
                # v_up(+1), h_left(-1) -> sign +1 (CCW)
                # v_down(-1), h_right(+1) -> sign +1 (CCW)
                # v_down(-1), h_left(-1) -> sign -1 (CW)
                # This is equivalent to sign = -v_dir * h_dir
                sign = -v_dir[i] * h_dir[j]
                writhe += sign
                crossings.append(int(sign))

    return writhe, crossings

def count_components(p):
    """Counts the number of cycles in the permutation p."""
    if not p:
        return 0
    n = len(p)
    visited = set()
    cycles = 0
    for i in range(1, n + 1):
        if i not in visited:
            cycles += 1
            j = i
            while j not in visited:
                visited.add(j)
                j = p[j]
    return cycles

def destabilize(p, i_pos):
    """Performs a positive destabilization on permutation p at index i_pos."""
    n = len(p)
    
    # Identify items to remove and remap
    i_rem = i_pos
    p_rem = p[i_rem]
    
    old_indices = sorted(p.keys())
    old_indices.remove(i_rem)
    new_indices_map = {old: new for new, old in enumerate(old_indices, 1)}

    old_values = sorted(p.values())
    old_values.remove(p_rem)
    new_values_map = {old: new for new, old in enumerate(old_values, 1)}

    new_p = {}
    for old_idx in old_indices:
        new_idx = new_indices_map[old_idx]
        old_val = p[old_idx]
        new_val = new_values_map[old_val]
        new_p[new_idx] = new_val
        
    return new_p

def find_pos_destab(p):
    """Finds an index for positive destabilization."""
    n = len(p)
    for i in range(1, n):
        if p.get(i + 1) == p.get(i) + 1:
            return i
    return None

def solve():
    """
    Calculates the maximal Thurston-Bennequin number for the given grid diagram.
    """
    # Step 1: Define the grid diagram permutation
    # O's at (i, i) for i=1..5
    # X's at (1,4), (2,5), (3,1), (4,2), (5,3)
    # Permutation p maps column i to row of X in that column
    p_initial = {1: 4, 2: 5, 3: 1, 4: 2, 5: 3}

    # Step 2: Calculate initial writhe and number of components
    writhe, crossings = calculate_writhe_and_crossings(p_initial)
    num_components = count_components(p_initial)
    
    # Step 3: Calculate the initial Thurston-Bennequin number
    tb_initial = writhe - num_components

    print("Step 1: Calculate the writhe of the initial grid diagram.")
    crossings_str = [f"({c})" for c in crossings]
    print(f"The writhe is the sum of signs of all crossings.")
    print(f"w = {' + '.join(crossings_str)} = {writhe}")
    print("")

    print("Step 2: Calculate the initial Thurston-Bennequin number (tb).")
    print(f"The number of components, c, is the number of cycles in the permutation, which is {num_components}.")
    print(f"The formula for tb is: tb = w - c.")
    print(f"tb_initial = {writhe} - {num_components} = {tb_initial}")
    print("")

    # Step 4: Iteratively destabilize the grid to find the maximal tb
    num_destabilizations = 0
    tb = tb_initial
    current_p = p_initial.copy()
    
    print("Step 3: Maximize tb by destabilizing the grid diagram.")
    print("A grid can be destabilized if p(i+1) = p(i)+1 for some i (positive) or if its transpose has this property (negative). Each destabilization increases tb by 1.")

    while len(current_p) > 1:
        i_pos = find_pos_destab(current_p)
        if i_pos is not None:
            num_destabilizations += 1
            tb += 1
            print(f"Found a positive destabilization at i={i_pos}. Destabilizing... (tb increases to {tb})")
            current_p = destabilize(current_p, i_pos+1) # The rule simplifies i and i+1
        else:
            p_inv = {v: k for k, v in current_p.items()}
            i_neg = find_pos_destab(p_inv)
            if i_neg is not None:
                num_destabilizations += 1
                tb += 1
                print(f"Found a negative destabilization. Destabilizing... (tb increases to {tb})")
                p_inv_new = destabilize(p_inv, i_neg+1)
                current_p = {v: k for k,v in p_inv_new.items()}
            else:
                print("No more destabilizations found. The diagram is minimal.")
                break
    
    tb_max = tb_initial + num_destabilizations
    print("")
    print("Step 4: Final calculation.")
    print("The knot associated with the diagram is the one represented by the minimal grid.")
    print(f"Initial tb was {tb_initial}. There were {num_destabilizations} destabilizations.")
    print(f"Maximal tb = tb_initial + number_of_destabilizations")
    print(f"Maximal tb = {tb_initial} + {num_destabilizations} = {tb_max}")
    print(f"\nThe maximal Thurston-Bennequin number of the associated knot is {tb_max}.")
    print(f"<<<{tb_max}>>>")

solve()