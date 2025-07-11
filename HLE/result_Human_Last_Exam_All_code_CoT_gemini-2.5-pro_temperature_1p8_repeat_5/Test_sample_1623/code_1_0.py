import numpy as np

def calculate_writhe(o_pos, x_pos_by_col):
    """Calculates the writhe of the knot from the grid diagram."""
    n = len(o_pos)
    
    # Create map from row to column for O and X markers
    o_pos_rev = {r: c for c, r in o_pos.items()}
    x_pos_by_row = {r: c for c, r in x_pos_by_col.items()}

    # Trace the path to find orientations
    h_dirs = {}
    v_dirs = {}
    
    # Path tracing starts at o(1,1)
    # The knot is one component: 1->3->5->2->4->1 for cols
    path_map = {1:3, 3:5, 5:2, 2:4, 4:1}
    
    for start_col, end_col in path_map.items():
        o_row = o_pos[start_col]
        x_row = o_pos[end_col]
        # horizontal segment from o(start_col) to x on same row
        x_col_horiz = x_pos_by_row[o_row]
        h_dirs[o_row] = 1 if x_col_horiz > start_col else -1

        # vertical segment from x(end_col) to o on same col
        v_dirs[end_col] = 1 if x_row > x_pos_by_col[end_col] else -1

    writhe = 0
    crossings = []
    for i in range(1, n + 1): # column for vertical segment
        for j in range(1, n + 1): # row for horizontal segment
            o_row_v = o_pos[i]
            x_row_v = x_pos_by_col[i]
            
            o_col_h = o_pos_rev[j]
            x_col_h = x_pos_by_row[j]
            
            # Check for crossing
            is_vert_segment = (j > min(o_row_v, x_row_v) and j < max(o_row_v, x_row_v))
            is_horiz_segment = (i > min(o_col_h, x_col_h) and i < max(o_col_h, x_col_h))

            if is_vert_segment and is_horiz_segment:
                sign = v_dirs[i] * h_dirs[j]
                writhe += sign
                crossings.append((i,j,sign))

    return writhe, crossings

def count_destabilizations(p_tuple):
    """Counts the number of possible destabilizations."""
    p = list(p_tuple)
    n = len(p)
    count = 0
    while n > 1:
        found_destab = False
        # Rule: p(i+1) = p(i) + 1 (using 1-based indexing)
        for i in range(n - 1):
            if p[i+1] == p[i] + 1:
                count += 1
                
                # Perform destabilization
                j = p[i]
                p.pop(i + 1)
                
                new_p = []
                for val in p:
                    if val > j:
                        new_p.append(val - 1)
                    else:
                        new_p.append(val)
                p = new_p
                n = len(p)
                found_destab = True
                break # Restart scan on the new smaller permutation
        if not found_destab:
            break
            
    return count

def main():
    n = 5
    o_positions = {i: i for i in range(1, n + 1)}
    x_positions_by_col = {1: 4, 2: 5, 3: 1, 4: 2, 5: 3}
    
    # The permutation p maps column i to the row of the X marker
    p_tuple = tuple(x_positions_by_col[i] for i in range(1, n + 1))
    
    writhe, crossings = calculate_writhe(o_positions, x_positions_by_col)
    
    print(f"Grid Size (n): {n}")
    print(f"Permutation p: {p_tuple}")
    
    print("\nStep 1: Calculate Writhe (w)")
    print(f"Found {len(crossings)} crossings:")
    for i, j, sign in crossings:
        print(f"  - At ({i},{j}) with sign {sign}")
    print(f"Total Writhe w = {writhe}")

    print("\nStep 2: Calculate Initial Thurston-Bennequin Number (tb)")
    tb_initial = writhe - n
    print(f"tb_initial = w - n = {writhe} - {n} = {tb_initial}")

    print("\nStep 3: Count Destabilizations")
    num_destabilizations = count_destabilizations(p_tuple)
    print(f"The grid can be destabilized {num_destabilizations} times.")

    print("\nStep 4: Calculate Maximal Thurston-Bennequin Number (TB_max)")
    tb_max = tb_initial + num_destabilizations
    print(f"TB_max = tb_initial + (number of destabilizations)")
    print(f"TB_max = {tb_initial} + {num_destabilizations} = {tb_max}")
    
    # Final answer format
    print(f"\nThe equation is: ({crossings[0][2]}) + ({crossings[1][2]}) + ({crossings[2][2]}) - {n} + {num_destabilizations} = {tb_max}")

main()