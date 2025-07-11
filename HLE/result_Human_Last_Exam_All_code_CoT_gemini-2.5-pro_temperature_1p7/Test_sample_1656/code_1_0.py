def get_permutations_from_coords(o_coords, x_coords, n):
    """Create permutation dictionaries from coordinates."""
    p_o = {i: 0 for i in range(1, n + 1)}
    p_x = {i: 0 for i in range(1, n + 1)}
    for col, row in o_coords:
        p_o[col] = row
    for col, row in x_coords:
        p_x[col] = row
    return p_o, p_x

def invert_permutation(p):
    """Invert a permutation dictionary."""
    return {v: k for k, v in p.items()}

def compose_permutations(p2, p1):
    """Compose two permutations (p2 o p1)."""
    return {k: p2[p1[k]] for k in p1}

def count_cycles(p, n):
    """Count the number of cycles in a permutation."""
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

def count_extrema(p, n):
    """Count local minima and maxima in a permutation sequence, cyclically."""
    extrema_count = 0
    for i in range(1, n + 1):
        # Get previous and next indices with wrap-around
        prev_i = (i - 2 + n) % n + 1
        next_i = i % n + 1
        
        y_prev = p[prev_i]
        y_curr = p[i]
        y_next = p[next_i]
        
        # Check for local max or min
        is_max = y_prev < y_curr and y_curr > y_next
        is_min = y_prev > y_curr and y_curr < y_next
        
        if is_max or is_min:
            extrema_count += 1
            
    return extrema_count

def main():
    """
    Main function to calculate the braid index from grid diagram data.
    """
    grid_number = 7
    
    # O markers at (column, row)
    o_coords = [(1,1), (2,7), (3,4), (4,5), (5,3), (6,6), (7,2)]
    
    # X markers at (column, row)
    x_coords = [(1,2), (2,6), (3,3), (4,1), (5,7), (6,5), (7,4)]

    # 1. Define permutations p_o and p_x
    p_o, p_x = get_permutations_from_coords(o_coords, x_coords, grid_number)

    # 2. Calculate the number of Seifert circles (s)
    p_x_inv = invert_permutation(p_x)
    p_seifert = compose_permutations(p_x_inv, p_o)
    s = count_cycles(p_seifert, grid_number)
    
    # 3. Count vertical extrema (n_o, n_x)
    n_o = count_extrema(p_o, grid_number)
    n_x = count_extrema(p_x, grid_number)
    
    # 4. Calculate the braid index
    braid_index = n_o + n_x - s
    
    print(f"The number of Seifert circles (s) is: {s}")
    print(f"The number of 'O' extrema (n_O) is: {n_o}")
    print(f"The number of 'X' extrema (n_X) is: {n_x}")
    print("The braid index is calculated by the formula: b = n_O + n_X - s")
    print(f"Braid Index = {n_o} + {n_x} - {s} = {braid_index}")

if __name__ == "__main__":
    main()