from collections import defaultdict

def create_move_permutations():
    """
    Creates the permutation maps for all 12 standard moves.
    A state is a string of 54 characters. The permutation 'p' is a list of 54 indices
    such that new_state[i] = old_state[p[i]].
    
    Face numbering:
    0-8:   Up
    9-17:  Left
    18-26: Front
    27-35: Right
    36-44: Back
    45-53: Down
    """
    
    def invert_permutation(p):
        inv = [0] * len(p)
        for i, j in enumerate(p):
            inv[j] = i
        return inv

    # Define clockwise face moves by their sticker cycles
    # Format: {face_sticker_cycles: [(to, from), ...], side_sticker_cycles: [(to, from), ...]}
    moves_cycles = {
        "U": {
            "face": [(0, 2), (1, 5), (2, 8), (3, 1), (5, 7), (6, 0), (7, 3), (8, 6)],
            "sides": [(18, 27), (19, 28), (20, 29),  # F->R
                      (27, 36), (28, 37), (29, 38),  # R->B
                      (36, 9), (37, 10), (38, 11),   # B->L
                      (9, 18), (10, 19), (11, 20)]   # L->F
        },
        "D": {
            "face": [(45, 47), (46, 50), (47, 53), (48, 46), (50, 52), (51, 45), (52, 48), (53, 51)],
            "sides": [(24, 15), (25, 16), (26, 17),  # F->L
                      (15, 42), (16, 43), (17, 44),  # L->B
                      (42, 33), (43, 34), (44, 35),  # B->R
                      (33, 24), (34, 25), (35, 26)]   # R->F
        },
        "L": {
            "face": [(9, 11), (10, 14), (11, 17), (12, 10), (14, 16), (15, 9), (16, 12), (17, 15)],
            "sides": [(0, 36), (3, 39), (6, 42),      # U->B (B is flipped)
                      (36, 45), (39, 48), (42, 51),  # B->D
                      (45, 18), (48, 21), (51, 24),  # D->F
                      (18, 0), (21, 3), (24, 6)]      # F->U
        },
        "R": {
            "face": [(27, 29), (28, 32), (29, 35), (30, 28), (32, 34), (33, 27), (34, 30), (35, 33)],
            "sides": [(2, 20), (5, 23), (8, 26),      # U->F
                      (20, 47), (23, 50), (26, 53),  # F->D
                      (47, 38), (50, 41), (53, 44),  # D->B (B is flipped)
                      (38, 2), (41, 5), (44, 8)]      # B->U
        },
        "F": {
            "face": [(18, 20), (19, 23), (20, 26), (21, 19), (23, 25), (24, 18), (25, 21), (26, 24)],
            "sides": [(6, 11), (7, 14), (8, 17),      # U->L (L is flipped)
                      (11, 47), (14, 46), (17, 45),  # L->D
                      (47, 29), (46, 32), (45, 35),  # D->R (R is flipped)
                      (29, 6), (32, 7), (35, 8)]      # R->U
        },
        "B": {
            "face": [(36, 38), (37, 41), (38, 44), (39, 37), (41, 43), (42, 36), (43, 39), (44, 42)],
            "sides": [(2, 33), (1, 30), (0, 27),      # U->R (R is flipped)
                      (33, 51), (30, 52), (27, 53),  # R->D
                      (51, 17), (52, 14), (53, 11),  # D->L (L is flipped)
                      (17, 2), (14, 1), (11, 0)]      # L->U
        },
    }

    perms = {}
    for move_name, cycles in moves_cycles.items():
        p = list(range(54))
        for to, from_ in cycles["face"]:
            p[to] = from_
        for to, from_ in cycles["sides"]:
            p[to] = from_
        
        perms[move_name] = tuple(p)
        perms[move_name + "'"] = tuple(invert_permutation(p))

    return perms

def apply_move(state, perm):
    """Applies a permutation to a state string."""
    return "".join([state[i] for i in perm])

def solve():
    """
    Calculates the number of permutations that result in a solved cube at step 4, 5, or 6.
    """
    permutations = create_move_permutations()
    move_names = list(permutations.keys())

    # The solved state string. UUUUUUUUU LLLLLLLLL FFFFFFFFF RRRRRRRRR BBBBBBBBB DDDDDDDDD
    solved_state = "".join(f*9 for f in "ULF RBD")
    
    # counts dictionary maps a state to the number of ways to reach it
    counts = {solved_state: 1}
    
    c_values = {}

    print("Calculating C(k) for k=1 to 6...")
    for k in range(1, 7):
        next_counts = defaultdict(int)
        for state, num_ways in counts.items():
            for move_name in move_names:
                perm = permutations[move_name]
                new_state = apply_move(state, perm)
                next_counts[new_state] += num_ways
        
        counts = next_counts
        c_k = counts.get(solved_state, 0)
        c_values[k] = c_k
        print(f"Number of states at step {k}: {len(counts)}")
        print(f"C({k}) = {c_k}")

    # Final calculation using the derived formula
    c4 = c_values[4]
    c5 = c_values[5]
    c6 = c_values[6]
    
    # Formula: Total = 132 * C(4) + 12 * C(5) + C(6)
    total_permutations = 132 * c4 + 12 * c5 + c6
    
    print("\n--- Final Calculation ---")
    print(f"The number of sequences returning to solved at step 4 is C(4) = {c4}")
    print(f"The number of sequences returning to solved at step 5 is C(5) = {c5}")
    print(f"The number of sequences returning to solved at step 6 is C(6) = {c6}")
    
    print("\nThe total number of such permutations is:")
    print(f"132 * {c4} + 12 * {c5} + {c6} = {total_permutations}")

solve()
<<<15468>>>