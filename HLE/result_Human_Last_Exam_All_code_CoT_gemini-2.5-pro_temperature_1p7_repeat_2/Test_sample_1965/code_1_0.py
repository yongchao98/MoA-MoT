import collections
try:
    import pycuber as pc
except ImportError:
    print("This script requires the 'pycuber' library.")
    print("Please install it using: pip install pycuber")
    exit()

def solve_cube_permutations():
    """
    Calculates the number of 6-move permutations that solve the cube at step 4, 5, or 6.
    """
    # There are 12 basic 90-degree moves
    moves_str = ["R", "R'", "L", "L'", "U", "U'", "D", "D'", "F", "F'", "B", "B'"]
    moves = [pc.Formula(m) for m in moves_str]
    
    # The identity permutation represents a solved cube
    identity_perm = pc.Formula()

    # s_k_I[k] will store the number of k-move sequences returning to Identity
    s_k_I = {}

    # Use dynamic programming to count reachable states
    # `counts` maps a permutation to the number of ways to reach it
    counts = collections.defaultdict(int)
    counts[identity_perm] = 1  # Start with 1 way to be in the solved state (0 moves)

    # Iterate from 1 to 6 moves
    for k in range(1, 7):
        next_counts = collections.defaultdict(int)
        for perm, num in counts.items():
            for move in moves:
                # Apply next move
                next_perm = perm * move
                next_counts[next_perm] += num
        
        counts = next_counts
        # Store the number of ways to return to identity after k moves
        s_k_I[k] = counts.get(identity_perm, 0)
        
    s4_I = s_k_I.get(4, 0)
    s5_I = s_k_I.get(5, 0)
    s6_I = s_k_I.get(6, 0)

    # Final formula from Principle of Inclusion-Exclusion:
    # 132 * s_4(I) + 12 * s_5(I) + s_6(I)
    total_count = 132 * s4_I + 12 * s5_I + s6_I

    # Print the equation with the calculated values
    print(f"Number of 4-move sequences returning to solved state (s_4(I)): {s4_I}")
    print(f"Number of 5-move sequences returning to solved state (s_5(I)): {s5_I}")
    print(f"Number of 6-move sequences returning to solved state (s_6(I)): {s6_I}")
    print("\nThe total number of permutations is calculated by the formula: 132 * s_4(I) + 12 * s_5(I) + s_6(I)")
    print(f"Total = 132 * {s4_I} + 12 * {s5_I} + {s6_I} = {total_count}")
    print(f"\nFinal Answer: {total_count}")

solve_cube_permutations()