import pycuber as pc
import collections

def solve():
    """
    Calculates the number of 6-move sequences that solve the cube during the last 3 moves.
    """
    # Define the 12 basic 90-degree moves
    move_strings = ["R", "R'", "L", "L'", "U", "U'", "D", "D'", "F", "F'", "B", "B'"]
    basic_moves = [pc.Formula(s) for s in move_strings]

    # M_perms is the set of 12 cube states resulting from a single 90-degree move
    # from the solved state.
    solved_cube = pc.Cube()
    M_perms = set()
    for move in basic_moves:
        c = pc.Cube()
        c(move)
        M_perms.add(c)

    # counts[n] will be a dictionary mapping a cube state to the number of
    # n-move sequences that produce that state.
    # We start with the solved cube at n=0, reachable in 1 way (the empty sequence).
    counts = {0: collections.defaultdict(int, {solved_cube: 1})}
    
    # We need to compute up to k_5, so we need to explore paths up to length 5.
    for n in range(5):
        counts[n+1] = collections.defaultdict(int)
        # For each state reachable in n moves, explore all 12 next moves.
        for cube_state, num_sequences in counts[n].items():
            for move in basic_moves:
                next_cube = cube_state.copy()
                next_cube(move)
                # Add the number of paths to the previous state to the new state's count.
                counts[n+1][next_cube] += num_sequences

    # Now calculate k_n for n=3, 4, 5.
    # k_n is the number of n-move sequences whose final state is in M_perms.
    k_values = {}
    for n in range(1, 6):
        k_n = 0
        if n in counts:
            for p in M_perms:
                k_n += counts[n][p]
        k_values[n] = k_n
        
    k3 = k_values.get(3, 0)
    k4 = k_values.get(4, 0)
    k5 = k_values.get(5, 0)

    # Final calculation based on the derived formula.
    total = 132 * k3 + 12 * k4 + k5
    
    print("Step 1: The problem is solved using the formula: Total = 132 * k_3 + 12 * k_4 + k_5")
    print("Step 2: Calculate the values of k_n, the number of n-move sequences that equal a single move.")
    print(f"k_3 (computed) = {k3}")
    print(f"k_4 (computed) = {k4}")
    print(f"k_5 (computed) = {k5}")
    print("\nStep 3: Substitute these values into the formula.")
    print(f"Total permutations = 132 * {k3} + 12 * {k4} + {k5}")
    print(f"Total permutations = {132 * k3} + {12 * k4} + {k5}")
    print(f"Total permutations = {total}")

solve()
<<<28428>>>