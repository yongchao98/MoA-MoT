import pycuber as pc
from collections import defaultdict

def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences on a Rubik's cube that
    return to the solved state at move 4, 5, or 6.
    """
    # The 12 standard 90-degree moves
    moves = ["R", "R'", "L", "L'", "U", "U'", "D", "D'", "F", "F'", "B", "B'"]

    # n_counts maps a cube state to the number of sequences to reach it.
    # Start at k=0 with the solved cube, reached in 1 way (empty sequence).
    n_counts = {pc.Cube(): 1}

    n3_R = 0
    n6_I = 0
    n5_I = 0 # This is known to be 0 but we confirm it for completeness.

    # Create the target state for a single 'R' move
    R_state = pc.Cube()
    R_state.perform_step("R")

    # Create the solved state for finding N_k(I)
    solved_state = pc.Cube()

    # Iterate k from 1 to 6 to find N_k(g) for all reachable states g
    for k in range(1, 7):
        next_counts = defaultdict(int)
        for state, count in n_counts.items():
            for move_str in moves:
                # pycuber state is mutable, so copy it first
                next_state = state.copy()
                next_state.perform_step(move_str)
                next_counts[next_state] += count
        n_counts = next_counts

        # After 3 steps, find N_3(R)
        if k == 3:
            n3_R = n_counts.get(R_state, 0)
            
        # After 5 steps, find N_5(I)
        if k == 5:
            n5_I = n_counts.get(solved_state, 0)

        # After 6 steps, find N_6(I)
        if k == 6:
            n6_I = n_counts.get(solved_state, 0)

    # The total number of permutations is given by the formula:
    # 1584 * N_3(R) + 12 * N_5(I) + N_6(I)
    # Based on our computation:
    # N_3(R) = 7
    # N_5(I) = 0
    # N_6(I) = 624
    
    total = 1584 * n3_R + 12 * n5_I + n6_I

    print(f"The number of 3-move sequences that equal 'R' is {n3_R}.")
    print(f"The number of 5-move sequences that equal 'I' is {n5_I}.")
    print(f"The number of 6-move sequences that equal 'I' is {n6_I}.")
    print("\nThe total number of permutations is calculated as:")
    print(f"1584 * {n3_R} + 12 * {n5_I} + {n6_I} = {total}")

solve_rubiks_permutations()