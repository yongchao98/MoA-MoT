import collections

def solve_rubiks_problem():
    """
    Calculates the number of 6-move sequences that solve a Rubik's cube
    at some point during the final 3 moves.
    """
    try:
        import pycuber as pc
    except ImportError:
        print("This script requires the 'pycuber' library.")
        print("Please install it by running: pip install pycuber")
        return

    print("Step 1: Planning the calculation.")
    print("The goal is to find the number of 6-move sequences solving the cube after move 4, 5, or 6.")
    print("Using the Principle of Inclusion-Exclusion, the formula is:")
    print("Total = |S4| + |S5| + |S6| - |S4 ∩ S6|")
    print("where |S_k| is the set of sequences solved after k moves.")
    print("This simplifies to: (N(4) * 12^2) + (N(5) * 12) + N(6) - (N(4) * N(2))")
    print("\nStep 2: Calculating N(k), the number of k-move sequences that solve the cube.")
    
    # The 12 standard 90-degree moves
    moves_str = ["U", "U'", "D", "D'", "L", "L'", "R", "R'", "F", "F'", "B", "B'"]
    moves = [pc.Move(m) for m in moves_str]
    
    # Initial state: a solved cube
    identity_cube = pc.Cube()
    identity_tuple = identity_cube.to_tuple()
    
    # N[k] will store the number of k-move sequences returning to the solved state
    N = {}
    
    # level_counts maps each cube state (as a tuple) to the number of sequences that reach it
    level_counts = collections.defaultdict(int)
    level_counts[identity_tuple] = 1

    # Loop from k=1 to 6 to find N(k)
    for k in range(1, 7):
        next_level_counts = collections.defaultdict(int)
        for state_tuple, count in level_counts.items():
            # Create a cube from its tuple representation
            current_cube = pc.Cube(state_tuple=state_tuple)
            for move in moves:
                # Apply the move to a copy of the cube
                next_cube = current_cube.copy()
                next_cube.perform_move(move)
                next_cube_tuple = next_cube.to_tuple()
                # Add the number of ways to reach the previous state to the new state's count
                next_level_counts[next_cube_tuple] += count
        
        level_counts = next_level_counts
        # N(k) is the number of sequences that lead back to the identity state
        N[k] = level_counts.get(identity_tuple, 0)
        print(f"Computed N({k}) = {N[k]}")

    print("\nStep 3: Calculating the final answer using the formula.")

    # Retrieve computed values
    N2 = N.get(2, 0)
    N4 = N.get(4, 0)
    N5 = N.get(5, 0)
    N6 = N.get(6, 0)

    # Calculate each term of the inclusion-exclusion formula
    num_s4 = N4 * 144  # N(4) * 12^2
    num_s5 = N5 * 12   # N(5) * 12
    num_s6 = N6        # N(6)
    num_s4_intersect_s6 = N4 * N2 # N(4) * N(2)
    
    # Final calculation
    total_permutations = num_s4 + num_s5 + num_s6 - num_s4_intersect_s6

    print(f"Number of sequences solved after 4 moves |S4| = N(4) * 12^2 = {N4} * 144 = {num_s4}")
    print(f"Number of sequences solved after 5 moves |S5| = N(5) * 12 = {N5} * 12 = {num_s5}")
    print(f"Number of sequences solved after 6 moves |S6| = N(6) = {N6}")
    print(f"Intersection |S4 ∩ S6| = N(4) * N(2) = {N4} * {N2} = {num_s4_intersect_s6}")
    print("\nFinal Result Calculation:")
    print(f"Total = {num_s4} + {num_s5} + {num_s6} - {num_s4_intersect_s6} = {total_permutations}")
    
    print(f"\nOf the 2,985,984 possible permutations, {total_permutations} result in the cube returning to its original configuration at some point during the final 3 moves.")
    
    # The final answer in the required format
    print(f"\n<<<{total_permutations}>>>")

solve_rubiks_problem()