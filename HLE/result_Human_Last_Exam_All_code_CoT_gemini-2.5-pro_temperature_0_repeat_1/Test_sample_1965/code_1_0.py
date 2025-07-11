import collections

def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that solve a Rubik's cube
    at some point during the final 3 moves.

    This script requires the 'pycuber' library.
    Install it using: pip install pycuber
    """
    try:
        import pycuber as pc
    except ImportError:
        print("Error: The 'pycuber' library is required.")
        print("Please install it using: pip install pycuber")
        return

    # The 12 standard 90-degree moves
    moves = ["U", "U'", "D", "D'", "L", "L'", "R", "R'", "F", "F'", "B", "B'"]

    # Get the string representation of a solved cube
    identity_cube = pc.Cube()
    identity_string = str(identity_cube)

    # Dictionary to store the number of k-move sequences returning to identity
    N = {}

    # We start at k=0 with one "path" to the identity state (the empty sequence)
    # paths[state_string] = number of ways to reach this state
    paths = collections.defaultdict(int)
    paths[identity_string] = 1
    
    print("Calculating Nk: the number of k-move sequences that return to the solved state...")

    # Calculate N_k for k from 1 to 6
    for k in range(1, 7):
        next_paths = collections.defaultdict(int)
        # For each state reachable in k-1 moves, explore the next 12 moves
        for state_str, count in paths.items():
            # It's faster to create the cube from string once per state
            current_cube = pc.Cube().from_string(state_str)
            for move in moves:
                # Apply the move to a copy to not alter the base for the next move
                next_cube = current_cube.copy()
                next_cube(move)
                next_paths[str(next_cube)] += count
        
        paths = next_paths
        N[k] = paths.get(identity_string, 0)
        print(f"N_{k} = {N[k]}")

    N_2 = N[2]
    N_4 = N[4]
    N_6 = N[6]

    # Using the simplified formula from the plan: Total = |A| + |C| - |A n C|
    # |A| = N_4 * 12 * 12
    # |C| = N_6
    # |A n C| = N_4 * N_2
    # Total = (N_4 * 144) + N_6 - (N_4 * N_2)
    # Total = N_4 * (144 - N_2) + N_6
    
    term1 = N_4 * (144 - N_2)
    result = term1 + N_6

    print("\nBased on the Principle of Inclusion-Exclusion, the formula is:")
    print("Total = N_4 * (144 - N_2) + N_6")
    print("\nSubstituting the calculated values:")
    print(f"Total = {N_4} * (144 - {N_2}) + {N_6}")
    print(f"Total = {N_4} * {144 - N_2} + {N_6}")
    print(f"Total = {term1} + {N_6}")
    print(f"Total = {result}")

solve_rubiks_permutations()