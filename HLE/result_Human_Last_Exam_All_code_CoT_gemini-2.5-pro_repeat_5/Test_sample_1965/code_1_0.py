import collections

def solve_rubiks_permutations():
    """
    This script calculates the number of 6-move sequences that solve a Rubik's cube
    at some point during the final 3 moves.

    It requires the 'pycuber' library. You can install it using:
    pip install pycuber
    """
    try:
        import pycuber as pc
    except ImportError:
        print("Error: The 'pycuber' library is required.")
        print("Please install it using: pip install pycuber")
        return

    # Define the 12 possible 90-degree face turns
    moves_str = "U U' D D' L L' R R' F F' B B'"
    moves = [pc.Move(s) for s in moves_str.split()]
    
    print("Starting computation... (This may take up to a minute)")

    # BFS to find N(k), the number of k-move sequences that solve the cube.
    # We use a cube's tuple representation as a hashable state identifier.
    solved_cube = pc.Cube()
    solved_id = solved_cube.f_tuple

    # counts: maps state_id -> number of paths to reach it
    counts = collections.defaultdict(int)
    counts[solved_id] = 1 # Start with 1 path of length 0 to the solved state

    # known_cubes: caches cube objects to avoid recreating them from their IDs
    known_cubes = {solved_id: solved_cube}

    N_values = {}
    max_k = 6

    for k in range(1, max_k + 1):
        next_counts = collections.defaultdict(int)
        for cube_id, count in counts.items():
            current_cube = known_cubes[cube_id]
            for move in moves:
                # The '+' operator on a cube and a move creates a new cube object
                next_cube = current_cube + move
                next_cube_id = next_cube.f_tuple
                
                # Cache the new cube object if it's the first time we've seen it
                if next_cube_id not in known_cubes:
                    known_cubes[next_cube_id] = next_cube
                
                # Add the number of paths from the previous state to the new state's total
                next_counts[next_cube_id] += count
        
        counts = next_counts
        # N(k) is the number of paths of length k that lead back to the solved state
        N_values[k] = counts[solved_id]
        print(f"Computed N({k}) = {N_values[k]}")

    N2 = N_values[2]
    N4 = N_values[4]
    N6 = N_values[6]

    # Based on the inclusion-exclusion principle, the total number of permutations is:
    # N(A) + N(C) - N(A intersect C)
    # = (N(4) * 12 * 12) + N(6) - (N(4) * N(2))
    # = 144 * N(4) + N(6) - N(2) * N(4)
    
    total_permutations = (144 * N4) + N6 - (N2 * N4)
    
    print("\n--- Calculation ---")
    print(f"Number of 2-move solving sequences: N(2) = {N2}")
    print(f"Number of 4-move solving sequences: N(4) = {N4}")
    print(f"Number of 6-move solving sequences: N(6) = {N6}")
    
    print("\nThe final formula is: 144 * N(4) + N(6) - N(2) * N(4)")
    print(f"Plugging in the numbers: 144 * {N4} + {N6} - {N2} * {N4}")
    print(f"= {144 * N4} + {N6} - {N2 * N4}")
    print(f"= {144 * N4 + N6} - {N2 * N4}")
    print(f"= {total_permutations}")
    
    print("\n--- Final Answer ---")
    print("The number of permutations that result in the cube returning to its original configuration is:")
    print(total_permutations)

if __name__ == '__main__':
    solve_rubiks_permutations()
<<<193248>>>