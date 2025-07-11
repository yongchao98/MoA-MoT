# First, you need to install the pycuber library:
# pip install pycuber

import pycuber as pc
import collections

def count_solved_f2l(cube):
    """
    Checks the cube state and counts how many of the four F2L pairs are correctly solved.
    An F2L pair is solved if its white corner and corresponding edge are in the
    correct location and orientation in the first two layers (relative to the white face).
    """
    solved_count = 0
    
    # Pair 1: White-Green-Orange
    try:
        if (cube[("W", "G", "O")].pos == ("F", "D", "L") and cube[("W", "G", "O")].orientation == 0 and
            cube[("G", "O")].pos == ("F", "L") and cube[("G", "O")].orientation == 0):
            solved_count += 1
    except KeyError: pass # Piece not identifiable (only for malformed cubes)

    # Pair 2: White-Orange-Blue
    try:
        if (cube[("W", "O", "B")].pos == ("B", "D", "L") and cube[("W", "O", "B")].orientation == 0 and
            cube[("O", "B")].pos == ("B", "L") and cube[("O", "B")].orientation == 0):
            solved_count += 1
    except KeyError: pass
            
    # Pair 3: White-Blue-Red
    try:
        if (cube[("W", "B", "R")].pos == ("B", "D", "R") and cube[("W", "B", "R")].orientation == 0 and
            cube[("B", "R")].pos == ("B", "R") and cube[("B", "R")].orientation == 0):
            solved_count += 1
    except KeyError: pass

    # Pair 4: White-Red-Green
    try:
        if (cube[("W", "R", "G")].pos == ("F", "D", "R") and cube[("W", "R", "G")].orientation == 0 and
            cube[("R", "G")].pos == ("F", "R") and cube[("R", "G")].orientation == 0):
            solved_count += 1
    except KeyError: pass
            
    return solved_count

def find_f2l_solution():
    """
    Solves the problem by setting up the cube, defining the search parameters,
    and running a Breadth-First Search to find the optimal move sequence.
    """
    # The scramble given by Johnny
    scramble = pc.Formula("R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D")

    # Create a new cube and apply the scramble
    cube = pc.Cube()
    cube(scramble)

    # Moves are performed in the new orientation (Yellow top, Orange front).
    # We map these new moves to their equivalent moves in the standard (White top, Green front) orientation.
    move_map = {
        "U": pc.Formula("D"), "U'": pc.Formula("D'"), "U2": pc.Formula("D2"),
        "D": pc.Formula("U"), "D'": pc.Formula("U'"), "D2": pc.Formula("U2"),
        "F": pc.Formula("L"), "F'": pc.Formula("L'"), "F2": pc.Formula("L2"),
        "B": pc.Formula("R"), "B'": pc.Formula("R'"), "B2": pc.Formula("R2"),
        "L": pc.Formula("F"), "L'": pc.Formula("F'"), "L2": pc.Formula("F2"),
        "R": pc.Formula("B"), "R'": pc.Formula("B'"), "R2": pc.Formula("B2"),
    }
    all_new_moves = list(move_map.keys())

    # Initial state for BFS
    start_cube = cube

    # Check if 2 pairs are already solved (they aren't, but for completeness)
    if count_solved_f2l(start_cube) >= 2:
        print("Two F2L pairs are already solved before any moves.")
        print("Minimum length: 0")
        print("Sequence: (none)")
        print("\n<<< 0 >>>")
        return

    # BFS setup
    # The queue stores tuples of (cube_state, path_so_far)
    queue = collections.deque([(start_cube, [])]) 
    # Visited set stores string representations of cubes to prevent cycles
    visited = {str(start_cube)}

    print("Searching for the shortest sequence to solve two F2L pairs...")

    # Start BFS
    while queue:
        current_cube, path = queue.popleft()

        # Iterate through all possible next moves in the new orientation
        for new_move_str in all_new_moves:
            
            # Get the equivalent standard move from our map
            standard_move_formula = move_map[new_move_str]
            
            # Create the next cube state
            next_cube = current_cube.copy()
            next_cube(standard_move_formula)
            
            next_cube_str = str(next_cube)

            if next_cube_str not in visited:
                new_path = path + [new_move_str]
                
                # Check if this new state meets our goal
                if count_solved_f2l(next_cube) >= 2:
                    print(f"\nSuccess! Found a solution of length {len(new_path)}.")
                    print("The minimum length sequence of moves is:")
                    
                    final_sequence = " ".join(new_path)
                    # Print each character of the final sequence as requested
                    for char in final_sequence:
                        print(char, end="")
                    print()
                    
                    # Return the final answer in the specified format
                    print(f"\n<<< {len(new_path)} >>>")
                    return

                # If not a solution, add to the queue and visited set to explore later
                visited.add(next_cube_str)
                queue.append((next_cube, new_path))
    
    print("No solution found.")

if __name__ == '__main__':
    find_f2l_solution()
