import collections

# Attempt to import pycuber, providing instructions if not found.
try:
    import pycuber
except ImportError:
    print("This script requires the 'pycuber' library.")
    print("Please install it by running: pip install pycuber")
    exit()

def solve_f2l_pairs():
    """
    This function calculates the shortest move sequence to solve two F2L pairs
    from a given scrambled state.
    """
    # Step 1: Define the scramble and reorientation
    # The scramble sequence provided by the user.
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    # The reorientation from (White-Top, Green-Front) to (Yellow-Top, Orange-Front) is x2 y.
    reorient_moves = "x2 y"

    # Step 2: Set up the initial and target states.
    # We use a standard WCA color scheme: U-White, F-Green, R-Red, L-Orange, D-Yellow, B-Blue
    # After reorientation (x2 y), the cube faces will be:
    # U: Yellow, D: White, F: Orange, B: Red, R: Blue, L: Green
    
    # Create the scrambled cube
    my_cube = pycuber.Cube()
    my_cube(scramble)

    # Verify the white cross is solved on the bottom (D-face) as a sanity check.
    # Note: A full check is complex; we trust the problem statement. The code proceeds regardless.
    # print("Initial cube after scramble (White Top, Green Front):")
    # print(my_cube)
    
    # Re-orient the cube to match how Johnny picks it up.
    my_cube(reorient_moves)

    # Create a "target" solved cube in the new orientation to identify correct piece colors.
    target_cube = pycuber.Cube()
    target_cube(reorient_moves)
    
    # Map face colors for the new orientation.
    colors = {
        "F": target_cube.get_face("F")[1][1].color,  # Orange
        "B": target_cube.get_face("B")[1][1].color,  # Red
        "R": target_cube.get_face("R")[1][1].color,  # Blue
        "L": target_cube.get_face("L")[1][1].color,  # Green
        "D": target_cube.get_face("D")[1][1].color   # White
    }

    # Step 3: Define the goal-checking function
    def check_goal(cube):
        """
        Checks if at least two F2L pairs are solved.
        An F2L pair is solved if its edge and corner are in the correct
        position with the correct orientation in the bottom two layers.
        """
        solved_pairs = 0
        
        # Pair 1: Front-Right (Orange-Blue edge, White-Orange-Blue corner)
        if (cube.get_face("F")[1][2].color == colors["F"] and
            cube.get_face("R")[1][0].color == colors["R"] and
            cube.get_face("F")[2][2].color == colors["F"] and
            cube.get_face("R")[2][0].color == colors["R"] and
            cube.get_face("D")[2][2].color == colors["D"]):
            solved_pairs += 1
            
        # Pair 2: Front-Left (Orange-Green edge, White-Orange-Green corner)
        if (cube.get_face("F")[1][0].color == colors["F"] and
            cube.get_face("L")[1][2].color == colors["L"] and
            cube.get_face("F")[2][0].color == colors["F"] and
            cube.get_face("L")[2][2].color == colors["L"] and
            cube.get_face("D")[2][0].color == colors["D"]):
            solved_pairs += 1
            
        # Pair 3: Back-Right (Red-Blue edge, White-Red-Blue corner)
        if (cube.get_face("B")[1][0].color == colors["B"] and
            cube.get_face("R")[1][2].color == colors["R"] and
            cube.get_face("B")[2][0].color == colors["B"] and
            cube.get_face("R")[2][2].color == colors["R"] and
            cube.get_face("D")[0][2].color == colors["D"]):
            solved_pairs += 1
            
        # Pair 4: Back-Left (Red-Green edge, White-Red-Green corner)
        if (cube.get_face("B")[1][2].color == colors["B"] and
            cube.get_face("L")[1][0].color == colors["L"] and
            cube.get_face("B")[2][2].color == colors["B"] and
            cube.get_face("L")[2][0].color == colors["L"] and
            cube.get_face("D")[0][0].color == colors["D"]):
            solved_pairs += 1

        return solved_pairs >= 2

    # Step 4: Perform Breadth-First Search (BFS)
    # The set of all possible moves.
    moves = ["U", "U'", "U2", "D", "D'", "D2", "F", "F'", "F2",
             "B", "B'", "B2", "R", "R'", "R2", "L", "L'", "L2"]
             
    # A mapping to prevent redundant moves (e.g., U after U', R after R).
    # Key is the face, value is the set of move prefixes on that face.
    face_map = {move[0] for move in moves}
    
    # Queue for BFS: stores tuples of (cube_state, path_so_far)
    queue = collections.deque([(my_cube, "")])
    # Visited set to avoid cycles and redundant computations. Stores string representation of cubes.
    visited = {str(my_cube)}
    
    print("Searching for the shortest solution...")

    while queue:
        current_cube, path = queue.popleft()
        
        # Check if the current state is the goal
        if check_goal(current_cube):
            path = path.strip()
            length = len(path.split())
            print("\n--- Solution Found ---")
            print(f"The minimum number of moves to solve two F2L pairs is: {length}")
            # The problem asks to output the numbers in the final equation,
            # which we interpret as printing the final sequence.
            print(f"The sequence is: {path}")
            print("--- End of Solution ---")
            
            # Return final answer in specified format.
            print(f"\n<<<{length}>>>")
            return

        # Generate next possible states
        last_move_face = None
        if path:
            last_move_face = path.split()[-1][0]
            
        for move in moves:
            # Pruning: Don't turn the same face twice in a row (e.g., R R' is redundant).
            if last_move_face and move[0] == last_move_face:
                continue
                
            next_cube = current_cube.copy()
            next_cube(move)
            
            if str(next_cube) not in visited:
                visited.add(str(next_cube))
                new_path = f"{path} {move}"
                queue.append((next_cube, new_path))
    
    print("No solution was found. This should not happen with a valid cube state.")

if __name__ == '__main__':
    solve_f2l_pairs()