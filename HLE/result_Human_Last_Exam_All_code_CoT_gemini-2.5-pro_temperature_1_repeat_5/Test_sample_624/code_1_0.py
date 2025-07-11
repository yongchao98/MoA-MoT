import pycuber as pc
from collections import deque

def solve_f2l_pairs():
    """
    This function sets up the cube, solves for two F2L pairs, and prints the result.
    """
    # The scramble sequence provided by the user.
    scramble_str = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"

    # Create a pristine cube. By default, pycuber is oriented with White on top (U) and Green on front (F).
    cube = pc.Cube()

    # Apply the scramble to the cube.
    cube(scramble_str)

    # Re-orient the cube to have Yellow on top and Orange on the front.
    # This corresponds to an x2 y' rotation from the standard orientation.
    # After this, U is Yellow, D is White, F is Orange, R is Blue, L is Green, B is Red.
    cube("x2 y'")

    # Define the four F2L pairs by the colors of their stickers.
    # These are the "identities" of the pieces we are looking for.
    # The corner faces are specified by their axis (e.g., D-face is White, F-face is Orange).
    F2L_PAIRS = {
        # Front-Right Slot (Orange-Blue)
        "FR": (pc.Corner(D='W', F='O', R='B'), pc.Edge(F='O', R='B')),
        # Front-Left Slot (Orange-Green)
        "FL": (pc.Corner(D='W', F='O', L='G'), pc.Edge(F='O', L='G')),
        # Back-Right Slot (Red-Blue)
        "BR": (pc.Corner(D='W', B='R', R='B'), pc.Edge(B='R', R='B')),
        # Back-Left Slot (Red-Green)
        "BL": (pc.Corner(D='W', B='R', L='G'), pc.Edge(B='R', L='G')),
    }

    # Define the locations where the F2L pairs should be to be considered solved.
    F2L_SLOTS = {
        "FR": ("DFR", "FR"),
        "FL": ("DFL", "FL"),
        "BR": ("DBR", "BR"),
        "BL": ("DBL", "BL"),
    }

    def count_solved_f2l_pairs(c):
        """ Checks a cube state and returns the number of solved F2L pairs. """
        solved_count = 0
        for slot_name, (corner_identity, edge_identity) in F2L_PAIRS.items():
            corner_pos, edge_pos = F2L_SLOTS[slot_name]

            # Get the pieces currently in the target slots.
            corner_piece_in_place = c[corner_pos]
            edge_piece_in_place = c[edge_pos]

            # pycuber's Cubie equality '==' checks for both correct piece and correct orientation.
            if corner_piece_in_place == corner_identity and edge_piece_in_place == edge_identity:
                solved_count += 1
        return solved_count

    # --- Breadth-First Search (BFS) to find the shortest solution ---

    # A queue to hold states to visit: (cube_object, path_of_moves)
    queue = deque([(cube, [])])
    # A set to keep track of visited cube states to avoid cycles and redundant work.
    visited = {str(cube)}
    
    # Check if the initial state already has 2+ pairs solved.
    initial_solved_count = count_solved_f2l_pairs(cube)
    if initial_solved_count >= 2:
        print("The cube already has 2 or more F2L pairs solved.")
        print("The shortest sequence is empty.")
        print("The length of this sequence is: 0")
        print("<<<0>>>")
        return

    # Define all possible moves to explore.
    all_moves = "U U' U2 R R' R2 F F' F2 D D' D2 L L' L2 B B' B2".split()

    solution_path = None
    
    # Start the search.
    while queue:
        current_cube, path = queue.popleft()

        # Stop searching if paths get too long (optimization). F2L solutions are usually short.
        if len(path) >= 10:
            continue
            
        # Explore all possible next moves.
        for move_str in all_moves:
            next_cube = current_cube.copy()
            next_cube(move_str)
            
            # Check if this state has been visited.
            next_cube_id = str(next_cube)
            if next_cube_id not in visited:
                visited.add(next_cube_id)
                new_path = path + [move_str]

                # Check if the new state is the goal state.
                if count_solved_f2l_pairs(next_cube) >= 2:
                    solution_path = new_path
                    break # Found the shortest solution
                
                queue.append((next_cube, new_path))
        
        if solution_path:
            break

    # Print the final result.
    if solution_path:
        solution_str = " ".join(solution_path)
        solution_len = len(solution_path)
        print(f"The shortest sequence to solve two F2L pairs is: {solution_str}")
        print(f"The length of this sequence is: {solution_len}")
        print(f"<<<{solution_len}>>>")
    else:
        print("Could not find a solution within the search depth.")

# Run the solver.
solve_f2l_pairs()