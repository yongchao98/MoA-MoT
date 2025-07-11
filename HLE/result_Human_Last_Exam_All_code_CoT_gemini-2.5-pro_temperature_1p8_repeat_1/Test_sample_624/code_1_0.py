import sys
from collections import deque

try:
    import pycuber as pc
except ImportError:
    print("Error: This script requires the 'pycuber' library.")
    print("Please run 'pip install pycuber' to install it and then run the script again.")
    sys.exit(1)

def is_f2l_goal_state(cube):
    """
    Checks if any of the three target F2L pairs are solved.
    We assume the White-Orange-Green pair is already solved from the start.
    The cube orientation is U:Yellow, D:White, F:Orange, B:Red, L:Green, R:Blue.
    """
    # Define the colors of the face centers for piece identification.
    W, O, B, G, R = (
        cube.get_face("D")[1][1], cube.get_face("F")[1][1],
        cube.get_face("R")[1][1], cube.get_face("L")[1][1],
        cube.get_face("B")[1][1]
    )

    def is_pair_solved(slot_cnr_pos, slot_edg_pos, cnr_colors, edg_colors):
        """Helper function to check if a specific F2L pair is solved."""
        try:
            # A piece is solved if it's in the correct slot with the correct orientation.
            # pycuber defines orientation (0,0,0) as solved.
            cnr = cube.select_pieces(cnr_colors)[0]
            edg = cube.select_pieces(edg_colors)[0]
            return (cnr.pos == slot_cnr_pos and cnr.orientation == (0, 0, 0) and
                    edg.pos == slot_edg_pos and edg.orientation == (0, 0, 0))
        except (IndexError, KeyError):
            return False

    # Check Pair 1: White-Orange-Blue (Front-Right slot)
    if is_pair_solved("DFR", "FR", {W, O, B}, {O, B}):
        return True

    # Check Pair 2: White-Red-Blue (Back-Right slot)
    if is_pair_solved("DBR", "BR", {W, R, B}, {R, B}):
        return True

    # Check Pair 3: White-Red-Green (Back-Left slot)
    if is_pair_solved("DBL", "BL", {W, R, G}, {R, G}):
        return True

    return False

def solve_f2l():
    """
    Finds the shortest sequence to solve a second F2L pair.
    """
    # 1. Set up the initial cube state from the problem description.
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"
    reorientation = "x2 y'"  # Reorients to Yellow top, Orange front.
    
    # All 18 possible moves in the Half-Turn Metric.
    moves = ["U", "U'", "U2", "D", "D'", "D2", "L", "L'", "L2", "R", "R'", "R2", "F", "F'", "F2", "B", "B'", "B2"]
    
    cube = pc.Cube()
    cube(scramble)
    cube(reorientation)
    
    # 2. Set up the Breadth-First Search (BFS).
    initial_formula = cube.formula()
    queue = deque([(initial_formula, [])]) # Queue stores (cube_formula, path_to_state)
    visited = {initial_formula}            # Set to track visited states.

    print("Searching for the shortest sequence to solve a second F2L pair...")
    
    # 3. Run the BFS.
    while queue:
        current_formula, path = queue.popleft()
        
        # Recreate the cube object for analysis.
        current_cube = pc.Cube()
        current_cube(current_formula)
        
        # Check if the current state is a goal state.
        if is_f2l_goal_state(current_cube):
            print("\n--- Solution Found ---")
            equation = " ".join(path)
            length = len(path)
            print(f"The minimum number of moves is: {length}")
            print(f"The optimal sequence is: {equation}")
            return length

        # If not a goal state, generate and enqueue all valid next states.
        for move in moves:
            # Optimization: avoid nonsensical moves like R followed by R'.
            if path:
                last_move_face = path[-1][0]
                if last_move_face == move[0]:
                    continue
            
            next_cube = current_cube.copy()
            next_cube(move)
            next_formula = next_cube.formula()

            if next_formula not in visited:
                visited.add(next_formula)
                new_path = path + [move]
                queue.append((next_formula, new_path))
    
    print("A solution could not be found.")
    return -1

if __name__ == '__main__':
    final_length = solve_f2l()
    # The final answer format is specified by the user.
    # It will be included at the very end of the response.
