# Note: This script requires the 'pycuber' library.
# You can install it by running: pip install pycuber-next
# The original pycuber is unmaintained, pycuber-next is a maintained fork.

import pycuber as pc
from collections import deque

def find_min_moves_for_two_f2l_pairs():
    """
    This function calculates the minimum number of moves to solve two F2L pairs
    from a specific scrambled state.
    """
    # The scramble sequence given in the problem
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"

    # Reorientation from "White top, Green front" to "Yellow top, Orange front"
    # corresponds to a z2 whole-cube rotation followed by a y' whole-cube rotation.
    reorientation = "z2 y'"

    # Create a solved cube. Default pycuber color scheme:
    # U: White, D: Yellow, F: Green, B: Blue, L: Orange, R: Red
    c = pc.Cube()

    # Apply the scramble and then the reorientation to get the initial state
    c(scramble)
    c(reorientation)

    # After reorientation "z2 y'", the face identities change. We define the colors
    # for our checker functions based on the new orientation.
    # New Up: Yellow (Original D), New Front: Orange (Original L), New Right: Green (Original F), etc.
    WHITE, YELLOW = pc.Square("U"), pc.Square("D")
    GREEN, BLUE = pc.Square("F"), pc.Square("B")
    ORANGE, RED = pc.Square("L"), pc.Square("R")

    # Helper functions to check if each of the four F2L pairs is solved.
    # A pair is solved if its edge and corner pieces are in the correct slot
    # with the correct orientation. This is verified by checking 5 specific stickers.

    def is_fr_solved(cube):
        # Front-Right pair: Orange(F)/Green(R) edge, Orange/Green/White(D) corner
        # `pycuber` colors for this slot: Front=Orange("L"), Right=Green("F"), Down=White("U")
        return (cube.get_face("F")[1][2] == ORANGE and
                cube.get_face("R")[1][0] == GREEN and
                cube.get_face("F")[2][2] == ORANGE and
                cube.get_face("R")[2][0] == GREEN and
                cube.get_face("D")[0][2] == WHITE)

    def is_fl_solved(cube):
        # Front-Left pair: Orange(F)/Blue(L) edge, Orange/Blue/White(D) corner
        # `pycuber` colors: Front=Orange("L"), Left=Blue("B"), Down=White("U")
        return (cube.get_face("F")[1][0] == ORANGE and
                cube.get_face("L")[1][2] == BLUE and
                cube.get_face("F")[2][0] == ORANGE and
                cube.get_face("L")[2][2] == BLUE and
                cube.get_face("D")[0][0] == WHITE)

    def is_br_solved(cube):
        # Back-Right pair: Red(B)/Green(R) edge, Red/Green/White(D) corner
        # `pycuber` colors: Back=Red("R"), Right=Green("F"), Down=White("U")
        return (cube.get_face("B")[1][0] == RED and
                cube.get_face("R")[1][2] == GREEN and
                cube.get_face("B")[2][0] == RED and
                cube.get_face("R")[2][2] == GREEN and
                cube.get_face("D")[2][2] == WHITE)

    def is_bl_solved(cube):
        # Back-Left pair: Red(B)/Blue(L) edge, Red/Blue/White(D) corner
        # `pycuber` colors: Back=Red("R"), Left=Blue("B"), Down=White("U")
        return (cube.get_face("B")[1][2] == RED and
                cube.get_face("L")[1][0] == BLUE and
                cube.get_face("B")[2][2] == RED and
                cube.get_face("L")[2][0] == BLUE and
                cube.get_face("D")[2][0] == WHITE)

    def is_goal(cube):
        # The goal is achieved if at least two F2L pairs are solved.
        solved_count = sum([
            is_fr_solved(cube),
            is_fl_solved(cube),
            is_br_solved(cube),
            is_bl_solved(cube)
        ])
        return solved_count >= 2

    # --- Breadth-First Search (BFS) to find the shortest path ---
    
    # Check if the initial state is already a solution (0 moves).
    if is_goal(c):
        print(0)
        return

    # All possible moves in the Half-Turn Metric (18 moves).
    moves = [
        "U", "U'", "U2", "D", "D'", "D2",
        "L", "L'", "L2", "R", "R'", "R2",
        "F", "F'", "F2", "B", "B'", "B2"
    ]

    # Queue for BFS stores tuples of (cube_state, path_of_moves).
    queue = deque([(c, [])])
    # A set to store string representations of visited cube states to avoid cycles.
    visited = {str(c)}

    while queue:
        current_cube, path = queue.popleft()
        
        # Explore all possible next moves.
        for move in moves:
            next_cube = current_cube.copy()
            next_cube(move)
            
            next_cube_str = str(next_cube)
            
            if next_cube_str not in visited:
                visited.add(next_cube_str)
                new_path = path + [move]
                
                # Check if the new state is a goal state.
                if is_goal(next_cube):
                    # If so, we've found the shortest solution. Print its length.
                    print(len(new_path))
                    return
                
                # Otherwise, add the new state to the queue to explore later.
                queue.append((next_cube, new_path))

if __name__ == '__main__':
    find_min_moves_for_two_f2l_pairs()