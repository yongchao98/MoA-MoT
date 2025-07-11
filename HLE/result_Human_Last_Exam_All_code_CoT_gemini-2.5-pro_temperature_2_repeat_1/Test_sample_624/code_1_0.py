import pycuber as pc
from collections import deque
import sys

def solve_f2l_pairs():
    """
    Finds the shortest move sequence to solve two F2L pairs from a given scramble.
    """
    # Define the scramble and setup the initial cube state
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    initial_cube = pc.Cube()
    initial_cube(scramble)

    # Reorient the cube: "yellow on the top and orange on the front"
    # This corresponds to a x2 y' rotation from the standard W-top, G-front orientation.
    initial_cube("x2 y'")

    # After reorientation, the face centers are:
    # U: Yellow, D: White, F: Orange, R: Green, L: Blue, B: Red
    # We define the 4 F2L pairs for the White (Down) face.
    pairs_info = {
        "FR": {"corner": pc.Corner('W','O','G'), "edge": pc.Edge('O','G')},
        "FL": {"corner": pc.Corner('W','O','B'), "edge": pc.Edge('O','B')},
        "BR": {"corner": pc.Corner('W','R','G'), "edge": pc.Edge('R','G')},
        "BL": {"corner": pc.Corner('W','R','B'), "edge": pc.Edge('R','B')},
    }

    def is_pair_solved(cube, pair_id):
        """ Checks if a specific F2L pair is solved correctly. """
        info = pairs_info[pair_id]
        corner_slot = "D" + pair_id
        edge_slot = pair_id
        
        # Check if the correct pieces are in the correct slots
        if cube.get_piece(corner_slot).colors != info["corner"].colors:
            return False
        if cube.get_piece(edge_slot).colors != info["edge"].colors:
            return False
        
        # Check orientation of the pieces
        if cube[corner_slot].sticker_in("D").color != 'W':
            return False
        main_face_char = edge_slot[0]
        if cube[edge_slot].sticker_in(main_face_char).color != cube.get_face(main_face_char)[1][1].color:
            return False
            
        return True

    def count_solved_pairs(cube):
        """ Counts the total number of solved F2L pairs. """
        return sum(1 for pair_id in pairs_info if is_pair_solved(cube, pair_id))

    # Initialize the Breadth-First Search (BFS)
    # The queue will store tuples of (cube_state, path_of_moves)
    queue = deque([(initial_cube, [])])
    # The visited set stores string representations of cube states to avoid re-visiting
    visited = {str(initial_cube)}

    # Check if the initial state already solves the problem
    if count_solved_pairs(initial_cube) >= 2:
        print(0)
        print("The cube already has two F2L pairs solved.")
        return

    # Define all possible moves in the Half-Turn Metric (HTM)
    all_moves = [
        "U", "U'", "U2", "D", "D'", "D2",
        "L", "L'", "L2", "R", "R'", "R2",
        "F", "F'", "F2", "B", "B'", "B2"
    ]
    
    print("Searching for the shortest sequence... (this might take a moment)", file=sys.stderr)

    # Start the BFS loop
    while queue:
        current_cube, path = queue.popleft()

        # Explore all possible next moves
        for move in all_moves:
            next_cube = current_cube.copy()
            next_cube(move)
            
            next_cube_str = str(next_cube)
            if next_cube_str in visited:
                continue

            # If the new state has 2 or more pairs solved, we found the solution
            if count_solved_pairs(next_cube) >= 2:
                final_path = path + [move]
                print(len(final_path))
                print(" ".join(final_path))
                return

            visited.add(next_cube_str)
            queue.append((next_cube, path + [move]))

solve_f2l_pairs()