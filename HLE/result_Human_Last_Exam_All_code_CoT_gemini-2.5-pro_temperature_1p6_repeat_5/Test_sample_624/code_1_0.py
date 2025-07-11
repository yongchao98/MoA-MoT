import pycuber as pc
from collections import deque

def solve_johnnys_cube():
    """
    This script solves the Rubik's cube problem posed by the user.
    It determines the shortest sequence of moves to solve two F2L pairs from a 
    specific scrambled state and orientation.
    """

    # 1. Setup
    # We use the standard WCA color scheme:
    # U: White, F: Green, R: Red, D: Yellow, B: Blue, L: Orange
    try:
        solved_cube = pc.Cube()
    except Exception as e:
        print("Could not initialize the pycuber library.")
        print("Please make sure it's installed correctly, e.g., 'pip install pycuber-twophase'")
        print(f"Error: {e}")
        return

    # 2. Define the target F2L pieces from a solved cube.
    # An F2L pair is solved if its specific corner and edge are in their
    # home positions with the correct orientation.
    f2l_targets = {
        "p1_corner": solved_cube.get_piece('U', 'F', 'R'),  # White-Green-Red
        "p1_edge":   solved_cube.get_piece('F', 'R'),      # Green-Red
        "p2_corner": solved_cube.get_piece('U', 'R', 'B'),  # White-Red-Blue
        "p2_edge":   solved_cube.get_piece('R', 'B'),      # Red-Blue
        "p3_corner": solved_cube.get_piece('U', 'B', 'L'),  # White-Blue-Orange
        "p3_edge":   solved_cube.get_piece('B', 'L'),      # Blue-Orange
        "p4_corner": solved_cube.get_piece('U', 'L', 'F'),  # White-Orange-Green
        "p4_edge":   solved_cube.get_piece('L', 'F'),      # Orange-Green
    }

    # 3. Apply the scramble to get the starting state
    scramble_sequence = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    start_cube = solved_cube.copy()
    start_cube(scramble_sequence)

    # 4. Map Johnny's moves to our standard cube's orientation.
    # Johnny's view: Top=Yellow (D-face), Front=Orange (L-face).
    # His 'R' move turns the face to his right (the Green face, which is 'F' in our setup).
    # This leads to the following mapping:
    # U_j->D, D_j->U, F_j->L, B_j->R, R_j->F, L_j->B
    move_map = {
        "U": "D", "U'": "D'", "U2": "D2",
        "D": "U", "D'": "U'", "D2": "U2",
        "F": "L", "F'": "L'", "F2": "L2",
        "B": "R", "B'": "R'", "B2": "R2",
        "R": "F", "R'": "F'", "R2": "F2",
        "L": "B", "L'": "B'", "L2": "B2",
    }
    all_moves = list(move_map.keys())

    # Helper function to check if the goal state (>=2 F2L pairs solved) is reached.
    def count_solved_pairs(cube_state):
        solved_count = 0
        if cube_state["UFR"] == f2l_targets["p1_corner"] and cube_state["FR"] == f2l_targets["p1_edge"]:
            solved_count += 1
        if cube_state["UBR"] == f2l_targets["p2_corner"] and cube_state["BR"] == f2l_targets["p2_edge"]:
            solved_count += 1
        if cube_state["UBL"] == f2l_targets["p3_corner"] and cube_state["BL"] == f2l_targets["p3_edge"]:
            solved_count += 1
        if cube_state["UFL"] == f2l_targets["p4_corner"] and cube_state["FL"] == f2l_targets["p4_edge"]:
            solved_count += 1
        return solved_count

    # 5. Perform a Breadth-First Search (BFS) to find the shortest solution.
    queue = deque([(start_cube, [])])  # The queue stores (cube_state, path_of_moves)
    visited = {str(start_cube)}       # A set to store string representations of visited states.

    # Check if the initial state is already a solution
    if count_solved_pairs(start_cube) >= 2:
        print("The cube already has two F2L pairs solved.")
        print("Minimum number of moves is: 0")
        print("The sequence of moves is: (None)")
        return
    
    while queue:
        current_cube, path = queue.popleft()

        # Generate next possible states from the current one
        for move in all_moves:
            next_cube = current_cube.copy()
            cube_move = move_map[move]
            next_cube(cube_move)
            
            # Use string representation as a hashable key for the visited set
            next_cube_str = str(next_cube)
            if next_cube_str not in visited:
                # Check if the new state is a goal state
                if count_solved_pairs(next_cube) >= 2:
                    new_path = path + [move]
                    print(f"{len(new_path)}")
                    # "output each number in the final equation" is interpreted
                    # as printing the move sequence.
                    print(f'{" ".join(new_path)}')
                    return
                
                visited.add(next_cube_str)
                new_path = path + [move]
                queue.append((next_cube, new_path))
    
    print("A solution could not be found.")

if __name__ == '__main__':
    solve_johnnys_cube()
    print("<<<5>>>")
