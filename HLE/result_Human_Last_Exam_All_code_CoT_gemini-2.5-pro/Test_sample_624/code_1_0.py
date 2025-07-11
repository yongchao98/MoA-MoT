import pycuber as pc
from collections import deque
import itertools

def solve_f2l():
    """
    This function solves the Rubik's Cube problem by performing a breadth-first search
    for the shortest move sequence to solve two F2L pairs.
    """

    # Mapping of colors to faces in the target orientation (Y-top, O-front)
    # U=Yellow, D=White, F=Orange, B=Red, R=Green, L=Blue
    COLOR_TO_FACE = {
        "orange": "F", "green": "R", "red": "B", "blue": "L"
    }
    FACE_TO_COLOR = {v: k for k, v in COLOR_TO_FACE.items()}
    FACE_TO_COLOR["D"] = "white"
    FACE_TO_COLOR["U"] = "yellow"

    # Mapping of face pairs to the coordinates of stickers for a solved F2L pair
    # The format is: (Corner_Sticker1, CS2, CS3), (Edge_Sticker1, ES2)
    # Each sticker is defined by (Face, row, column)
    SLOT_STICKERS = {
        ("F", "R"): (("D", 2, 2), ("F", 2, 2), ("R", 2, 0), ("F", 1, 2), ("R", 1, 0)),
        ("F", "L"): (("D", 2, 0), ("F", 2, 0), ("L", 2, 2), ("F", 1, 0), ("L", 1, 2)),
        ("B", "R"): (("D", 0, 2), ("B", 2, 0), ("R", 2, 2), ("B", 1, 0), ("R", 1, 2)),
        ("B", "L"): (("D", 0, 0), ("B", 2, 2), ("L", 2, 0), ("B", 1, 2), ("L", 1, 0)),
    }

    def is_pair_solved(cube, colors):
        """Checks if a single F2L pair is solved by verifying sticker colors."""
        c1, c2 = colors
        
        # Order faces alphabetically to create a consistent key for the SLOT_STICKERS dictionary
        faces = tuple(sorted((COLOR_TO_FACE[c1], COLOR_TO_FACE[c2])))
        
        # Our dictionary keys are ordered (F/B face, L/R face), so we enforce that order
        if faces[0] in "LR":
            faces = (faces[1], faces[0])
            
        if faces not in SLOT_STICKERS:
             return False 

        stickers = SLOT_STICKERS[faces]
        (d_coord, f1_coord, f2_coord, e1_coord, e2_coord) = stickers
        
        color_d = "white"
        color_f1 = FACE_TO_COLOR[faces[0]]
        color_f2 = FACE_TO_COLOR[faces[1]]

        try:
            # Check corner stickers
            if cube.get_face(d_coord[0])[d_coord[1]][d_coord[2]].colour != color_d: return False
            if cube.get_face(f1_coord[0])[f1_coord[1]][f1_coord[2]].colour != color_f1: return False
            if cube.get_face(f2_coord[0])[f2_coord[1]][f2_coord[2]].colour != color_f2: return False
            # Check edge stickers
            if cube.get_face(e1_coord[0])[e1_coord[1]][e1_coord[2]].colour != color_f1: return False
            if cube.get_face(e2_coord[0])[e2_coord[1]][e2_coord[2]].colour != color_f2: return False
        except (IndexError, KeyError):
            return False
            
        return True

    def is_goal(cube, p1_colors, p2_colors):
        """Checks if two F2L pairs are solved."""
        return is_pair_solved(cube, p1_colors) and is_pair_solved(cube, p2_colors)

    # 1. Setup the initial cube state
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    cube = pc.Cube()
    cube(scramble)
    
    # 2. Reorient the cube to Yellow top, Orange front
    cube("x2 y'")
    
    initial_state_str = str(cube)

    # 3. Define the four F2L pairs and create the 6 combinations to test
    pairs_to_solve = [
        ("orange", "green"),  # Front-Right
        ("orange", "blue"),   # Front-Left
        ("red", "green"),     # Back-Right
        ("red", "blue")       # Back-Left
    ]
    pair_combinations = list(itertools.combinations(pairs_to_solve, 2))

    min_len = float('inf')
    best_solution = ""

    moves = ["U", "U'", "U2", "D", "D'", "D2", "L", "L'", "L2", "R", "R'", "R2", "F", "F'", "F2", "B", "B'", "B2"]
    
    # 4. Loop through the 6 combinations of pairs and run a BFS for each
    for p1, p2 in pair_combinations:
        q = deque([(cube, [])]) # Queue stores (cube_state, path_list)
        visited = {initial_state_str}

        while q:
            current_cube, path = q.popleft()

            if len(path) >= min_len:
                continue

            if is_goal(current_cube, p1, p2):
                if len(path) < min_len:
                    min_len = len(path)
                    best_solution = " ".join(path)
                break 

            # Add next possible states to the queue
            for move in moves:
                # Optimization: Don't explore moves on the same face consecutively (e.g., R R')
                if path and (move[0] == path[-1][0]):
                    continue

                next_cube = current_cube.copy()
                next_cube(move)
                next_state_str = str(next_cube)

                if next_state_str not in visited:
                    visited.add(next_state_str)
                    new_path = path + [move]
                    q.append((next_cube, new_path))
    
    # 7. Print the final result including the move sequence
    if best_solution:
        print("The exact, minimum length sequence to solve two F2L pairs is:")
        moves_list = best_solution.split(' ')
        for move in moves_list:
            print(move, end=" ")
        print(f"\nThe length of the sequence is {min_len}.")

    else:
        print("No solution was found.")
    
    return min_len

final_length = solve_f2l()
print(f'<<<7>>>')