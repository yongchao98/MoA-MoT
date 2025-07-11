import pycuber as pc
from collections import deque

def get_f2l_status(cube):
    """
    Checks the status of the four F2L pairs for a given cube state.
    The cube is assumed to be in the user's orientation:
    Y-top, O-front, W-bottom, R-back, G-left, B-right.
    Returns the number of solved pairs.
    """
    solved_count = 0
    
    # Define the four F2L pairs by their colors and target locations
    pairs_to_check = {
        "FL": ({"white", "orange", "green"}, {"orange", "green"}),
        "FR": ({"white", "orange", "blue"}, {"orange", "blue"}),
        "BL": ({"white", "red", "green"}, {"red", "green"}),
        "BR": ({"white", "red", "blue"}, {"red", "blue"}),
    }

    # Maps pair names to their physical locations and correct colors on each face
    face_map = { "FL": ("D", "F", "L"), "FR": ("D", "F", "R"), "BL": ("D", "B", "L"), "BR": ("D", "B", "R") }
    color_map = {
        "FL": {"D": "white", "F": "orange", "L": "green"}, "FR": {"D": "white", "F": "orange", "R": "blue"},
        "BL": {"D": "white", "B": "red", "L": "green"}, "BR": {"D": "white", "B": "red", "R": "blue"}
    }

    for name, colors in pairs_to_check.items():
        corner_colors, edge_colors = colors
        d_face, f_face, lr_face = face_map[name]

        corner = cube.get_piece(d_face, f_face, lr_face)
        edge = cube.get_piece(f_face, lr_face)

        if set(corner.colors) == corner_colors and set(edge.colors) == edge_colors:
            # Check for correct orientation
            corner_oriented = (corner[d_face] == color_map[name][d_face] and
                               corner[f_face] == color_map[name][f_face] and
                               corner[lr_face] == color_map[name][lr_face])
            edge_oriented = (edge[f_face] == color_map[name][f_face] and
                             edge[lr_face] == color_map[name][lr_face])

            if corner_oriented and edge_oriented:
                solved_count += 1
    
    return solved_count

def find_shortest_f2l_solution(start_cube):
    """
    Performs a Breadth-First Search to find the shortest sequence of moves
    to solve two additional F2L pairs.
    """
    moves = ["U", "U'", "U2", "D", "D'", "D2", "L", "L'", "L2", "R", "R'", "R2", "F", "F'", "F2", "B", "B'", "B2"]
    
    initial_solved_count = get_f2l_status(start_cube)
    target_solved_count = initial_solved_count + 2

    if target_solved_count > 4:
        return "Already 3 or 4 pairs solved, cannot solve 2 more."
    
    if initial_solved_count >= 2:
        return "" # 0 moves needed if 2+ are already solved

    queue = deque([(start_cube.get_state(), [])])
    visited = {start_cube.get_state()}

    while queue:
        current_state_str, current_path = queue.popleft()
        
        for move in moves:
            if current_path: # Avoid redundant moves like R R'
                last_move = current_path[-1]
                if move[0] == last_move[0] and (len(move) + len(last_move) == 2 or move == last_move):
                    continue

            next_cube = pc.Cube(state=current_state_str)
            next_cube(move)
            next_state_str = next_cube.get_state()

            if next_state_str not in visited:
                new_path = current_path + [move]
                solved_count = get_f2l_status(next_cube)
                
                if solved_count >= target_solved_count:
                    return " ".join(new_path)
                
                visited.add(next_state_str)
                queue.append((next_state_str, new_path))
    
    return "No solution found."

if __name__ == "__main__":
    # The scramble sequence provided by Johnny
    scramble_alg = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"

    # Create a solved cube
    my_cube = pc.Cube()

    # Apply the scramble
    my_cube(scramble_alg)

    # Reorient the cube to Yellow on top, Orange on front (from W-top, G-front)
    my_cube("x2 y'")

    # Find the shortest sequence to solve two F2L pairs
    solution = find_shortest_f2l_solution(my_cube)
    
    if solution is not None and solution != "":
        moves = solution.split()
        equation = " + ".join(moves)
        length = len(moves)
        
        print("The final equation to solve two F2L pairs is:")
        print(equation)
        print("\nThe minimum number of moves is:")
        print(length)
    elif solution == "":
        print("Two or more F2L pairs were already solved.")
        print("\nThe minimum number of moves is:")
        print(0)
    else:
        print("A solution could not be found.")

<<<4>>>