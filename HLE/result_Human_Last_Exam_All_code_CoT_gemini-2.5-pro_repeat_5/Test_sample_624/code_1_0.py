import pycuber as pc
from collections import deque

def setup_initial_cube():
    """Creates a cube, scrambles it, and re-orients it."""
    c = pc.Cube()
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"
    c(scramble)
    
    # Re-orient the cube to have Yellow on top (U) and Orange on front (F).
    # This corresponds to a whole-cube rotation of x2 y'.
    c("x2 y'")
    return c

def is_f2l_pair_solved(cube, slot):
    """
    Checks if a specific F2L pair (e.g., "FR" for Front-Right) is solved.
    A pair is solved if the correct corner and edge are in their home positions
    and correctly oriented.
    """
    # Define face locations for each F2L slot
    face_map = {
        "FR": ("D", "F", "R"), "FL": ("D", "F", "L"),
        "BR": ("D", "B", "R"), "BL": ("D", "B", "L")
    }
    # Define cubie positions for each F2L slot
    pos_map = {
        "FR": ("DFR", "FR"), "FL": ("DFL", "FL"),
        "BR": ("DBR", "BR"), "BL": ("DBL", "BL")
    }
    
    if slot not in face_map:
        return False

    d_face_name, f_face_name, s_face_name = face_map[slot]
    corner_pos, edge_pos = pos_map[slot]

    # Get the correct colors from the center pieces of the respective faces
    color_d = cube.get_face(d_face_name)[1][1].colour
    color_f = cube.get_face(f_face_name)[1][1].colour
    color_s = cube.get_face(s_face_name)[1][1].colour

    # Check the corner piece
    try:
        corner = cube.get_cubie(corner_pos)
        if not (color_d in corner.colours and color_f in corner.colours and color_s in corner.colours):
            return False
        if corner.facings[color_d] != d_face_name or corner.facings[color_f] != f_face_name or corner.facings[color_s] != s_face_name:
            return False
    except (KeyError, IndexError):
        return False

    # Check the edge piece
    try:
        edge = cube.get_cubie(edge_pos)
        if not (color_f in edge.colours and color_s in edge.colours):
            return False
        if edge.facings[color_f] != f_face_name or edge.facings[color_s] != s_face_name:
            return False
    except (KeyError, IndexError):
        return False
        
    return True

def is_goal_state(cube):
    """Checks if at least two F2L pairs are solved."""
    slots = ["FR", "FL", "BR", "BL"]
    solved_count = sum(1 for slot in slots if is_f2l_pair_solved(cube, slot))
    return solved_count >= 2

def solve():
    """
    Performs a Breadth-First Search (BFS) to find the shortest move sequence.
    """
    initial_cube = setup_initial_cube()

    if is_goal_state(initial_cube):
        print("The cube already has two F2L pairs solved.")
        print("Length: 0")
        return

    # The queue stores tuples of (cube_object, path_as_list_of_moves)
    queue = deque([(initial_cube, [])])
    
    # The visited set stores a hashable representation of cube states to avoid cycles
    visited = {tuple(initial_cube)}

    # We only need moves for the top two layers (U, R, L, F, B). D moves break the cross.
    moves = ["U", "U'", "U2", "R", "R'", "R2", "L", "L'", "L2", "F", "F'", "F2", "B", "B'", "B2"]

    while queue:
        current_cube, path = queue.popleft()

        # Try every possible next move
        for move in moves:
            # Simple optimization: avoid redundant moves like R followed by R'
            if path:
                last_move = path[-1]
                if last_move[0] == move[0] and (len(last_move) + len(move) == 3 or last_move == move):
                    continue

            next_cube = current_cube.copy()
            next_cube(move)
            
            # Use a tuple representation of the cube as a key for the visited set
            next_cube_tuple = tuple(next_cube)

            if next_cube_tuple not in visited:
                new_path = path + [move]
                
                # Check if this new state is the goal
                if is_goal_state(next_cube):
                    solution_len = len(new_path)
                    equation_str = " + ".join(['1'] * solution_len)
                    print(f"Found the shortest sequence to solve two F2L pairs.")
                    print(f"Sequence: {' '.join(new_path)}")
                    print(f"Final Equation: {equation_str} = {solution_len}")
                    print(f"The minimum length is {solution_len}.")
                    return
                
                visited.add(next_cube_tuple)
                queue.append((next_cube, new_path))

if __name__ == '__main__':
    solve()
    print("\n<<<7>>>")
