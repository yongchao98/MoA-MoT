import collections
import pycuber as pc

def get_oriented_cube():
    """
    Returns a cube object that has been scrambled and re-oriented
    as per the problem description.
    """
    # Start with a solved cube
    cube = pc.Cube()

    # The scramble sequence
    scramble_str = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"
    scramble_formula = pc.Formula(scramble_str)
    
    # Apply the scramble
    cube(scramble_formula)
    
    # Reorient the cube: Yellow on top (was Bottom, D), Orange on front (was Left, L).
    # This corresponds to an x2 y whole-cube rotation.
    cube("x2 y")
    
    return cube

def check_f2l_pairs(cube):
    """
    Counts the number of solved F2L pairs in the given cube state.
    The cube is assumed to be in the Yellow-top, Orange-front orientation.
    """
    # Define colors
    W, Y = pc.màu.WHITE, pc.màu.YELLOW
    O, R = pc.màu.ORANGE, pc.màu.RED
    G, B = pc.màu.GREEN, pc.màu.BLUE
    
    # Define the colors for each of the 4 F2L pairs and their slot positions
    f2l_slots = {
        "Front-Right": {"corner_colors": {W, O, G}, "edge_colors": {O, G},
                        "corner_pos": ["D", "F", "R"], "edge_pos": ["F", "R"],
                        "edge_face_check": ("F", O)},
        "Front-Left":  {"corner_colors": {W, O, B}, "edge_colors": {O, B},
                        "corner_pos": ["D", "F", "L"], "edge_pos": ["F", "L"],
                        "edge_face_check": ("F", O)},
        "Back-Right":  {"corner_colors": {W, R, G}, "edge_colors": {R, G},
                        "corner_pos": ["D", "B", "R"], "edge_pos": ["B", "R"],
                        "edge_face_check": ("B", R)},
        "Back-Left":   {"corner_colors": {W, R, B}, "edge_colors": {R, B},
                        "corner_pos": ["D", "B", "L"], "edge_pos": ["B", "L"],
                        "edge_face_check": ("B", R)},
    }
    
    solved_count = 0
    for name, slot in f2l_slots.items():
        # Check if correct corner is in the correct slot with correct orientation
        corner = cube.get_piece(*slot["corner_pos"])
        if not (set(corner.colors) == slot["corner_colors"] and corner.facings["D"] == W):
            continue

        # Check if correct edge is in the correct slot with correct orientation
        edge = cube.get_piece(*slot["edge_pos"])
        face, color = slot["edge_face_check"]
        if not (set(edge.colors) == slot["edge_colors"] and edge.facings[face] == color):
            continue
            
        solved_count += 1
        
    return solved_count


def find_shortest_f2l_solution():
    """
    Uses Breadth-First Search (BFS) to find the shortest move sequence
    to solve two F2L pairs.
    """
    start_cube = get_oriented_cube()
    
    # Check if the initial state is already a solution
    if check_f2l_pairs(start_cube) >= 2:
        return ""
        
    # Queue for BFS: stores (cube_state, path_of_moves)
    queue = collections.deque([(start_cube, "")])
    
    # Visited set to avoid redundant computations and cycles
    visited = {start_cube.to_tuple()}
    
    # All 18 possible quarter and half turns
    moves = [face + mod for face in "RUFLDB" for mod in ["", "'", "2"]]

    while queue:
        current_cube, path = queue.popleft()
        
        for move in moves:
            # Create a copy to apply the next move
            next_cube = current_cube.copy()
            next_cube(move)

            # Avoid re-visiting states
            next_cube_tuple = next_cube.to_tuple()
            if next_cube_tuple in visited:
                continue

            # Check if this new state is the goal
            if check_f2l_pairs(next_cube) >= 2:
                # Goal reached, return the sequence of moves
                solution_path = path + " " + move if path else move
                return solution_path

            visited.add(next_cube_tuple)
            new_path = path + " " + move if path else move
            queue.append((next_cube, new_path))
            
    return "Solution not found."

if __name__ == "__main__":
    solution = find_shortest_f2l_solution()
    print("The shortest sequence to solve two F2L pairs is:")
    print(solution)
