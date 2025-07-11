import collections
import sys

# The user needs to have the 'pycuber' library installed.
# You can install it with: pip install pycuber
try:
    import pycuber as pc
except ImportError:
    print("Error: This script requires the 'pycuber' library. Please install it using 'pip install pycuber' and try again.", file=sys.stderr)
    sys.exit(1)

def check_f2l_solved_count(cube):
    """
    Checks a pycuber.Cube object for solved F2L pairs in the standard orientation (White on Top).
    An F2L pair is a correctly placed and oriented White corner and its corresponding middle-layer edge.
    """
    # Define the 4 F2L pairs (corner and edge) for the white layer.
    pairs_to_check = [
        # Pair 1: White-Green-Red Corner & Green-Red Edge
        {"c_pos": ("U", "F", "R"), "c_colors": {"U": "white", "F": "green", "R": "red"},
         "e_pos": ("F", "R"), "e_colors": {"F": "green", "R": "red"}},
        # Pair 2: White-Green-Orange Corner & Green-Orange Edge
        {"c_pos": ("U", "F", "L"), "c_colors": {"U": "white", "F": "green", "L": "orange"},
         "e_pos": ("F", "L"), "e_colors": {"F": "green", "L": "orange"}},
        # Pair 3: White-Blue-Red Corner & Blue-Red Edge
        {"c_pos": ("U", "B", "R"), "c_colors": {"U": "white", "B": "blue", "R": "red"},
         "e_pos": ("B", "R"), "e_colors": {"B": "blue", "R": "red"}},
        # Pair 4: White-Blue-Orange Corner & Blue-Orange Edge
        {"c_pos": ("U", "L", "B"), "c_colors": {"U": "white", "L": "orange", "B": "blue"},
         "e_pos": ("L", "B"), "e_colors": {"L": "orange", "B": "blue"}},
    ]

    solved_count = 0
    for p in pairs_to_check:
        try:
            # Get the cubie pieces at the target locations
            corner = cube.get_piece(*p["c_pos"])
            edge = cube.get_piece(*p["e_pos"])

            # Check if the corner is the correct piece and oriented correctly
            corner_ok = all(corner[face].color == color for face, color in p["c_colors"].items())
            
            # Check if the edge is the correct piece and oriented correctly
            edge_ok = all(edge[face].color == color for face, color in p["e_colors"].items())

            if corner_ok and edge_ok:
                solved_count += 1
        except (KeyError, AttributeError):
            # This block is entered if the piece at the target location is of the wrong type (e.g., an edge in a corner spot),
            # which means the pair is not solved.
            continue
            
    return solved_count

def get_cube_identifier(cube):
    """
    Creates a unique, hashable string representation for a cube state.
    This is used to keep track of visited states in the search.
    """
    return str(cube.get_all_pieces())


def find_f2l_solution():
    """
    Initializes the cube, runs the BFS search, and prints the solution.
    """
    # 1. Create a solved cube and apply the scramble.
    # Initial orientation: White on top, Green on front (pycuber default)
    my_cube = pc.Cube()
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"
    my_cube(scramble)

    # 2. Define the mapping for moves due to reorientation.
    # User View: Yellow top, Orange front. Cube Model View: White top, Green front.
    # This reorientation corresponds to a physical `x2 y'` rotation.
    move_map = {
        "U": "D", "D": "U", "F": "L", "B": "R", "L": "B", "R": "F",
    }
    
    # 3. Set up the Breadth-First Search (BFS).
    queue = collections.deque([(my_cube, [])]) 
    visited = {get_cube_identifier(my_cube)}
    
    # All possible moves in Half-Turn Metric
    user_moves = [face + mod for face in "UDFBLR" for mod in ["", "'", "2"]]

    # Check if initial state is already a solution
    if check_f2l_solved_count(my_cube) >= 2:
        print("0")
        return

    # 4. Run the BFS search.
    while queue:
        current_cube, path = queue.popleft()

        if len(path) > 7: # Pruning the search for performance. F2L solutions are typically short.
             continue

        # Try every possible next move
        for user_move in user_moves:
            next_cube = current_cube.copy()
            
            face = user_move[0]
            mod = user_move[1:]
            actual_move = move_map[face] + mod
            next_cube(actual_move)
            
            identifier = get_cube_identifier(next_cube)
            if identifier not in visited:
                visited.add(identifier)
                new_path = path + [user_move]
                
                if check_f2l_solved_count(next_cube) >= 2:
                    # Solution found!
                    solution_str = " + ".join(new_path)
                    solution_len = len(new_path)
                    print(f"{solution_str} = {solution_len}")
                    return
                
                queue.append((next_cube, new_path))
    
    print("No solution found within the searched depth.")

if __name__ == '__main__':
    find_f2l_solution()