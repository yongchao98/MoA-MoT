import sys
from collections import deque
# The 'pycuber' library is required to run this code.
# You can install it with: pip install pycuber
try:
    import pycuber as pc
except ImportError:
    print("Error: The 'pycuber' library is not installed.", file=sys.stderr)
    print("Please install it by running: pip install pycuber", file=sys.stderr)
    sys.exit(1)

def get_piece_by_colours(cube, colours_str):
    """Finds a piece in the cube by its colours by using its fingerprint."""
    sorted_colours = "".join(sorted(list(colours_str)))
    target_piece_signature = None
    if len(sorted_colours) == 2:
        target_piece_signature = pc.Edge(colours=sorted_colours).fp
    elif len(sorted_colours) == 3:
        target_piece_signature = pc.Corner(colours=sorted_colours).fp
    else:
        # This case should not be reached with valid inputs
        raise ValueError("Invalid colour string for a piece.")
    # The cube.children dictionary maps piece fingerprints to piece objects
    return cube.children[target_piece_signature]

def count_solved_f2l_pairs(cube_state):
    """
    Counts the number of solved F2L pairs for a cube oriented with
    Yellow on top (U) and Orange on front (F).
    This corresponds to a z2 rotation from standard W-top, G-front.
    Centers: U=Y, F=O, R=G, D=W, B=R, L=B.
    """
    solved_count = 0
    # Pair definitions (Edge Colours, Corner Colours) mapped to their solved (Position, Position)
    # Colour strings for pycuber are sorted alphabetically (e.g., Green/Orange -> "GO")
    pairs = {
        # Green-Orange pair -> slot FR/DFR
        ("GO", "GOW"): (("F", "R"), ("D", "F", "R")),
        # Orange-Blue pair -> slot FL/DFL
        ("BO", "BOW"): (("F", "L"), ("D", "F", "L")),
        # Blue-Red pair -> slot BL/DBL
        ("BR", "BRW"): (("B", "L"), ("D", "B", "L")),
        # Red-Green pair -> slot BR/DBR
        ("GR", "GRW"): (("B", "R"), ("D", "B", "R")),
    }

    for (edge_colours, corner_colours), (edge_pos, corner_pos) in pairs.items():
        try:
            edge = get_piece_by_colours(cube_state, edge_colours)
            corner = get_piece_by_colours(cube_state, corner_colours)
            
            # A pair is solved if both pieces are in their correct slots with orientation 0.
            if edge.pos == edge_pos and edge.orient == 0 and \
               corner.pos == corner_pos and corner.orient == 0:
                solved_count += 1
        except (KeyError, ValueError):
            print(f"Error: Could not find piece for colours {edge_colours} or {corner_colours}.", file=sys.stderr)
            return -1 # Indicate an error
            
    return solved_count

def solve_f2l_pairs_task():
    """Main function to find the shortest F2L solution."""
    # 1. Initialise a solved cube
    my_cube = pc.Cube()

    # 2. Apply the scramble from the problem description
    scramble = pc.Formula("R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D")
    my_cube(scramble)

    # 3. Apply a z2 rotation for Yellow top, Orange front
    my_cube("z2")

    # 4. Define the goal (>= 2 pairs solved) and check the initial state
    goal_pairs_count = 2
    initial_solved_count = count_solved_f2l_pairs(my_cube)
    
    if initial_solved_count >= goal_pairs_count:
        print("The cube already has 2 or more F2L pairs solved.")
        print("Final Answer: 0 + = 0 moves")
        print("<<<0>>>")
        return

    # 5. Start a Breadth-First Search (BFS)
    queue = deque([(my_cube, pc.Formula())])
    visited = {my_cube.tostr()}

    moves = [pc.Formula(m) for m in ["U", "U'", "U2", "D", "D'", "D2", "F", "F'", "F2", "B", "B'", "B2", "L", "L'", "L2", "R", "R'", "R2"]]

    solution_path = None
    while queue:
        current_cube, path = queue.popleft()

        for move in moves:
            # Optimization: avoid obviously redundant moves like R followed by R'
            if len(path) > 0 and path[-1].surface == move.surface:
                continue

            next_cube = current_cube.copy()
            next_cube(move)
            
            cube_str_rep = next_cube.tostr()
            if cube_str_rep not in visited:
                if count_solved_f2l_pairs(next_cube) >= goal_pairs_count:
                    solution_path = path + move
                    break  # Solution found

                visited.add(cube_str_rep)
                queue.append((next_cube, path + move))
        
        if solution_path:
            break

    # 6. Output the result in the specified format
    if solution_path:
        moves_list = str(solution_path).split()
        length = len(moves_list)
        equation_str = " + ".join(moves_list)
        
        print("The shortest sequence to solve two F2L pairs is:")
        print(f"Final Answer: {equation_str} = {length} moves")
        print(f"<<<{length}>>>")
    else:
        # This case should not be reached for this problem
        print("A solution was not found.", file=sys.stderr)

solve_f2l_pairs_task()