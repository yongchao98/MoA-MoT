import collections
import pycuber as pc

def solve_f2l_pairs():
    """
    This script finds the shortest move sequence to solve two F2L pairs on a scrambled cube.

    The process is as follows:
    1. A standard cube is created and the specified scramble is applied.
    2. A goal state is defined based on the new orientation (Yellow top, Orange front).
       The F2L pairs for the White layer are the target.
    3. A Breadth-First Search (BFS) is used to find the shortest path from the
       scrambled state to any state where at least two of these F2L pairs are solved.
    4. The script prints the length and the sequence of the first solution found.
    """
    # Define colors for checking stickers. The pycuber default orientation is:
    # U: White, F: Green, R: Red, D: Yellow, B: Blue, L: Orange
    W = pc.Colour("white")
    G = pc.Colour("green")
    R = pc.Colour("red")
    B = pc.Colour("blue")
    O = pc.Colour("orange")
    Y = pc.Colour("yellow")

    # The scramble given by the user
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"

    # Create a solved cube and apply the scramble
    cube = pc.Cube()
    cube(scramble)

    # --- Goal State Definition ---
    # Johnny's View: U=Yellow, D=White, F=Orange, B=Red, R=Green, L=Blue
    # pycuber's internal representation: U=White, D=Yellow, F=Green, B=Blue, R=Red, L=Orange
    # We need to map Johnny's goal (solved White F2L pairs) to pycuber's fixed coordinates.
    # Johnny's White layer (D) is pycuber's White layer (U).
    #
    # Example mapping for one pair:
    # Johnny's Down-Front-Right (D-F-R) slot holds the White-Orange-Green pair.
    # J_D -> P_U (White Face)
    # J_F -> P_L (Orange Face)
    # J_R -> P_F (Green Face)
    # So, Johnny's DFR slot is the UFL slot in the pycuber model.

    TARGET_PAIRS = {
        # J_DFR (W-O-G) -> P_UFL slot
        "WOG_pair": {"c_pos": "UFL", "c_stickers": {"U": W, "F": G, "L": O},
                     "e_pos": "FL", "e_stickers": {"F": G, "L": O}},
        # J_DFL (W-O-B) -> P_UBL slot
        "WOB_pair": {"c_pos": "UBL", "c_stickers": {"U": W, "B": B, "L": O},
                     "e_pos": "BL", "e_stickers": {"B": B, "L": O}},
        # J_DBL (W-R-B) -> P_UBR slot
        "WRB_pair": {"c_pos": "UBR", "c_stickers": {"U": W, "B": B, "R": R},
                     "e_pos": "BR", "e_stickers": {"B": B, "R": R}},
        # J_DBR (W-R-G) -> P_UFR slot
        "WRG_pair": {"c_pos": "UFR", "c_stickers": {"U": W, "F": G, "R": R},
                     "e_pos": "FR", "e_stickers": {"F": G, "R": R}},
    }

    def count_solved_pairs(c):
        """Counts how many F2L pairs are solved according to the target definition."""
        solved_count = 0
        for pair_data in TARGET_PAIRS.values():
            corner = c[pair_data["c_pos"]]
            edge = c[pair_data["e_pos"]]

            # Check if all stickers of the corner piece match the required colors at that position
            is_corner_solved = all(
                corner[face].colour == colour for face, colour in pair_data["c_stickers"].items()
            )
            # Check if all stickers of the edge piece match
            is_edge_solved = all(
                edge[face].colour == colour for face, colour in pair_data["e_stickers"].items()
            )

            if is_corner_solved and is_edge_solved:
                solved_count += 1
        return solved_count

    # --- Breadth-First Search (BFS) ---
    # Check initial state
    initial_solved = count_solved_pairs(cube)
    if initial_solved >= 2:
        print("A solution of 0 moves was found: 2 or more F2L pairs are already solved.")
        print("Minimum moves: 0")
        print("Solution sequence:")
        return 0

    # Queue for BFS: stores (cube_formula, path_list)
    # A cube's formula is a string representation, perfect for the 'visited' set.
    queue = collections.deque([(cube.formula(), [])])
    visited = {cube.formula()}

    # All 18 possible moves (quarter and half turns)
    all_moves = [pc.Move(m) for m in ("U", "U'", "U2", "D", "D'", "D2", "F", "F'", "F2", "B", "B'", "B2", "R", "R'", "R2", "L", "L'", "L2")]
    
    solution_path = None
    
    while queue:
        current_formula, current_path = queue.popleft()
        
        # This cube state will be used to generate next moves
        parent_cube = pc.Cube()
        parent_cube.from_formula(current_formula)

        if len(current_path) > 8: # Safety break for very long searches
            break

        for move in all_moves:
            # Avoid redundant moves like R R' or U U2 U'
            if current_path:
                last_move_face = current_path[-1][0]
                if move.face == last_move_face:
                    continue
            
            next_cube = parent_cube.copy()
            next_cube(move)
            next_formula = next_cube.formula()

            if next_formula not in visited:
                visited.add(next_formula)
                next_path = current_path + [str(move)]

                if count_solved_pairs(next_cube) >= 2:
                    solution_path = next_path
                    break # Found the shortest solution
                
                queue.append((next_formula, next_path))
        
        if solution_path:
            break

    if solution_path:
        print(f"Minimum moves: {len(solution_path)}")
        # The prompt asks to output each "number" in the equation.
        # This is interpreted as printing each move in the sequence.
        print(f"Solution sequence: {' '.join(solution_path)}")
        return len(solution_path)
    else:
        print("No solution found within the search depth limit.")
        return None

# Run the solver and get the result length for the final answer
result_length = solve_f2l_pairs()
print(f"\n<<<__{result_length}__>>>")
