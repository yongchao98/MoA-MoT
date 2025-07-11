import pycuber as pc
from collections import deque

def check_f2l_pairs_solved(cube):
    """
    Checks the number of solved F2L pairs for a Yellow-Top, Orange-Front orientation.
    In this orientation: U=Yellow, D=White, F=Orange, B=Red, R=Green, L=Blue.
    """
    
    # Define colors
    W = pc.Colour("white")
    Y = pc.Colour("yellow")
    G = pc.Colour("green")
    B = pc.Colour("blue")
    O = pc.Colour("orange")
    R = pc.Colour("red")
    
    # Define target orientations of face centers for sanity check
    # This must be true for the cube to be in the specified orientation
    if not (cube.U.colour == Y and cube.D.colour == W and
            cube.F.colour == O and cube.B.colour == R and
            cube.R.colour == G and cube.L.colour == B):
        return 0

    solved_count = 0
    
    # Pair definitions: (Slot_Name, Corner_Colors, Edge_Colors)
    pairs_to_check = [
        ("DFR", "FR", {W, O, G}, {O, G}), # Slot 1: Down-Front-Right
        ("DFL", "FL", {W, O, B}, {O, B}), # Slot 2: Down-Front-Left
        ("DBR", "BR", {W, R, G}, {R, G}), # Slot 3: Down-Back-Right
        ("DBL", "BL", {W, R, B}, {R, B}), # Slot 4: Down-Back-Left
    ]

    for c_pos, e_pos, c_cols, e_cols in pairs_to_check:
        try:
            corner = cube.get_piece(c_pos)
            edge = cube.get_piece(e_pos)
            
            # 1. Check if the correct pieces are in the slots
            is_corner_correct = (set(c.colour for c in corner) == c_cols)
            is_edge_correct = (set(c.colour for c in edge) == e_cols)
            
            if is_corner_correct and is_edge_correct:
                # 2. Check if pieces are correctly oriented
                # The Down-facing sticker of the corner must be White.
                # The orientation of the edge is checked by its side stickers.
                # If the D-sticker on corner is white, edge can't be misoriented without also misorienting the corner.
                # A more rigorous check is for edge orientation based on its non-D/U colors.
                d_face, f_face, s_face = c_pos[0], c_pos[1], c_pos[2]

                if corner[d_face].colour == W and edge[f_face].colour == cube.F.colour and edge[s_face].colour == cube.R.colour :
                    # For back slots, the 'front' face is 'B' and 'side' is 'L' or 'R'
                    is_oriented = True
                    if f_face == "B":
                        if edge[f_face].colour != cube.B.colour: is_oriented=False
                    else: # f_face == "F"
                        if edge[f_face].colour != cube.F.colour: is_oriented=False
                    
                    if s_face == "L":
                         if edge[s_face].colour != cube.L.colour: is_oriented=False
                    else: # s_face == "R"
                         if edge[s_face].colour != cube.R.colour: is_oriented=False
                    
                    if is_oriented and corner[d_face].colour == W:
                         solved_count += 1

        except (KeyError, ValueError):
            continue

    return solved_count

def main():
    """ Main function to run the F2L solver """
    scramble_str = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    
    # 1. Setup the initial cube state
    my_cube = pc.Cube()
    my_cube(scramble_str)
    
    # 2. Re-orient the cube (Yellow Top, Orange Front)
    my_cube("x2 z'")
    
    # 3. Setup BFS
    q = deque()
    # state is (cube object, path list)
    q.append((my_cube, []))
    
    # Use get_state() tuple representation for the visited set for efficiency
    visited = {my_cube.get_state()}
    
    # Define all possible moves in Half-Turn Metric
    moves = [
        "U", "U'", "U2", "D", "D'", "D2",
        "R", "R'", "R2", "L", "L'", "L2",
        "F", "F'", "F2", "B", "B'", "B2"
    ]
    
    print("Searching for the shortest sequence to solve two F2L pairs...")

    initial_solved = check_f2l_pairs_solved(my_cube)
    print(f"Initial state has {initial_solved} F2L pair(s) solved.")
    if initial_solved >= 2:
        print("\nTwo F2L pairs are already solved.")
        print("Minimum length: 0")
        print("Move sequence: (None)")
        print("\n<<<0>>>")
        return

    # 4. Run BFS
    while q:
        current_cube, path = q.popleft()
        
        # For every possible move
        for move in moves:
            next_cube = current_cube.copy()
            next_cube(move)
            
            state_key = next_cube.get_state()

            if state_key not in visited:
                
                # Check for goal state
                if check_f2l_pairs_solved(next_cube) >= 2:
                    final_path = path + [move]
                    final_length = len(final_path)
                    final_sequence = " ".join(final_path)
                    
                    print(f"\nSolution found!")
                    print(f"The minimum length is: {final_length}")
                    print(f"The move sequence is: {final_sequence}")
                    print(f"\n<<<{final_length}>>>")
                    return
                
                visited.add(state_key)
                new_path = path + [move]
                q.append((next_cube, new_path))
    
    print("No solution found within a reasonable depth.")

if __name__ == '__main__':
    main()