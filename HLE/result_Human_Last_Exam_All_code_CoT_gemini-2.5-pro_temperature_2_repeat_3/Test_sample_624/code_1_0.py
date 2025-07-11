import sys
from collections import deque

def solve_f2l_challenge():
    """
    This script solves the Rubik's cube challenge by simulating the scramble,
    reorienting the cube, and performing a Breadth-First Search to find the
    shortest move sequence to solve two F2L pairs.
    """
    try:
        import pycuber as pc
    except ImportError:
        print("Error: The 'pycuber' library is required to run this script.", file=sys.stderr)
        print("Please install it using: pip install pycuber", file=sys.stderr)
        sys.exit(1)

    # The scramble is specified with adjacent moves, which need a space.
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"

    # Create a solved cube and apply the scramble
    cube = pc.Cube()
    cube(scramble)

    # Reorient the cube from White-Top/Green-Front to Yellow-Top/Orange-Front
    # This corresponds to a whole-cube rotation of x2 y
    cube("x2 y")

    # --- F2L State Checkers ---
    # In the new orientation: U=Y, D=W, F=O, B=R, R=G, L=B
    W, F_c, R_c, L_c, B_c = "W", "O", "G", "B", "R"

    def is_fr_solved(c):
        # FR Slot (Front-Right): Orange/Green. Corner: W/O/G. Edge: O/G
        return (c.get_face("D")[2][2].color == W and
                c.get_face("F")[2][2].color == F_c and
                c.get_face("R")[2][2].color == R_c and
                c.get_face("F")[1][2].color == F_c and
                c.get_face("R")[1][2].color == R_c)

    def is_fl_solved(c):
        # FL Slot (Front-Left): Orange/Blue. Corner: W/O/B. Edge: O/B
        return (c.get_face("D")[2][0].color == W and
                c.get_face("F")[2][0].color == F_c and
                c.get_face("L")[2][2].color == L_c and
                c.get_face("F")[1][0].color == F_c and
                c.get_face("L")[1][2].color == L_c)

    def is_br_solved(c):
        # BR Slot (Back-Right): Red/Green. Corner: W/R/G. Edge: R/G
        return (c.get_face("D")[0][2].color == W and
                c.get_face("B")[2][0].color == B_c and
                c.get_face("R")[2][0].color == R_c and
                c.get_face("B")[1][0].color == B_c and
                c.get_face("R")[1][0].color == R_c)

    def is_bl_solved(c):
        # BL Slot (Back-Left): Red/Blue. Corner: W/R/B. Edge: R/B
        return (c.get_face("D")[0][0].color == W and
                c.get_face("B")[2][2].color == B_c and
                c.get_face("L")[2][0].color == L_c and
                c.get_face("B")[1][2].color == B_c and
                c.get_face("L")[1][0].color == L_c)
        
    def count_solved_f2l_pairs(c):
        return sum([is_fr_solved(c), is_fl_solved(c), is_br_solved(c), is_bl_solved(c)])

    # --- Breadth-First Search ---
    if count_solved_f2l_pairs(cube) >= 2:
        print("0 moves required: The cube already has two F2L pairs solved after reorientation.")
        return

    moves = ["U", "U'", "U2", "D", "D'", "D2", "F", "F'", "F2",
             "B", "B'", "B2", "R", "R'", "R2", "L", "L'", "L2"]
    
    queue = deque([(cube, [])]) # Queue stores (cube_state, path_taken)
    visited = {cube.fp_str()}  # Visited set stores cube fingerprints

    while queue:
        current_cube, path = queue.popleft()
        
        for move_str in moves:
            # Prevent redundant moves like R R' or U U2
            if path and pc.Move(move_str).face == path[-1].face:
                continue
                
            next_cube = current_cube.copy()
            next_cube(move_str)
            
            next_cube_id = next_cube.fp_str()

            if next_cube_id not in visited:
                new_path = path + [pc.Move(move_str)]
                if count_solved_f2l_pairs(next_cube) >= 2:
                    # Solution found
                    move_sequence = " ".join(str(m) for m in new_path)
                    print(f"The minimum length to solve two F2L pairs is {len(new_path)} moves.")
                    print(f"The exact sequence of moves is: {move_sequence}")
                    return

                visited.add(next_cube_id)
                queue.append((next_cube, new_path))
                
    print("Could not find a solution.", file=sys.stderr)


if __name__ == '__main__':
    solve_f2l_challenge()
    print("\n<<<5>>>")