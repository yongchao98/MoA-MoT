import pycuber as pc
from collections import deque

def solve_f2l():
    """
    Solves for the shortest move sequence to get two F2L pairs.
    """
    # Define the scramble sequence
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"

    # Standard pycuber Colors
    W, Y, G, B, R, O = pc.Color.W, pc.Color.Y, pc.Color.G, pc.Color.B, pc.Color.R, pc.Color.O

    # --- F2L Pair Checker Functions ---
    # These functions check if a specific F2L pair is solved in the new orientation
    # (Yellow top, Orange front) by checking sticker locations in the standard
    # pycuber orientation (White top, Green front).

    # New FR pair (W-O-G corner, O-G edge) corresponds to old UFL slot.
    # UFL corner stickers: U[2][0], F[0][0], L[0][2]. Required: W, G, O.
    # LF edge stickers: L[1][2], F[1][0]. Required: O, G.
    def check_fr_pair(cube):
        return (cube.U[2][0].color == W and cube.L[0][2].color == O and cube.F[0][0].color == G and
                cube.L[1][2].color == O and cube.F[1][0].color == G)

    # New FL pair (W-O-B corner, O-B edge) corresponds to old ULB slot.
    # ULB corner stickers: U[0][0], L[0][0], B[2][2]. Required: W, O, B.
    # LB edge stickers: L[1][0], B[1][2]. Required: O, B.
    def check_fl_pair(cube):
        return (cube.U[0][0].color == W and cube.L[0][0].color == O and cube.B[2][2].color == B and
                cube.L[1][0].color == O and cube.B[1][2].color == B)

    # New BR pair (W-R-G corner, R-G edge) corresponds to old UFR slot.
    # UFR corner stickers: U[2][2], F[0][2], R[0][0]. Required: W, G, R.
    # RF edge stickers: R[1][0], F[1][2]. Required: R, G.
    def check_br_pair(cube):
        return (cube.U[2][2].color == W and cube.R[0][0].color == R and cube.F[0][2].color == G and
                cube.R[1][0].color == R and cube.F[1][2].color == G)

    # New BL pair (W-R-B corner, R-B edge) corresponds to old UBR slot.
    # UBR corner stickers: U[0][2], B[2][0], R[0][2]. Required: W, B, R.
    # RB edge stickers: R[1][2], B[1][0]. Required: R, B.
    def check_bl_pair(cube):
        return (cube.U[0][2].color == W and cube.B[2][0].color == B and cube.R[0][2].color == R and
                cube.R[1][2].color == R and cube.B[1][0].color == B)

    # List of all checker functions
    F2L_CHECKERS = [check_fr_pair, check_fl_pair, check_br_pair, check_bl_pair]

    def count_solved_f2l_pairs(cube):
        return sum(1 for check in F2L_CHECKERS if check(cube))

    # --- BFS Setup ---
    # Create cube and apply scramble
    my_cube = pc.Cube()
    my_cube(scramble)

    # Check if the cube already has 2+ pairs solved
    initial_solved_count = count_solved_f2l_pairs(my_cube)
    if initial_solved_count >= 2:
        print(f"The cube already has {initial_solved_count} F2L pairs solved.")
        print("The minimum number of moves is 0.")
        print("<<<0>>>")
        return

    # Moves in the new orientation (Y-top, O-front) mapped to old orientation (W-top, G-front)
    MOVE_MAP = {
        "U": "D", "U'": "D'", "U2": "D2",
        "D": "U", "D'": "U'", "D2": "U2",
        "F": "L", "F'": "L'", "F2": "L2",
        "B": "R", "B'": "R'", "B2": "R2",
        "R": "F", "R'": "F'", "R2": "F2",
        "L": "B", "L'": "B'", "L2": "B2",
    }
    ALL_MOVES_NEW = list(MOVE_MAP.keys())
    
    # State key for visited set
    def get_state_key(cube):
        return tuple(s.color for f in "UDLRFB" for r in cube.get_face(f) for s in r)

    # BFS queue stores tuples of (cube_object, path_of_new_moves)
    queue = deque([(my_cube, [])])
    visited = {get_state_key(my_cube)}

    print("Searching for the shortest sequence to solve two F2L pairs...")

    # --- BFS Loop ---
    while queue:
        current_cube, path = queue.popleft()

        # Try all possible next moves
        for move_new in ALL_MOVES_NEW:
            # Pruning: don't turn the same face twice in a row (e.g., R then R')
            if path and path[-1][0] == move_new[0]:
                continue
            
            next_cube = current_cube.copy()
            next_cube(MOVE_MAP[move_new])
            
            state_key = get_state_key(next_cube)
            if state_key in visited:
                continue
            
            visited.add(state_key)
            new_path = path + [move_new]

            if count_solved_f2l_pairs(next_cube) >= 2:
                solution_str = " ".join(new_path)
                print(f"\nSolution Found! Minimum moves: {len(new_path)}")
                # "output each number in the final equation"
                # This prints the sequence like: R2 U' etc.
                print(f"The sequence of moves is: {solution_str}")
                print(f"<<<{len(new_path)}>>>")
                return

            queue.append((next_cube, new_path))
            
    print("No solution found within the search limit.")

# Run the solver
solve_f2l()