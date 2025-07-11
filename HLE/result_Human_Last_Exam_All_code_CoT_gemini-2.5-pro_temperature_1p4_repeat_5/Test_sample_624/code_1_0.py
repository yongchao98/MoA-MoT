import pycuber as pc
from collections import deque
import itertools

def solve_f2l_pairs():
    """
    This script finds the shortest move sequence to solve two F2L pairs
    from the given scrambled state.
    """

    # 1. SETUP THE CUBE

    # Create a cube with WCA standard orientation (White on top, Green on front)
    # pycuber's default is Yellow on top, so we rotate it by x2.
    my_cube = pc.Cube()
    my_cube.perform_move("x2")

    # Apply the user's scramble sequence
    scramble_seq = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    my_cube.perform_algo(scramble_seq)

    # Re-orient the cube to Johnny's perspective: Yellow top, Orange front.
    # This corresponds to a `y x2` rotation from the standard WCA orientation.
    my_cube.perform_move("y x2")

    # 2. DEFINE THE GOAL STATE (SOLVING F2L PAIRS)

    # In this new orientation, the centers define the faces:
    # U=Yellow, D=White, F=Orange, B=Red, R=Blue, L=Green
    F_center = my_cube.F.get_center_colour()
    B_center = my_cube.B.get_center_colour()
    L_center = my_cube.L.get_center_colour()
    R_center = my_cube.R.get_center_colour()
    D_center = my_cube.D.get_center_colour()

    # The four F2L pairs are defined by adjacent centers in the middle layer
    f2l_pairs = [
        (F_center, L_center, "FL"), (F_center, R_center, "FR"),
        (B_center, L_center, "BL"), (B_center, R_center, "BR")
    ]

    def is_pair_solved(cube, pair_info):
        """Checks if a single F2L pair is solved."""
        side1_center, side2_center, name = pair_info
        
        # Determine piece positions from name (e.g., "FL" -> edge at FL, corner at DFL)
        edge_pos = name
        corner_pos = "D" + name

        # Check edge piece identity and orientation
        edge = cube.get_cubie(edge_pos)
        if not {side1_center, side2_center} == {c.colour for c in edge.facings.values()}:
            return False
        if edge[name[0]] != side1_center:
            return False

        # Check corner piece identity and orientation
        corner = cube.get_cubie(corner_pos)
        if not {D_center, side1_center, side2_center} == {c.colour for c in corner.facings.values()}:
            return False
        if corner["D"] != D_center:
            return False
            
        return True

    # 3. FIND THE SHORTEST PATH USING BREADTH-FIRST SEARCH (BFS)
    
    # Check how many pairs are solved initially
    initial_solved_count = sum(1 for p in f2l_pairs if is_pair_solved(my_cube, p))

    if initial_solved_count >= 2:
        print("Two F2L pairs are already solved.")
        print("Move sequence: (no moves needed)")
        print("Length: 0")
        return

    # A queue for the BFS, storing (cube_state, path_of_moves)
    queue = deque([(my_cube, "")])
    # A set to keep track of visited states to avoid cycles and redundant work
    visited = {str(my_cube)}

    print("Searching for the shortest move sequence...")
    while queue:
        current_cube, path = queue.popleft()

        # Try every possible move from the current state
        for move in pc.moves.MOVES:
            next_cube = current_cube.copy()
            next_cube.perform_move(move)
            
            # Use the string representation of the cube for efficient visited checking
            state_str = str(next_cube)
            if state_str not in visited:
                
                # Check if the new state meets the goal
                solved_count = sum(1 for p in f2l_pairs if is_pair_solved(next_cube, p))
                
                if solved_count >= 2:
                    final_path = (path + " " + move).strip()
                    print("Found a solution!")
                    print(f"The shortest sequence to solve two F2L pairs is: {final_path}")
                    print("The sequence is:")
                    for m in final_path.split():
                        print(m)
                    print(f"\nThe minimum length is: {len(final_path.split())}")
                    return

                visited.add(state_str)
                queue.append((next_cube, (path + " " + move).strip()))

if __name__ == '__main__':
    solve_f2l_pairs()
