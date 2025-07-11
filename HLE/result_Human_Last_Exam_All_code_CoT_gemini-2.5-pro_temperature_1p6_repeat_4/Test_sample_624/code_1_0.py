import pycuber as pc
import collections

def solve_rubiks_f2l_challenge():
    """
    This script solves for the shortest move sequence to complete two F2L pairs
    from a specific scrambled state.
    """

    # Step 1: Define F2L pair checking functions for the target orientation.
    # The orientation is Yellow on top (U), White on bottom (D), Orange on front (F),
    # Red on back (B), Green on right (R), and Blue on left (L).

    def is_fr_solved(c):
        """Checks if the Front-Right (Orange-Green) F2L pair is solved."""
        try:
            corner = c.get_cubie(1, -1, 1) # DFR position
            edge = c.get_cubie(1, 0, 1)   # FR position
            return (isinstance(corner, pc.Corner) and corner["D"].colour == "white" and corner["F"].colour == "orange" and corner["R"].colour == "green" and
                    isinstance(edge, pc.Edge) and edge["F"].colour == "orange" and edge["R"].colour == "green")
        except (KeyError, TypeError):
            return False

    def is_fl_solved(c):
        """Checks if the Front-Left (Orange-Blue) F2L pair is solved."""
        try:
            corner = c.get_cubie(-1, -1, 1) # DFL position
            edge = c.get_cubie(-1, 0, 1)    # FL position
            return (isinstance(corner, pc.Corner) and corner["D"].colour == "white" and corner["F"].colour == "orange" and corner["L"].colour == "blue" and
                    isinstance(edge, pc.Edge) and edge["F"].colour == "orange" and edge["L"].colour == "blue")
        except (KeyError, TypeError):
            return False

    def is_br_solved(c):
        """Checks if the Back-Right (Red-Green) F2L pair is solved."""
        try:
            corner = c.get_cubie(1, -1, -1) # DBR position
            edge = c.get_cubie(1, 0, -1)   # BR position
            return (isinstance(corner, pc.Corner) and corner["D"].colour == "white" and corner["B"].colour == "red" and corner["R"].colour == "green" and
                    isinstance(edge, pc.Edge) and edge["B"].colour == "red" and edge["R"].colour == "green")
        except (KeyError, TypeError):
            return False

    def is_bl_solved(c):
        """Checks if the Back-Left (Red-Blue) F2L pair is solved."""
        try:
            corner = c.get_cubie(-1, -1, -1) # DBL position
            edge = c.get_cubie(-1, 0, -1)    # BL position
            return (isinstance(corner, pc.Corner) and corner["D"].colour == "white" and corner["B"].colour == "red" and corner["L"].colour == "blue" and
                    isinstance(edge, pc.Edge) and edge["B"].colour == "red" and edge["L"].colour == "blue")
        except (KeyError, TypeError):
            return False

    def count_solved_pairs(c):
        """Counts how many F2L pairs are currently solved."""
        return sum([is_fr_solved(c), is_fl_solved(c), is_br_solved(c), is_bl_solved(c)])

    # Step 2: Set up the initial cube state.
    # The scramble is applied with White on top, Green on front.
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"
    cube = pc.Cube()
    cube(scramble)

    # Step 3: Reorient the cube to Yellow on top, Orange on front.
    # A z2 rotation puts Yellow on top, and a y' rotation puts Orange on front.
    cube.perform_algo("z2 y'")

    # Step 4: Use Breadth-First Search (BFS) to find the shortest solution.
    # The Half-Turn Metric (HTM) is used, where any turn of a face (e.g., U, U', U2) counts as 1 move.
    moves = ["U", "U'", "U2", "D", "D'", "D2", "L", "L'", "L2", "R", "R'", "R2", "F", "F'", "F2", "B", "B'", "B2"]
    
    # A queue for BFS, storing tuples of (cube_state, path_taken)
    queue = collections.deque([(cube, [])])
    # A set to keep track of visited cube states to avoid redundant work.
    visited = {str(cube)}
    
    # Check if the initial state already satisfies the condition.
    # The problem implies this isn't the case, but it is good practice to check.
    if count_solved_pairs(cube) >= 2:
        solution = []
    else:
        while queue:
            current_cube, path = queue.popleft()
            last_move_face = path[-1][0] if path else None

            for move in moves:
                # Pruning: Avoid inverse or redundant moves like R R' or R R2.
                if last_move_face and last_move_face == move[0]:
                    continue

                # Create the next state by applying the move.
                next_cube = current_cube.copy()
                next_cube(move)
                state_str = str(next_cube)

                if state_str not in visited:
                    new_path = path + [move]
                    
                    # Check if the new state is a goal state.
                    if count_solved_pairs(next_cube) >= 2:
                        solution = new_path
                        queue.clear() # End the search
                        break
                    
                    visited.add(state_str)
                    queue.append((next_cube, new_path))
    
    # Step 5: Print the final answer.
    print(f"The minimum length to solve two F2L pairs is {len(solution)} moves.")
    # The instruction "output each number in the final equation" is interpreted
    # as printing each move in the final sequence.
    print(f"The shortest sequence of moves is: {' '.join(solution)}")

solve_rubiks_f2l_challenge()
<<<4>>>