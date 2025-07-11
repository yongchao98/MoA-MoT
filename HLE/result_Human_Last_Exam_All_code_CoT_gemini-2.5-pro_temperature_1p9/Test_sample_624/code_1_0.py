import pycuber as pc
from collections import deque

def solve_rubiks_f2l_problem():
    """
    Solves the Rubik's cube F2L problem as described by the user.
    This function encapsulates the entire logic to avoid any global state issues
    and can be run as a self-contained script.
    Note: This process can be computationally intensive and may take a minute to run.
    """

    # Step 1: Create a solved cube and apply the scramble
    try:
        cube = pc.Cube()
    except Exception as e:
        print(f"Failed to initialize pycuber. Please ensure it's installed ('pip install pycuber'). Error: {e}")
        return

    scramble_str = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    
    # Using pc.Formula to handle potentially unspaced moves
    scramble_formula = pc.Formula(scramble_str)
    cube(scramble_formula)

    # Step 2: Re-orient the cube to have Yellow on top (U) and Orange on front (F)
    # This corresponds to a whole-cube rotation of x2 y'
    cube("x2 y'")
    
    start_cube = cube.copy()

    # Step 3: Define what a "solved" F2L pair is in the new orientation.
    # We create a target solved cube in the same orientation to get piece definitions.
    target_cube = pc.Cube()("x2 y'")

    # In the new orientation (Y-U, O-F):
    # F(Front)=Orange, B(Back)=Red, U(Up)=Yellow, D(Down)=White, L(Left)=Green, R(Right)=Blue
    F2L_PAIRS = {
        "Front-Left (Orange-Green)": {
            "edge": target_cube.get_edge('O', 'G'),
            "corner": target_cube.get_corner('W', 'O', 'G')
        },
        "Front-Right (Orange-Blue)": {
            "edge": target_cube.get_edge('O', 'B'),
            "corner": target_cube.get_corner('W', 'O', 'B')
        },
        "Back-Left (Red-Green)": {
            "edge": target_cube.get_edge('R', 'G'),
            "corner": target_cube.get_corner('W', 'R', 'G')
        },
        "Back-Right (Red-Blue)": {
            "edge": target_cube.get_edge('R', 'B'),
            "corner": target_cube.get_corner('W', 'R', 'B')
        }
    }

    def get_solved_pairs(c):
        """Checks a cube state and returns a list of names of solved F2L pairs."""
        solved_list = []
        for name, pieces in F2L_PAIRS.items():
            target_edge = pieces["edge"]
            target_corner = pieces["corner"]
            
            # A piece is identified by its sticker colors.
            # `get_edge_by_sticker` is a hypothetical clearer name for `get_edge`.
            current_edge = c.get_edge(*target_edge.colours)
            current_corner = c.get_corner(*target_corner.colours)

            # Equality check in pycuber verifies both position and orientation.
            if current_edge == target_edge and current_corner == target_corner:
                solved_list.append(name)
        return solved_list

    # Step 4: Find the optimal solution using Breadth-First Search (BFS)
    def find_shortest_f2l_solution(initial_cube):
        # Check initial state
        initial_solved = get_solved_pairs(initial_cube)
        if len(initial_solved) >= 2:
            return "", initial_solved

        # Queue stores tuples of (cube_state, path_of_moves)
        q = deque([(initial_cube, "")])
        # Visited set stores string representations of cube states to avoid cycles
        visited = {str(initial_cube)}

        moves = ["R", "R'", "R2", "L", "L'", "L2", "U", "U'", "U2", 
                 "D", "D'", "D2", "F", "F'", "F2", "B", "B'", "B2"]

        while q:
            current_cube, path = q.popleft()
            
            last_move_face = path.split()[-1][0] if path else None

            for move in moves:
                # Optimization: Don't turn the same face twice in a row (e.g., R then R')
                if last_move_face and move[0] == last_move_face:
                    continue
                
                next_cube = current_cube.copy()
                next_cube(move)
                
                cube_str = str(next_cube)
                if cube_str not in visited:
                    visited.add(cube_str)
                    new_path = f"{path} {move}" if path else move
                    
                    solved_now = get_solved_pairs(next_cube)
                    if len(solved_now) >= 2:
                        # Goal reached: found the shortest sequence
                        return new_path, solved_now
                        
                    q.append((next_cube, new_path))
        
        return "No solution found", []

    # Run the solver and print the result
    solution_path, solved_pairs = find_shortest_f2l_solution(start_cube)

    if "No solution" in solution_path:
        print(solution_path)
    else:
        solution_length = len(solution_path.split())
        print(f"The minimum length sequence to solve two F2L pairs is {solution_length} moves.")
        print(f"The sequence is: {solution_path}")
        print(f"This sequence solves the following two pairs: {', '.join(solved_pairs)}")
        print(f"<<<{solution_length}>>>")


if __name__ == "__main__":
    solve_rubiks_f2l_problem()
