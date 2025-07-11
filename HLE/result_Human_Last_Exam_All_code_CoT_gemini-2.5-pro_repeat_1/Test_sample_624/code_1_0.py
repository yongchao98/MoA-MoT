import collections
import pycuber as pc

def solve_f2l_problem():
    """
    This script solves the Rubik's cube problem by finding the minimum moves
    to solve two F2L pairs from a given scramble and orientation.
    
    It requires the 'pycuber-ng' library. You can install it using:
    pip install pycuber-ng
    """

    # Define the standard colors for a WCA-oriented cube
    U_COLOR, D_COLOR = pc.Color("white"), pc.Color("yellow")
    F_COLOR, B_COLOR = pc.Color("green"), pc.Color("blue")
    R_COLOR, L_COLOR = pc.Color("red"), pc.Color("orange")

    def get_solved_f2l_count(cube):
        """
        Checks the cube state and returns the number of solved F2L pairs.
        A pair is solved if all its stickers are in the correct positions and have the correct colors.
        This is defined relative to the original orientation (White top, Green front).
        """
        solved_count = 0
        
        # F2L Pair 1: Front-Right (White-Green-Red)
        if (cube['U'][2][2].colour == U_COLOR and cube['F'][0][2].colour == F_COLOR and cube['R'][0][0].colour == R_COLOR and
            cube['F'][1][2].colour == F_COLOR and cube['R'][1][0].colour == R_COLOR):
            solved_count += 1
            
        # F2L Pair 2: Front-Left (White-Green-Orange)
        if (cube['U'][2][0].colour == U_COLOR and cube['F'][0][0].colour == F_COLOR and cube['L'][0][2].colour == L_COLOR and
            cube['F'][1][0].colour == F_COLOR and cube['L'][1][2].colour == L_COLOR):
            solved_count += 1

        # F2L Pair 3: Back-Right (White-Blue-Red)
        if (cube['U'][0][2].colour == U_COLOR and cube['B'][0][0].colour == B_COLOR and cube['R'][0][2].colour == R_COLOR and
            cube['B'][1][0].colour == B_COLOR and cube['R'][1][2].colour == R_COLOR):
            solved_count += 1
            
        # F2L Pair 4: Back-Left (White-Blue-Orange)
        if (cube['U'][0][0].colour == U_COLOR and cube['B'][0][2].colour == B_COLOR and cube['L'][0][0].colour == L_COLOR and
            cube['B'][1][2].colour == B_COLOR and cube['L'][1][0].colour == L_COLOR):
            solved_count += 1
            
        return solved_count

    # The scramble sequence given by the user
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"

    # Create a cube and apply the scramble
    initial_cube = pc.Cube()
    initial_cube(scramble)

    # The user picks up the cube with Yellow on top and Orange on the front.
    # We must map the moves from this new orientation to the cube's original orientation.
    # Solve-Orientation -> Scramble-Orientation
    # U (Yellow top) -> D' 
    # D (White bottom) -> U'
    # F (Orange front) -> L'
    # B (Red back) -> R'
    # R (Green right) -> F'
    # L (Blue left) -> B'
    move_map = {
        "U": "D'", "U'": "D", "U2": "D2",
        "D": "U'", "D'": "U", "D2": "U2",
        "F": "L'", "F'": "L", "F2": "L2",
        "B": "R'", "B'": "R", "B2": "R2",
        "R": "F'", "R'": "F", "R2": "F2",
        "L": "B'", "L'": "B", "L2": "B2",
    }
    all_solve_moves = list(move_map.keys())

    # Check if 2 pairs are already solved before starting the search
    if get_solved_f2l_count(initial_cube) >= 2:
        print("The minimum number of moves to solve two F2L pairs is 0.")
        print("The final equation is: ")
        print("<<<0>>>")
        return

    # Setup for Breadth-First Search (BFS)
    # The queue stores tuples of (cube_state_formula, path_of_moves)
    queue = collections.deque([(initial_cube.formula(), [])])
    # The visited set stores cube formulas to avoid redundant searches
    visited = {initial_cube.formula()}

    while queue:
        current_formula, path = queue.popleft()

        for solve_move in all_solve_moves:
            # Create a new cube state from the current formula
            next_cube = pc.Cube(formula=current_formula)
            
            # Get the equivalent move in the original scramble orientation
            scramble_equiv_move = move_map[solve_move]
            
            # Apply the move
            next_cube(scramble_equiv_move)
            
            new_formula = next_cube.formula()

            if new_formula not in visited:
                new_path = path + [solve_move]
                
                # Check if this new state meets the goal
                if get_solved_f2l_count(next_cube) >= 2:
                    solution_length = len(new_path)
                    solution_path_str = " ".join(new_path)
                    
                    print(f"The minimum number of moves to solve two F2L pairs is {solution_length}.")
                    print(f"The final equation is: {solution_path_str}")
                    print(f"<<<{solution_length}>>>")
                    return
                
                visited.add(new_formula)
                queue.append((new_formula, new_path))

if __name__ == '__main__':
    solve_f2l_problem()