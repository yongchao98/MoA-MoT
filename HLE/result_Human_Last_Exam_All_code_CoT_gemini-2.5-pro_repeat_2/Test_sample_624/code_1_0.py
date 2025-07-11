import collections
import sys

# The pycuber library is required to run this script.
# If you don't have it, you can install it by running:
# pip install pycuber

try:
    import pycuber as pc
except ImportError:
    print("Error: The 'pycuber' library is not found.", file=sys.stderr)
    print("Please install it by running: pip install pycuber", file=sys.stderr)
    sys.exit(1)

def solve_johnnys_cube():
    """
    This script models Johnny's Rubik's cube problem and finds the shortest
    move sequence to solve two F2L pairs.
    """

    # The scramble sequence given in the problem
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"

    # Johnny's orientation: yellow on top, orange on the front.
    # In pycuber's standard model: U=white, D=yellow, F=green, B=blue, R=red, L=orange
    # To get from standard to Johnny's, we perform a whole-cube rotation:
    # x2 (puts yellow on top) then y' (rotates to bring orange to front).
    # The resulting face orientation is:
    # Top(U)=yellow, Bottom(D)=white, Front(F)=orange, Back(B)=red, Right(R)=blue, Left(L)=green
    orientation_rotation = "x2 y'"

    # --- Goal State Checking Functions ---
    # These functions check if an F2L pair is solved by verifying sticker colours.

    # Pair 1: White-Orange-Green (Slot: DFL/FL)
    def is_wog_solved(c):
        return (c.colour_at('DFL', 'D') == 'white' and
                c.colour_at('DFL', 'F') == 'orange' and
                c.colour_at('DFL', 'L') == 'green' and
                c.colour_at('FL', 'F') == 'orange' and
                c.colour_at('FL', 'L') == 'green')

    # Pair 2: White-Orange-Blue (Slot: DFR/FR)
    def is_wob_solved(c):
        return (c.colour_at('DFR', 'D') == 'white' and
                c.colour_at('DFR', 'F') == 'orange' and
                c.colour_at('DFR', 'R') == 'blue' and
                c.colour_at('FR', 'F') == 'orange' and
                c.colour_at('FR', 'R') == 'blue')

    # Pair 3: White-Red-Blue (Slot: DBR/BR)
    def is_wrb_solved(c):
        return (c.colour_at('DBR', 'D') == 'white' and
                c.colour_at('DBR', 'B') == 'red' and
                c.colour_at('DBR', 'R') == 'blue' and
                c.colour_at('BR', 'B') == 'red' and
                c.colour_at('BR', 'R') == 'blue')

    # Pair 4: White-Red-Green (Slot: DBL/BL)
    def is_wrg_solved(c):
        return (c.colour_at('DBL', 'D') == 'white' and
                c.colour_at('DBL', 'B') == 'red' and
                c.colour_at('DBL', 'L') == 'green' and
                c.colour_at('BL', 'B') == 'red' and
                c.colour_at('BL', 'L') == 'green')

    # The goal is reached if at least two pairs are solved.
    def goal_reached(c):
        check_functions = [is_wog_solved, is_wob_solved, is_wrb_solved, is_wrg_solved]
        solved_count = sum(f(c) for f in check_functions)
        return solved_count >= 2

    # --- Main Logic ---
    # 1. Set up the initial cube state
    my_cube = pc.Cube()
    my_cube(scramble)
    my_cube(orientation_rotation)

    # 2. Check if the goal is already met (0 moves)
    if goal_reached(my_cube):
        print("Two F2L pairs are already solved.")
        print("Minimum length of the move sequence: 0")
        print("The final equation is:")
        print(0)
        return

    # 3. Perform a Breadth-First Search (BFS) for the shortest path
    moves = [pc.Move(m) for m in ["U", "U'", "U2", "D", "D'", "D2", "F", "F'", "F2", "B", "B'", "B2", "R", "R'", "R2", "L", "L'", "L2"]]
    
    queue = collections.deque([(my_cube, [])]) 
    visited = {str(my_cube)}

    while queue:
        current_cube, path = queue.popleft()

        for move in moves:
            # Optimization: avoid moves like R R' or R R2
            if path and move.name[0] == path[-1].name[0]:
                continue

            next_cube = current_cube.copy()
            next_cube(move)
            
            next_cube_str = str(next_cube)
            if next_cube_str not in visited:
                new_path = path + [move]
                
                if goal_reached(next_cube):
                    solution_path = " ".join([str(m) for m in new_path])
                    solution_length = len(new_path)
                    
                    print(f"The minimum number of moves to solve two F2L pairs is {solution_length}.")
                    print(f"An example solution sequence is: {solution_path}")
                    print("The final equation is:")
                    print(solution_length)
                    return

                visited.add(next_cube_str)
                queue.append((next_cube, new_path))

if __name__ == '__main__':
    solve_johnnys_cube()