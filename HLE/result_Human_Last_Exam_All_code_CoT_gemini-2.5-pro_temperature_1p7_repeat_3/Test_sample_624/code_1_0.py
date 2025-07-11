import sys
from collections import deque

try:
    import pycuber as pc
except ImportError:
    print("This script requires the 'pycuber' library.")
    print("Please install it using: pip install pycuber")
    sys.exit(1)

# Helper functions to check if a specific F2L pair is solved
# These are for the target orientation: U=Yellow, F=Orange, R=Green, L=Blue, B=Red, D=White

def is_og_solved(cube):
    """Checks if the Orange-Green pair is solved in the FR slot."""
    try:
        return (cube.get_face("F")[1][2].colour == "orange" and
                cube.get_face("R")[1][0].colour == "green" and
                cube.get_face("D")[2][2].colour == "white" and
                cube.get_face("F")[2][2].colour == "orange" and
                cube.get_face("R")[2][0].colour == "green")
    except (IndexError, AttributeError):
        return False

def is_gr_solved(cube):
    """Checks if the Green-Red pair is solved in the BR slot."""
    try:
        return (cube.get_face("R")[1][2].colour == "green" and
                cube.get_face("B")[1][0].colour == "red" and
                cube.get_face("D")[0][2].colour == "white" and
                cube.get_face("R")[2][2].colour == "green" and
                cube.get_face("B")[2][0].colour == "red")
    except (IndexError, AttributeError):
        return False

def is_rb_solved(cube):
    """Checks if the Red-Blue pair is solved in the BL slot."""
    try:
        return (cube.get_face("B")[1][2].colour == "red" and
                cube.get_face("L")[1][0].colour == "blue" and
                cube.get_face("D")[0][0].colour == "white" and
                cube.get_face("B")[2][2].colour == "red" and
                cube.get_face("L")[2][0].colour == "blue")
    except (IndexError, AttributeError):
        return False

def is_bo_solved(cube):
    """Checks if the Blue-Orange pair is solved in the FL slot."""
    try:
        return (cube.get_face("L")[1][2].colour == "blue" and
                cube.get_face("F")[1][0].colour == "orange" and
                cube.get_face("D")[2][0].colour == "white" and
                cube.get_face("L")[2][2].colour == "blue" and
                cube.get_face("F")[2][0].colour == "orange")
    except (IndexError, AttributeError):
        return False

def count_solved_pairs(cube):
    """Counts the number of solved F2L pairs."""
    return sum([
        is_og_solved(cube),
        is_gr_solved(cube),
        is_rb_solved(cube),
        is_bo_solved(cube)
    ])

def solve_f2l_two_pairs():
    """
    Finds the shortest move sequence to solve two F2L pairs from the given scramble.
    """
    # 1. Set up the initial cube state
    c = pc.Cube()
    scramble_str = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"
    c(scramble_str)
    
    # Reorient the cube: Yellow on top, Orange on front (from W-top, G-front)
    c("x2 y'")
    
    # 2. Setup for Breadth-First Search (BFS)
    # The queue stores tuples of (cube_state, path_of_moves)
    queue = deque([(c, [])]) 
    
    # The visited set stores string representations of cubes to avoid redundant checks
    visited = {str(c)}
    
    # Moves in Half-Turn Metric (HTM)
    moves = [pc.Move(m) for m in "U U' U2 R R' R2 F F' F2 D D' D2 L L' L2 B B' B2".split()]

    print("Searching for the shortest move sequence to solve two F2L pairs...")
    print("This may take some time depending on the solution length...")

    # 3. Perform BFS
    while queue:
        current_cube, path = queue.popleft()
        
        # Check if goal is met
        if count_solved_pairs(current_cube) >= 2:
            solution_str = " ".join(str(m) for m in path)
            solution_len = len(path)
            
            print(f"\nSolution found!")
            print(f"The sequence of moves is: {solution_str}")
            for move in path:
                print(move, end=' ')
            print() # for newline
            print(f"The exact minimum length of the smallest sequence is {solution_len} moves.")

            return solution_len

        # To prevent the search from running too long, we can set a depth limit.
        # A typical F2L solution is under 10 moves. Let's cap it reasonably.
        if len(path) > 7:
            continue

        # Generate next states
        last_move_face = path[-1].face if path else None

        for move in moves:
            # Simple pruning: avoid sequences like R R' or R R2 or R R
            if move.face == last_move_face:
                continue

            next_cube = current_cube.copy()
            next_cube(move)
            
            next_cube_str = str(next_cube)
            if next_cube_str not in visited:
                visited.add(next_cube_str)
                new_path = path + [move]
                queue.append((next_cube, new_path))
                
    print("\nNo solution found within the defined search depth.")
    return None

if __name__ == '__main__':
    result = solve_f2l_two_pairs()
    if result is not None:
        print(f"\n<<< {result} >>>")
