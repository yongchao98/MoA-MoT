import collections
import pycuber as pc

def check_f2l_solved(cube):
    """
    Checks how many F2L pairs are solved for the white face (currently D face).
    The cube object is assumed to be in Johnny's orientation:
    U:Yellow, D:White, F:Orange, B:Red, R:Blue, L:Green
    """
    white_color = cube.D.color
    solved_pairs = 0

    # F2L Slot 1: DFR (Corner DFR, Edge FR) -> White-Orange-Blue pair
    try:
        corner = cube["DFR"]
        edge = cube["FR"]
        if (corner.D == white_color and corner.F == cube.F.color and corner.R == cube.R.color and
            edge.F == cube.F.color and edge.R == cube.R.color):
            solved_pairs += 1
    except (KeyError, AttributeError):
        pass

    # F2L Slot 2: DFL (Corner DFL, Edge FL) -> White-Orange-Green pair
    try:
        corner = cube["DFL"]
        edge = cube["FL"]
        if (corner.D == white_color and corner.F == cube.F.color and corner.L == cube.L.color and
            edge.F == cube.F.color and edge.L == cube.L.color):
            solved_pairs += 1
    except (KeyError, AttributeError):
        pass

    # F2L Slot 3: DBL (Corner DBL, Edge BL) -> White-Red-Green pair
    try:
        corner = cube["DBL"]
        edge = cube["BL"]
        if (corner.D == white_color and corner.B == cube.B.color and corner.L == cube.L.color and
            edge.B == cube.B.color and edge.L == cube.L.color):
            solved_pairs += 1
    except (KeyError, AttributeError):
        pass

    # F2L Slot 4: DBR (Corner DBR, Edge BR) -> White-Red-Blue pair
    try:
        corner = cube["DBR"]
        edge = cube["BR"]
        if (corner.D == white_color and corner.B == cube.B.color and corner.R == cube.R.color and
            edge.B == cube.B.color and edge.R == cube.R.color):
            solved_pairs += 1
    except (KeyError, AttributeError):
        pass

    return solved_pairs


def solve_f2l():
    """
    Finds the shortest sequence of moves to solve two F2L pairs.
    """
    # 1. Define the scramble and initial setup
    scramble = pc.Formula("R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D")
    cube = pc.Cube()
    cube(scramble)

    # 2. Re-orient the cube to Johnny's perspective (Yellow top, Orange front)
    # This corresponds to a whole-cube rotation of x2 y'
    cube.transform("x2 y'")

    # 3. Check initial state. If 2 pairs are solved, the length is 0.
    if check_f2l_solved(cube) >= 2:
        print("Minimum moves to solve two F2L pairs: 0")
        print("The equation is the sequence of moves:")
        print("(No moves needed)")
        return

    # 4. Perform a Breadth-First Search (BFS) for the shortest solution
    queue = collections.deque([(cube, [])])
    visited = {str(cube)}
    
    # All 18 possible moves in Half-Turn Metric
    moves = [pc.formula.Move(m) for m in ["U", "U'", "U2", "D", "D'", "D2", "F", "F'", "F2", "B", "B'", "B2", "R", "R'", "R2", "L", "L'", "L2"]]

    while queue:
        current_cube, path = queue.popleft()

        for move in moves:
            next_cube = current_cube.copy()
            next_cube(move)
            
            state_str = str(next_cube)
            if state_str not in visited:
                if check_f2l_solved(next_cube) >= 2:
                    # Solution found
                    solution_path = path + [str(move)]
                    print(f"Minimum moves to solve two F2L pairs: {len(solution_path)}")
                    print("The equation is the sequence of moves:")
                    # Print each move in the final sequence
                    for m in solution_path:
                        print(m, end=' ')
                    print() # for a final newline
                    return

                visited.add(state_str)
                queue.append((next_cube, path + [str(move)]))

if __name__ == '__main__':
    solve_f2l()