import pycuber as pc
import collections

def setup_cube():
    """
    Creates a cube, applies the scramble, and re-orients it as described.
    """
    # Create a standard solved cube (U:White, F:Green, R:Red)
    c = pc.Cube()

    # Apply the scramble sequence
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    c(scramble)

    # Re-orient the cube to have Yellow on Top and Orange on Front.
    # Initial: U=W, D=Y, F=G, L=O. Target: U=Y, F=O.
    # This corresponds to a whole-cube rotation of x2 followed by y.
    c("x2 y")
    return c

def get_solved_f2l_pairs(cube):
    """
    Counts the number of solved F2L pairs on the given cube.
    After re-orientation, the centers are:
    U:Yellow, D:White, F:Orange, R:Green, B:Red, L:Blue
    """
    count = 0
    # The pairs are defined by their colors, which correspond to the centers.
    # F2L Pair 1: DFR slot (White-Orange-Green)
    try:
        corner1 = cube["DFR"]
        edge1 = cube["FR"]
        if (corner1.colors == {"W", "O", "G"} and edge1.colors == {"O", "G"} and
                corner1["D"].color == "white" and corner1["F"].color == "orange" and
                edge1["F"].color == "orange"):
            count += 1
    except (KeyError, ValueError):
        pass # Piece might not be a corner/edge, ignore

    # F2L Pair 2: DRB slot (White-Green-Red)
    try:
        corner2 = cube["DRB"]
        edge2 = cube["RB"]
        if (corner2.colors == {"W", "G", "R"} and edge2.colors == {"G", "R"} and
                corner2["D"].color == "white" and corner2["R"].color == "green" and
                edge2["R"].color == "green"):
            count += 1
    except (KeyError, ValueError):
        pass

    # F2L Pair 3: DBL slot (White-Red-Blue)
    try:
        corner3 = cube["DBL"]
        edge3 = cube["BL"]
        if (corner3.colors == {"W", "R", "B"} and edge3.colors == {"R", "B"} and
                corner3["D"].color == "white" and corner3["B"].color == "red" and
                edge3["B"].color == "red"):
            count += 1
    except (KeyError, ValueError):
        pass

    # F2L Pair 4: DLF slot (White-Blue-Orange)
    try:
        corner4 = cube["DLF"]
        edge4 = cube["LF"]
        if (corner4.colors == {"W", "B", "O"} and edge4.colors == {"B", "O"} and
                corner4["D"].color == "white" and corner4["L"].color == "blue" and
                edge4["L"].color == "blue"):
            count += 1
    except (KeyError, ValueError):
        pass

    return count

def solve():
    """
    Performs a Breadth-First Search to find the shortest move sequence.
    """
    initial_cube = setup_cube()
    
    # Check the initial number of solved pairs
    initial_solved_count = get_solved_f2l_pairs(initial_cube)
    if initial_solved_count >= 2:
        print("The cube already has 2 or more F2L pairs solved.")
        print("Minimum length: 0")
        return 0

    # Initialize queue for BFS: (cube_state, path_of_moves)
    queue = collections.deque([(initial_cube, [])])
    # Keep track of visited states to avoid cycles and redundant work
    visited = {str(initial_cube.flat_str())}

    # All possible moves in Half-Turn Metric
    moves_str = "U U' U2 D D' D2 L L' L2 R R' R2 F F' F2 B B' B2".split()
    moves = [pc.move.Move(s) for s in moves_str]

    while queue:
        current_cube, path = queue.popleft()

        # Generate next possible states
        for move in moves:
            next_cube = current_cube.copy()
            next_cube(move)

            state_str = str(next_cube.flat_str())
            if state_str in visited:
                continue

            new_path = path + [move]
            
            # Check if we have reached the goal state
            if get_solved_f2l_pairs(next_cube) >= 2:
                solution_len = len(new_path)
                solution_str = " ".join(str(m) for m in new_path)
                print(f"The exact, minimum length of the smallest sequence of moves to solve two F2L pairs is: {solution_len}")
                print(f"A shortest sequence is: {solution_str}")
                return solution_len

            visited.add(state_str)
            queue.append((next_cube, new_path))
            
    print("A solution was not found within a reasonable depth.")
    return None

if __name__ == '__main__':
    solve()
