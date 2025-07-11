import collections
import pycuber as pc

def count_solved_f2l_pairs(cube):
    """
    Checks a cube state and returns the number of solved F2L pairs.
    Assumes Yellow-Up, Orange-Front orientation.
    Centers: U=Y, D=W, F=O, B=R(Red), R=G(Green), L=B(Blue)
    """
    solved_pairs = 0

    # F2L stickers locations, indexed by [face][row][col]
    # Pair 1: Front-Right slot (Orange/Green)
    # Corner: DFR, Edge: FR
    if (cube["D"][2][2].colour == "W" and cube["F"][2][2].colour == "O" and
        cube["R"][2][0].colour == "G" and cube["F"][1][2].colour == "O" and
        cube["R"][1][0].colour == "G"):
        solved_pairs += 1

    # Pair 2: Front-Left slot (Orange/Blue)
    # Corner: DFL, Edge: FL
    if (cube["D"][2][0].colour == "W" and cube["F"][2][0].colour == "O" and
        cube["L"][2][2].colour == "B" and cube["F"][1][0].colour == "O" and
        cube["L"][1][2].colour == "B"):
        solved_pairs += 1

    # Pair 3: Back-Right slot (Red/Green)
    # Corner: DBR, Edge: BR
    if (cube["D"][0][2].colour == "W" and cube["B"][2][0].colour == "R" and
        cube["R"][2][2].colour == "G" and cube["B"][1][0].colour == "R" and
        cube["R"][1][2].colour == "G"):
        solved_pairs += 1

    # Pair 4: Back-Left slot (Red/Blue)
    # Corner: DBL, Edge: BL
    if (cube["D"][0][0].colour == "W" and cube["B"][2][2].colour == "R" and
        cube["L"][2][0].colour == "B" and cube["B"][1][2].colour == "R" and
        cube["L"][1][0].colour == "B"):
        solved_pairs += 1
        
    return solved_pairs

def solve_f2l():
    """
    Finds the shortest move sequence to solve two F2L pairs.
    """
    # Johnny's scramble sequence
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    
    # Create a solved cube
    my_cube = pc.Cube()

    # Apply the scramble
    my_cube(scramble)

    # Reorient the cube: Yellow on top (U), Orange on front (F)
    # This corresponds to a z2 y' rotation from the standard W-up G-front
    my_cube.perform_alg("z2 y'")

    # Set of moves to explore (Half-Turn Metric)
    moves = ["R", "R'", "R2", "L", "L'", "L2", "U", "U'", "U2", 
             "D", "D'", "D2", "F", "F'", "F2", "B", "B'", "B2"]

    # Queue for BFS: stores tuples of (cube_object, path_string)
    queue = collections.deque([(my_cube, "")])
    # Visited set to store hashable cube states to avoid cycles
    visited = {my_cube.get_state()}

    print("Searching for the shortest sequence to solve two F2L pairs...")

    while queue:
        current_cube, path = queue.popleft()

        # Check if the current state is the goal
        if count_solved_f2l_pairs(current_cube) >= 2:
            path = path.strip()
            length = len(path.split())
            print("\n--- Solution Found! ---")
            print(f"Minimum length to solve two F2L pairs: {length}")
            print(f"The final equation (move sequence) is: {path}")
            return length

        # Generate next states
        for move in moves:
            next_cube = current_cube.copy()
            next_cube(move)
            
            # Check if state has been visited
            if next_cube.get_state() not in visited:
                visited.add(next_cube.get_state())
                new_path = path + move + " "
                queue.append((next_cube, new_path))
    
    print("No solution found (this should not happen for a valid cube).")
    return -1


if __name__ == "__main__":
    result_length = solve_f2l()
    if result_length != -1:
        print(f"\n<<< {result_length} >>>")
