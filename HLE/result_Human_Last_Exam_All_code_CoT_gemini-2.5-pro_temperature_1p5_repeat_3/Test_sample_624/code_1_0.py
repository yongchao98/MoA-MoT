import pycuber as pc
from collections import deque

def solve_f2l_puzzle():
    """
    This function sets up the Rubik's cube state from the problem description
    and performs a Breadth-First Search to find the shortest move sequence 
    to solve a second F2L pair.
    """

    # This helper function counts how many F2L pairs are correctly solved.
    # It assumes the cube is oriented with White on the Down ('D') face.
    # The center colors after re-orientation are:
    # U-Y, D-W, F-O, R-B, B-R, L-G.
    def count_solved_f2l_pairs(cube):
        solved_count = 0
        
        # A pair is "solved" if its corner and edge pieces are in their correct home slots
        # and the corner is correctly oriented (White sticker on the Down face).
        
        # Pair 1: White-Orange-Blue corner & Orange-Blue edge -> Home slot: DFR
        corner1 = cube.get_piece('W', 'O', 'B')
        edge1 = cube.get_piece('O', 'B')
        if corner1.pos == 'DFR' and edge1.pos == 'FR' and corner1.facings['D'] == 'white':
            solved_count += 1
            
        # Pair 2: White-Blue-Red corner & Blue-Red edge -> Home slot: DBR
        corner2 = cube.get_piece('W', 'B', 'R')
        edge2 = cube.get_piece('B', 'R')
        if corner2.pos == 'DBR' and edge2.pos == 'BR' and corner2.facings['D'] == 'white':
            solved_count += 1

        # Pair 3: White-Red-Green corner & Red-Green edge -> Home slot: DBL
        corner3 = cube.get_piece('W', 'R', 'G')
        edge3 = cube.get_piece('R', 'G')
        if corner3.pos == 'DBL' and edge3.pos == 'BL' and corner3.facings['D'] == 'white':
            solved_count += 1
            
        # Pair 4: White-Green-Orange corner & Green-Orange edge -> Home slot: DFL
        corner4 = cube.get_piece('W', 'G', 'O')
        edge4 = cube.get_piece('G', 'O')
        if corner4.pos == 'DFL' and edge4.pos == 'FL' and corner4.facings['D'] == 'white':
            solved_count += 1
            
        return solved_count

    # The scramble given in the problem
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"

    # Step 1: Create a cube and apply the scramble
    my_cube = pc.Cube()
    my_cube(scramble)

    # Step 2: Re-orient the cube to "yellow on top, orange on front"
    # This corresponds to a whole-cube rotation of 'x2 y\''
    my_cube("x2 y'")

    # Step 3: Check initial state. If 2 pairs are solved, length is 0.
    if count_solved_f2l_pairs(my_cube) >= 2:
        print("Two or more F2L pairs are already solved.")
        print("Sequence: (None)")
        print("Length: 0")
        return 0

    # Step 4: Perform a Breadth-First Search (BFS)
    # The queue will store tuples of (cube_state, path_list)
    queue = deque([(my_cube, [])])  
    # The visited set stores string representations of cube states to avoid re-visiting
    visited = {str(my_cube)}       

    # We will check all 18 possible quarter and half turns
    moves = ["U", "U'", "U2", "D", "D'", "D2", 
             "L", "L'", "L2", "R", "R'", "R2", 
             "F", "F'", "F2", "B", "B'", "B2"]

    while queue:
        current_cube, path = queue.popleft()

        # Explore all possible next moves from the current state
        for move in moves:
            next_cube = current_cube.copy()
            next_cube(move)
            
            cube_str = str(next_cube)
            if cube_str not in visited:
                visited.add(cube_str)
                new_path = path + [move]
                
                # Check if this new state meets the goal
                if count_solved_f2l_pairs(next_cube) >= 2:
                    # First solution found is guaranteed to be the shortest
                    final_sequence = " ".join(new_path)
                    print("The minimum length sequence to solve a second F2L pair is:")
                    print(f"Sequence: {final_sequence}")
                    print(f"Length: {len(new_path)}")
                    return len(new_path)

                # If not a solution, add it to the queue to explore from later
                queue.append((next_cube, new_path))
    
    # This part should not be reached if a solution exists
    return None

# Run the solver function
final_length = solve_f2l_puzzle()

if final_length is not None:
    print(f"<<<{final_length}>>>")