import sys
from collections import deque

try:
    import pycuber as pc
except ImportError:
    print("Error: The 'pycuber' library is required but not found.")
    print("Please install it by running the following command in your shell:")
    print("pip install pycuber")
    sys.exit(1)

def count_solved_f2l_pairs(cube):
    """
    Counts the number of solved F2L pairs in the standard orientation.
    A pair is solved if the correct corner and edge are in their home positions
    with the correct orientation. pycuber's default colors are U=White, F=Green, R=Red.
    """
    solved_count = 0
    # Define the correctly solved pieces for each of the four F2L slots.
    # The '==' comparison on pycuber pieces checks for both identity and orientation.
    
    # F2L Pair 1: Front-Right (White-Green-Red)
    if (cube.get_piece('DFR') == pc.Corner(W='W', G='G', R='R') and
        cube.get_piece('FR') == pc.Edge(G='G', R='R')):
        solved_count += 1
        
    # F2L Pair 2: Front-Left (White-Green-Orange)
    if (cube.get_piece('DFL') == pc.Corner(W='W', G='G', O='O') and
        cube.get_piece('FL') == pc.Edge(G='G', O='O')):
        solved_count += 1

    # F2L Pair 3: Back-Right (White-Blue-Red)
    if (cube.get_piece('DBR') == pc.Corner(W='W', B='B', R='R') and
        cube.get_piece('BR') == pc.Edge(B='B', R='R')):
        solved_count += 1
        
    # F2L Pair 4: Back-Left (White-Blue-Orange)
    if (cube.get_piece('DBL') == pc.Corner(W='W', B='B', O='O') and
        cube.get_piece('BL') == pc.Edge(B='B', O='O')):
        solved_count += 1
        
    return solved_count

def solve_f2l_pairs():
    """
    Main function to set up the cube and find the shortest solution to solve two F2L pairs.
    """
    # Create a solved cube and apply the given scramble.
    my_cube = pc.Cube()
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"
    my_cube(scramble)

    # Define the move mapping for the user's orientation (Yellow top, Orange front).
    # This translates user moves to the moves on our standard model.
    # User U (Yellow) -> Model D (Yellow)
    # User F (Orange) -> Model L (Orange)
    # User R (Green)  -> Model F (Green)
    move_map = {
        "U": "D", "U'": "D'", "U2": "D2",
        "D": "U", "D'": "U'", "D2": "U2",
        "F": "L", "F'": "L'", "F2": "L2",
        "B": "R", "B'": "R'", "B2": "R2",
        "R": "F", "R'": "F'", "R2": "F2",
        "L": "B", "L'": "B'", "L2": "B2",
    }
    user_moves = list(move_map.keys())

    # Use a Breadth-First Search (BFS) to find the shortest move sequence.
    # The queue stores tuples of: (cube_state_in_standard_orientation, path_in_user_orientation)
    queue = deque([(my_cube, [])]) 
    visited = {str(my_cube)}

    # Check if the initial state already satisfies the condition
    if count_solved_f2l_pairs(my_cube) >= 2:
        print("The minimum number of moves to solve two F2L pairs is 0.")
        print("One such sequence is: (already solved)")
        return 0

    while queue:
        current_cube, path = queue.popleft()

        # Explore next possible moves
        for user_move in user_moves:
            # Translate the user's move to the model's move.
            mapped_move = move_map[user_move]
            
            next_cube = current_cube.copy()
            next_cube(mapped_move)
            
            if str(next_cube) not in visited:
                # Check if this new state is the goal.
                if count_solved_f2l_pairs(next_cube) >= 2:
                    new_path = path + [user_move]
                    move_count = len(new_path)
                    move_sequence = " ".join(new_path)
                    
                    print(f"The minimum number of moves to solve two F2L pairs is {move_count}.")
                    print(f"One such sequence is: {move_sequence}")
                    return move_count

                visited.add(str(next_cube))
                new_path = path + [user_move]
                queue.append((next_cube, new_path))
    
    return None # Should not be reached if a solution exists

if __name__ == '__main__':
    final_answer = solve_f2l_pairs()
    if final_answer is not None:
        print(f"<<<{final_answer}>>>")
