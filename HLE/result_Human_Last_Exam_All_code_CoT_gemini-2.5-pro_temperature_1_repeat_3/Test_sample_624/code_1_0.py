import pycuber as pc
from collections import deque
import itertools

def solve_f2l_pairs():
    """
    Solves the Rubik's Cube problem by finding the shortest sequence of moves
    to solve two F2L pairs from a given scramble and orientation.
    """

    # The scramble provided by the user.
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"

    # The user performs moves from a "Yellow on top, Orange on front" perspective.
    # We need to map these moves to the pycuber's standard "White on top, Green on front" frame.
    # The transformation is equivalent to a whole-cube rotation of 'x2 y''.
    # User's Up (Yellow) -> Cube's Down
    # User's Down (White) -> Cube's Up
    # User's Front (Orange) -> Cube's Left
    # User's Back (Red) -> Cube's Right
    # User's Right (Blue) -> Cube's Back
    # User's Left (Green) -> Cube's Front
    move_map = {
        "U": "D", "U'": "D'", "U2": "D2",
        "D": "U", "D'": "U'", "D2": "U2",
        "F": "L", "F'": "L'", "F2": "L2",
        "B": "R", "B'": "R'", "B2": "R2",
        "R": "B", "R'": "B'", "R2": "B2",
        "L": "F", "L'": "F'", "L2": "F2",
    }
    # All possible moves the user can make.
    user_moves = list(move_map.keys())

    # 1. Set up the initial state of the cube
    cube = pc.Cube()
    cube(scramble)
    
    # 2. Define the goal state (what a solved F2L pair looks like)
    solved_cube = pc.Cube()
    f2l_pairs_definitions = {
        "Green-Orange": (solved_cube.get_cubie("DFL"), solved_cube.get_cubie("FL")),
        "Orange-Blue": (solved_cube.get_cubie("DBL"), solved_cube.get_cubie("BL")),
        "Blue-Red": (solved_cube.get_cubie("DBR"), solved_cube.get_cubie("BR")),
        "Red-Green": (solved_cube.get_cubie("DFR"), solved_cube.get_cubie("FR")),
    }
    
    f2l_pair_slots = {
        "Green-Orange": ("DFL", "FL"),
        "Orange-Blue": ("DBL", "BL"),
        "Blue-Red": ("DBR", "BR"),
        "Red-Green": ("DFR", "FR"),
    }

    def count_solved_f2l_pairs(current_cube):
        """Checks the cube state and returns the number of solved F2L pairs."""
        count = 0
        for pair_name, slots in f2l_pair_slots.items():
            corner_slot, edge_slot = slots
            target_corner, target_edge = f2l_pairs_definitions[pair_name]
            
            corner_is_solved = current_cube.get_cubie(corner_slot) == target_corner
            edge_is_solved = current_cube.get_cubie(edge_slot) == target_edge

            if corner_is_solved and edge_is_solved:
                count += 1
        return count

    # 3. Perform a Breadth-First Search (BFS) to find the shortest path
    # The queue will store tuples of (cube_state, path_of_moves)
    queue = deque([(cube, [])]) 
    # The visited set stores string representations of cube states to avoid cycles
    visited = {str(cube)}

    print("Searching for the shortest sequence to solve two F2L pairs...")

    while queue:
        current_cube, path = queue.popleft()

        # Check if we have reached the goal
        if count_solved_f2l_pairs(current_cube) >= 2:
            move_sequence_str = " ".join(path)
            length = len(path)
            print(f"\nFound a solution of length {length}.")
            print(f"The sequence of moves is: {move_sequence_str}")
            print(f"This solves the F2L pairs and the final answer for the length is {length}.")
            print("The final equation is the sequence of moves itself.")
            final_equation = " ".join(path)
            print(f"Final Equation: {final_equation}")
            return length, final_equation

        # If not, generate next possible states
        for move in user_moves:
            # Create a copy to apply the new move
            next_cube = current_cube.copy()
            
            # Apply the move based on the user's perspective
            cube_move = move_map[move]
            next_cube(cube_move)
            
            # If we haven't seen this cube state before, add it to the queue
            if str(next_cube) not in visited:
                visited.add(str(next_cube))
                new_path = path + [move]
                queue.append((next_cube, new_path))
    
    return None, "No solution found."

if __name__ == '__main__':
    length, sequence = solve_f2l_pairs()
    if length is not None:
        print(f"\nTo reiterate, the exact, minimum length of the smallest sequence is {length}.")
        print(f'<<<__{length}__>>>')
