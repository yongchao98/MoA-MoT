import pycuber as pc
from collections import deque

def solve_f2l():
    """
    Finds the shortest move sequence to solve two F2L pairs from a given scramble and orientation.
    """

    # The scramble sequence provided by the user.
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"

    # Create a new, solved cube.
    cube = pc.Cube()

    # Apply the scramble.
    cube(scramble)

    # Re-orient the cube: x2 y' rotation to get Yellow on top, Orange on front.
    # pycuber's rotate method performs whole-cube rotations.
    # rotate_x(2) puts Yellow top/Blue front.
    # rotate_y(-1) then puts Orange front.
    cube.rotate_x(2)
    cube.rotate_y(-1)

    # After the rotation, the new face color mapping is:
    # U: Yellow, D: White, F: Orange, B: Red, L: Green, R: Blue
    
    # Define the solved state for each of the 4 F2L pairs in the new orientation.
    # Each dictionary specifies the pieces' colors, positions, and sticker orientations.
    f2l_slots = [
        {
            "name": "Front-Right (Orange-Blue)",
            "corner_pos": ("D", "F", "R"), "edge_pos": ("F", "R"),
            "corner_colors": {"white", "orange", "blue"},
            "edge_colors": {"orange", "blue"},
            "corner_map": {"D": "white", "F": "orange", "R": "blue"},
            "edge_map": {"F": "orange", "R": "blue"}
        },
        {
            "name": "Front-Left (Orange-Green)",
            "corner_pos": ("D", "F", "L"), "edge_pos": ("F", "L"),
            "corner_colors": {"white", "orange", "green"},
            "edge_colors": {"orange", "green"},
            "corner_map": {"D": "white", "F": "orange", "L": "green"},
            "edge_map": {"F": "orange", "L": "green"}
        },
        {
            "name": "Back-Right (Red-Blue)",
            "corner_pos": ("D", "B", "R"), "edge_pos": ("B", "R"),
            "corner_colors": {"white", "red", "blue"},
            "edge_colors": {"red", "blue"},
            "corner_map": {"D": "white", "B": "red", "R": "blue"},
            "edge_map": {"B": "red", "R": "blue"}
        },
        {
            "name": "Back-Left (Red-Green)",
            "corner_pos": ("D", "B", "L"), "edge_pos": ("B", "L"),
            "corner_colors": {"white", "red", "green"},
            "edge_colors": {"red", "green"},
            "corner_map": {"D": "white", "B": "red", "L": "green"},
            "edge_map": {"B": "red", "L": "green"}
        }
    ]

    def count_solved_f2l_pairs(c):
        """Checks the cube and returns the number of solved F2L pairs."""
        solved_count = 0
        for slot in f2l_slots:
            try:
                corner = c.get_cubie(*slot["corner_pos"])
                edge = c.get_cubie(*slot["edge_pos"])

                # Check if the pieces at the target slots have the correct colors
                if set(c.get_piece(*slot["corner_pos"]).colors) != slot["corner_colors"]:
                    continue
                if set(c.get_piece(*slot["edge_pos"]).colors) != slot["edge_colors"]:
                    continue

                # If colors match, check if orientation is also correct
                is_corner_oriented = all(corner[face].color == color for face, color in slot["corner_map"].items())
                is_edge_oriented = all(edge[face].color == color for face, color in slot["edge_map"].items())

                if is_corner_oriented and is_edge_oriented:
                    solved_count += 1
            except KeyError:
                # This can happen if a center piece is in the slot somehow, though unlikely
                continue
        return solved_count

    # Define all 18 possible moves in the Half-Turn Metric.
    moves = [pc.Move(m) for m in ["U", "U'", "U2", "D", "D'", "D2", "F", "F'", "F2", "B", "B'", "B2", "L", "L'", "L2", "R", "R'", "R2"]]
    
    # Initialize the Breadth-First Search (BFS)
    # The queue stores tuples of (cube_state, path_of_moves)
    queue = deque([(cube, [])])
    # The visited set stores string representations of cube states to avoid cycles.
    visited = {str(cube)}
    
    # Check if the initial state already solves the problem
    initial_pairs = count_solved_f2l_pairs(cube)
    if initial_pairs >= 2:
        print(f"The cube already has {initial_pairs} F2L pairs solved.")
        print("The minimum number of moves is 0.")
        print("Final answer: 0")
        print("<<<0>>>")
        return

    # Start the search
    print("Searching for the shortest move sequence...")
    while queue:
        current_cube, path = queue.popleft()
        
        # Explore all possible next moves
        for move in moves:
            # Create a copy to apply the new move
            next_cube = current_cube.copy()
            next_cube(move)
            
            # Use a string representation for checking visited states
            next_cube_str = str(next_cube)
            
            if next_cube_str not in visited:
                visited.add(next_cube_str)
                new_path = path + [str(move)]
                
                # Check if this new state is the goal
                solved_pairs = count_solved_f2l_pairs(next_cube)
                if solved_pairs >= 2:
                    # Solution found!
                    solution_length = len(new_path)
                    solution_moves = " ".join(new_path)
                    
                    print(f"\nSolution found to solve {solved_pairs} F2L pairs!")
                    print(f"The exact, minimum length is: {solution_length}")
                    print(f"A possible sequence is: {solution_moves}")
                    
                    # Also output the final answer in the required format
                    print(f"<<<{solution_length}>>>")
                    return

                # If not the goal, add to the queue to explore later
                queue.append((next_cube, new_path))

if __name__ == '__main__':
    solve_f2l()