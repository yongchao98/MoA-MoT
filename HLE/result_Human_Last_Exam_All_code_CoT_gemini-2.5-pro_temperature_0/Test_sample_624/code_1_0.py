import pycuber as pc
from collections import deque

# This script requires the 'pycuber' library.
# You can install it by running: pip install pycuber

# 1. Define the scramble sequence and initial setup
scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B'F' U' R L' D R' B R F2 L' F2 L D"

# Create a solved cube and apply the scramble
cube = pc.Cube()
cube(scramble)

# 2. Reorient the cube: White top/Green front -> Yellow top/Orange front
# This corresponds to a z2 rotation (180 degrees around the Front-Back axis)
cube("z2")

# With Yellow on top and Orange on front, the face colors are:
# U: Yellow, D: White, F: Orange, B: Red, R: Blue, L: Green

# 3. Define the properties of the four solved F2L slots for this orientation
f2l_slots = {
    'Front-Right': {
        'corner_pos': 'DFR', 'edge_pos': 'FR',
        'corner_colors': {'White', 'Orange', 'Blue'}, 'edge_colors': {'Orange', 'Blue'},
        'corner_check': ('D', 'White'), 'edge_check': ('F', 'Orange')
    },
    'Front-Left': {
        'corner_pos': 'DFL', 'edge_pos': 'FL',
        'corner_colors': {'White', 'Orange', 'Green'}, 'edge_colors': {'Orange', 'Green'},
        'corner_check': ('D', 'White'), 'edge_check': ('F', 'Orange')
    },
    'Back-Right': {
        'corner_pos': 'DBR', 'edge_pos': 'BR',
        'corner_colors': {'White', 'Red', 'Blue'}, 'edge_colors': {'Red', 'Blue'},
        'corner_check': ('D', 'White'), 'edge_check': ('B', 'Red')
    },
    'Back-Left': {
        'corner_pos': 'DBL', 'edge_pos': 'BL',
        'corner_colors': {'White', 'Red', 'Green'}, 'edge_colors': {'Red', 'Green'},
        'corner_check': ('D', 'White'), 'edge_check': ('B', 'Red')
    }
}

def count_solved_f2l_pairs(c):
    """Counts how many F2L pairs are correctly solved."""
    solved_count = 0
    for slot_info in f2l_slots.values():
        try:
            corner = c.get_cubie(slot_info['corner_pos'])
            edge = c.get_cubie(slot_info['edge_pos'])

            # Check if the correct pieces are in the slot by comparing their colors
            if set(cubie.colour for cubie in corner.facings) != slot_info['corner_colors']:
                continue
            if set(cubie.colour for cubie in edge.facings) != slot_info['edge_colors']:
                continue

            # Check if the pieces are correctly oriented
            corner_face_axis, corner_face_color = slot_info['corner_check']
            edge_face_axis, edge_face_color = slot_info['edge_check']
            if corner[corner_face_axis].colour != corner_face_color:
                continue
            if edge[edge_face_axis].colour != edge_face_color:
                continue
            
            solved_count += 1
        except (KeyError, IndexError):
            continue
    return solved_count

# 4. Perform a Breadth-First Search (BFS) to find the shortest solution
def find_shortest_solution(start_cube):
    # All possible moves in the Half-Turn Metric
    moves = [
        "U", "U'", "U2", "D", "D'", "D2",
        "L", "L'", "L2", "R", "R'", "R2",
        "F", "F'", "F2", "B", "B'", "B2"
    ]

    # Check if the initial state is already a solution
    initial_solved_count = count_solved_f2l_pairs(start_cube)
    if initial_solved_count >= 2:
        print("The cube already has 2 or more F2L pairs solved.")
        print("Solution: (no moves needed)")
        print("Minimum length: 0")
        return 0

    # BFS queue stores tuples of (cube_state, path_of_moves)
    queue = deque([(start_cube, [])])
    # Visited set stores string representations of cube states to avoid cycles
    visited = {str(start_cube)}

    while queue:
        current_cube, path = queue.popleft()

        # Explore all possible next moves
        for move in moves:
            # Avoid redundant sequences like R followed by R'
            if path and pc.Move(move).is_inverse(pc.Move(path[-1])):
                continue

            next_cube = current_cube.copy()
            next_cube(move)
            
            next_cube_str = str(next_cube)
            if next_cube_str not in visited:
                visited.add(next_cube_str)
                new_path = path + [move]
                
                if count_solved_f2l_pairs(next_cube) >= 2:
                    # Solution found, print the results
                    solution_str = " ".join(new_path)
                    print(f"Shortest sequence to solve two F2L pairs: {solution_str}")
                    # As requested, outputting each "number" (move) in the final equation
                    print(f"The moves in the sequence are: {', '.join(new_path)}")
                    print(f"Minimum length: {len(new_path)}")
                    return len(new_path)
                
                queue.append((next_cube, new_path))
    
    print("No solution found.")
    return None

# Run the solver
find_shortest_solution(cube)