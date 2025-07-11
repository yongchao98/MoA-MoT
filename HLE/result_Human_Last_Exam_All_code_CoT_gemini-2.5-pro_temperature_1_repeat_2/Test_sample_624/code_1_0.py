import pycuber as pc
from collections import deque

def check_f2l_pairs(cube, target_cube):
    """
    Checks how many F2L pairs are solved on the given cube by comparing it
    to a solved target cube in the same orientation.
    An F2L pair consists of a corner on the D face and its corresponding
    edge in the middle layer, both correctly placed and oriented.
    """
    # Define the piece locations for the four F2L pairs
    f2l_slots = [
        (('F', 'R'), ('D', 'F', 'R')),  # Front-Right
        (('F', 'L'), ('D', 'F', 'L')),  # Front-Left
        (('B', 'R'), ('D', 'B', 'R')),  # Back-Right
        (('B', 'L'), ('D', 'B', 'L')),  # Back-Left
    ]
    
    solved_count = 0
    for edge_pos, corner_pos in f2l_slots:
        # pycuber's piece equality checks for identity, position, and orientation
        current_edge = cube.get_piece(*edge_pos)
        target_edge = target_cube.get_piece(*edge_pos)
        current_corner = cube.get_piece(*corner_pos)
        target_corner = target_cube.get_piece(*corner_pos)
        
        if current_edge == target_edge and current_corner == target_corner:
            solved_count += 1
            
    return solved_count

def solve_task():
    """
    Main function to solve the user's request.
    """
    # 1. Initial scramble sequence
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    
    # 2. Create a cube and apply the scramble
    cube = pc.Cube()
    cube(scramble)
    
    # 3. Apply the whole-cube rotation from (White Top, Green Front)
    # to (Yellow Top, Orange Front). This corresponds to an x2 y' rotation.
    cube("x2 y'")
    
    # 4. Create a solved target cube in the same final orientation for comparison
    target_cube = pc.Cube()
    target_cube("x2 y'")
    
    # 5. Set up the Breadth-First Search (BFS)
    # The queue stores tuples of (cube_state, path_of_moves)
    queue = deque([(cube, [])])
    # The visited set stores string representations of cube states to avoid cycles
    visited = {str(cube)}
    
    # Define the moves allowed for F2L (no D-layer moves to protect the cross)
    moves = [
        "U", "U'", "U2", "R", "R'", "R2", "L", "L'", "L2",
        "F", "F'", "F2", "B", "B'", "B2"
    ]
    
    # Check if the initial state already satisfies the condition
    if check_f2l_pairs(cube, target_cube) >= 2:
        print("Two F2L pairs were already solved after the rotation.")
        print("Final sequence: (No moves needed)")
        print("Minimum moves: 0")
        # For the final answer format
        print("<<<0>>>")
        return

    # 6. Run the BFS
    while queue:
        current_cube, current_path = queue.popleft()
        
        for move in moves:
            # Create a copy to apply the next move
            next_cube = current_cube.copy()
            next_cube(move)
            
            next_cube_str = str(next_cube)
            
            if next_cube_str not in visited:
                visited.add(next_cube_str)
                new_path = current_path + [move]
                
                # Check if the new state is a goal state
                if check_f2l_pairs(next_cube, target_cube) >= 2:
                    solution_str = " ".join(new_path)
                    print(f"Found shortest sequence to solve two F2L pairs.")
                    print(f"Final sequence: {solution_str}")
                    print(f"Minimum moves: {len(new_path)}")
                    # For the final answer format
                    print(f"<<<{len(new_path)}>>>")
                    return

# Execute the solution
solve_task()