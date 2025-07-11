import sys
from collections import deque

try:
    import pycuber as pc
except ImportError:
    print("Error: This script requires the 'pycuber' library.")
    print("Please install it by running: pip install pycuber-twophase")
    sys.exit(1)

def is_f2l_pair_solved(cube, slot_name):
    """
    Checks if a specific F2L pair is solved for the Yellow Top/Orange Front orientation.
    In this orientation, the center colors are:
    U=Yellow, D=White, F=Orange, B=Red, R=Blue, L=Green.
    """
    try:
        # Check Front-Right slot: Should contain White-Orange-Blue corner & Orange-Blue edge
        if slot_name == "FR":
            corner = cube.get_piece("DFR")
            edge = cube.get_piece("FR")
            return (corner['D'].colour == 'W' and corner['F'].colour == 'O' and corner['R'].colour == 'B' and
                    edge['F'].colour == 'O' and edge['R'].colour == 'B')
        # Check Front-Left slot: Should contain White-Orange-Green corner & Orange-Green edge
        elif slot_name == "FL":
            corner = cube.get_piece("DFL")
            edge = cube.get_piece("FL")
            return (corner['D'].colour == 'W' and corner['F'].colour == 'O' and corner['L'].colour == 'G' and
                    edge['F'].colour == 'O' and edge['L'].colour == 'G')
        # Check Back-Right slot: Should contain White-Red-Blue corner & Red-Blue edge
        elif slot_name == "BR":
            corner = cube.get_piece("DBR")
            edge = cube.get_piece("BR")
            return (corner['D'].colour == 'W' and corner['B'].colour == 'R' and corner['R'].colour == 'B' and
                    edge['B'].colour == 'R' and edge['R'].colour == 'B')
        # Check Back-Left slot: Should contain White-Red-Green corner & Red-Green edge
        elif slot_name == "BL":
            corner = cube.get_piece("DBL")
            edge = cube.get_piece("BL")
            return (corner['D'].colour == 'W' and corner['B'].colour == 'R' and corner['L'].colour == 'G' and
                    edge['B'].colour == 'R' and edge['L'].colour == 'G')
    except (KeyError, IndexError):
        # This occurs if a piece is in the wrong type of slot (e.g., a corner in an edge position).
        return False
    return False

def count_solved_f2l_pairs(cube):
    """Counts how many of the four F2L pairs are solved."""
    slots = ["FR", "FL", "BR", "BL"]
    return sum(1 for slot in slots if is_f2l_pair_solved(cube, slot))

def find_shortest_solution():
    """
    Main function to find the shortest move sequence to solve two F2L pairs.
    """
    # 1. Apply the scramble
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    cube = pc.Cube()
    cube(scramble)

    # 2. Reorient the cube to Yellow top, Orange front
    # Standard: U:W, F:G -> `x2` -> U:Y, F:B -> `y'` -> U:Y, F:O
    cube("x2 y'")

    # 3. Check initial state
    if count_solved_f2l_pairs(cube) >= 2:
        print("Two F2L pairs were already solved after the scramble and reorientation.")
        print("Move sequence: (none)")
        print("Length: 0")
        print("<<<0>>>")
        return

    # 4. Set up for Breadth-First Search (BFS)
    queue = deque([(cube, [])])  # Queue of (cube_state, path_list)
    visited = {str(cube)}
    
    # Standard moves to explore (Half-Turn Metric)
    moves = [
        "U", "U'", "U2", "D", "D'", "D2", "R", "R'", "R2", "L", "L'", "L2",
        "F", "F'", "F2", "B", "B'", "B2"
    ]
    
    # 5. Run BFS
    while queue:
        current_cube, path = queue.popleft()
        
        for move in moves:
            next_cube = current_cube.copy()
            next_cube(move)
            next_cube_state_str = str(next_cube)
            
            if next_cube_state_str not in visited:
                new_path = path + [move]
                visited.add(next_cube_state_str)
                
                # Check if we reached the goal state
                if count_solved_f2l_pairs(next_cube) >= 2:
                    solution_str = " ".join(new_path)
                    solution_len = len(new_path)
                    
                    # 6. Output the result
                    print(f"Minimum length to solve two F2L pairs: {solution_len}")
                    print(f"The shortest move sequence is: {solution_str}")
                    print(f"<<<{solution_len}>>>")
                    return
                
                queue.append((next_cube, new_path))

if __name__ == '__main__':
    find_shortest_solution()