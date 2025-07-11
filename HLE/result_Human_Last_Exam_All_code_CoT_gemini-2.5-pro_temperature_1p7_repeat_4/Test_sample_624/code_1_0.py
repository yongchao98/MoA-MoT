import pycuber as pc
from collections import deque

def solve_f2l_pairs():
    """
    This script finds the shortest move sequence to solve two F2L pairs from a given scramble.
    """

    # Define the 4 F2L pairs by their locations and colors in the target orientation:
    # Yellow top (U), White bottom (D), Orange front (F), Green right (R), Red back (B), Blue left (L)
    F2L_PAIRS_DEFINITION = {
        "Orange-Green": {
            "edge_pos": "FR", "corner_pos": "FRD",
            "edge_colors": {"orange", "green"}, "corner_colors": {"white", "orange", "green"},
            "solved_edge_orientation": {"F": "orange", "R": "green"},
            "solved_corner_orientation": {"D": "white", "F": "orange", "R": "green"}
        },
        "Green-Red": {
            "edge_pos": "BR", "corner_pos": "BRD",
            "edge_colors": {"green", "red"}, "corner_colors": {"white", "green", "red"},
            "solved_edge_orientation": {"B": "red", "R": "green"},
            "solved_corner_orientation": {"D": "white", "B": "red", "R": "green"}
        },
        "Red-Blue": {
            "edge_pos": "BL", "corner_pos": "BLD",
            "edge_colors": {"red", "blue"}, "corner_colors": {"white", "red", "blue"},
            "solved_edge_orientation": {"B": "red", "L": "blue"},
            "solved_corner_orientation": {"D": "white", "B": "red", "L": "blue"}
        },
        "Blue-Orange": {
            "edge_pos": "FL", "corner_pos": "FLD",
            "edge_colors": {"blue", "orange"}, "corner_colors": {"white", "blue", "orange"},
            "solved_edge_orientation": {"F": "orange", "L": "blue"},
            "solved_corner_orientation": {"D": "white", "F": "orange", "L": "blue"}
        }
    }

    def check_pair_solved(cube, pair_info):
        """ Checks if a specific F2L pair is solved """
        try:
            edge = cube.get_cubie(pair_info["edge_pos"])
            corner = cube.get_cubie(pair_info["corner_pos"])

            # Check if the correct pieces are in the slots by verifying their colors
            if not (set(c.color.name for c in edge) == pair_info["edge_colors"] and \
                    set(c.color.name for c in corner) == pair_info["corner_colors"]):
                return False

            # Check orientation of the edge piece
            for face_name, color_name in pair_info["solved_edge_orientation"].items():
                if edge[face_name].color.name != color_name:
                    return False
            
            # Check orientation of the corner piece
            for face_name, color_name in pair_info["solved_corner_orientation"].items():
                if corner[face_name].color.name != color_name:
                    return False
        except KeyError: # If a piece is in the wrong layer, get_cubie can fail
             return False
        
        return True

    def count_solved_f2l_pairs(cube):
        """ Counts how many F2L pairs are solved """
        count = 0
        for pair_name, pair_info in F2L_PAIRS_DEFINITION.items():
            if check_pair_solved(cube, pair_info):
                count += 1
        return count

    # 1. Create a cube and apply the scramble
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    cube = pc.Cube()
    cube(scramble)

    # 2. Re-orient the cube to Yellow top, Orange front
    # This corresponds to a y2 rotation followed by a z rotation.
    cube("y2 z")

    # 3. Perform a Breadth-First Search (BFS) to find the shortest solution
    
    # Check if 2 pairs are solved initially
    if count_solved_f2l_pairs(cube) >= 2:
        print("Two pairs are already solved.")
        print("The minimum number of moves is: 0")
        print("The final equation is: ") # Empty
        return 0

    # Possible moves (Half-Turn Metric)
    moves = [pc.formula.Move(m) for m in "U U' U2 D D' D2 R R' R2 L L' L2 F F' F2 B B' B2".split()]
    
    # Queue stores (cube_state, path_of_moves)
    queue = deque([(cube, pc.Formula())])
    # Visited set stores string representations of cube states to avoid re-visiting
    visited = {cube.tostring()}

    while queue:
        current_cube, path = queue.popleft()

        for move in moves:
            next_cube = current_cube.copy()
            next_cube(move)
            
            state_str = next_cube.tostring()
            if state_str not in visited:
                new_path = path + pc.Formula(str(move))
                
                if count_solved_f2l_pairs(next_cube) >= 2:
                    print(f"The minimum number of moves is: {len(new_path)}")
                    print("The final equation is:")
                    # Print each move in the solution
                    for m in new_path:
                        print(m, end=" ")
                    print() # for a newline
                    return len(new_path)

                visited.add(state_str)
                queue.append((next_cube, new_path))
    return -1 # Should not be reached if a solution exists

if __name__ == "__main__":
    solve_f2l_pairs()