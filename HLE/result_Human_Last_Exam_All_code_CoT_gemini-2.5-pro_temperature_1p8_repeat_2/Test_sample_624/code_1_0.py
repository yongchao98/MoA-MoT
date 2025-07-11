import pycuber as pc
from collections import deque

def get_cube_signature(cube):
    """Gets a unique, hashable string representation of the cube state."""
    return str(cube)

def is_f2l_pair_solved(cube, pair_info):
    """Checks if a single F2L pair is solved by verifying piece identity and orientation."""
    try:
        # Check corner piece identity by colors
        corner = cube.get_cubie(**pair_info["corner_pos"])
        if set(c.colour for c in corner) != pair_info["corner_colors"]:
            return False
        # Check corner orientation by sticker colors
        for face, color in pair_info["corner_stickers"].items():
            if corner[face].colour != color:
                return False

        # Check edge piece identity by colors
        edge = cube.get_cubie(**pair_info["edge_pos"])
        if set(c.colour for c in edge) != pair_info["edge_colors"]:
            return False
        # Check edge orientation by sticker colors
        for face, color in pair_info["edge_stickers"].items():
            if edge[face].colour != color:
                return False
    except (KeyError, ValueError):
        # A piece might not be in a standard cubie slot during moves
        return False

    return True

def count_solved_f2l_pairs(cube, all_pairs_info):
    """Counts how many F2L pairs are solved on the given cube."""
    count = 0
    for pair_name, pair_info in all_pairs_info.items():
        if is_f2l_pair_solved(cube, pair_info):
            count += 1
    return count

def solve_f2l():
    """
    Main function to set up the cube and find the solution for the first two F2L pairs.
    """
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"

    # Create a solved cube
    cube = pc.Cube()
    # Apply the scramble
    cube(scramble)

    # Reorient the cube to Johnny's view: Yellow top (x2), Orange front (y')
    # This results in: U=Y, D=W, F=O, B=R, R=G, L=B
    reorientation = pc.Formula("x2 y'")
    cube = reorientation.perform(cube)

    # Define F2L pairs for the current orientation (U=Y, D=W, F=O, R=G, L=B)
    f2l_pairs_definition = {
        "DFR": { # Down-Front-Right (White-Orange-Green)
            "corner_pos": pc.Cubie(x=1, y=1, z=-1), "edge_pos": pc.Cubie(x=1, y=1, z=0),
            "corner_colors": {"white", "orange", "green"}, "edge_colors": {"orange", "green"},
            "corner_stickers": {"D": "white", "F": "orange", "R": "green"}, "edge_stickers": {"F": "orange", "R": "green"},
        },
        "DFL": { # Down-Front-Left (White-Orange-Blue)
            "corner_pos": pc.Cubie(x=-1, y=1, z=-1), "edge_pos": pc.Cubie(x=-1, y=1, z=0),
            "corner_colors": {"white", "orange", "blue"}, "edge_colors": {"orange", "blue"},
            "corner_stickers": {"D": "white", "F": "orange", "L": "blue"}, "edge_stickers": {"F": "orange", "L": "blue"},
        },
        "DBL": { # Down-Back-Left (White-Red-Blue)
            "corner_pos": pc.Cubie(x=-1, y=-1, z=-1), "edge_pos": pc.Cubie(x=-1, y=-1, z=0),
            "corner_colors": {"white", "red", "blue"}, "edge_colors": {"red", "blue"},
            "corner_stickers": {"D": "white", "B": "red", "L": "blue"}, "edge_stickers": {"B": "red", "L": "blue"},
        },
        "DBR": { # Down-Back-Right (White-Red-Green)
            "corner_pos": pc.Cubie(x=1, y=-1, z=-1), "edge_pos": pc.Cubie(x=1, y=-1, z=0),
            "corner_colors": {"white", "red", "green"}, "edge_colors": {"red", "green"},
            "corner_stickers": {"D": "white", "B": "red", "R": "green"}, "edge_stickers": {"B": "red", "R": "green"},
        }
    }

    # BFS to find the shortest path to solving a second pair
    initial_pairs_count = count_solved_f2l_pairs(cube, f2l_pairs_definition)
    
    if initial_pairs_count >= 2:
        print(0)
        print("Two F2L pairs are already solved.")
        return

    q = deque([(cube, "")])  # Queue of (cube_state, path_string)
    visited = {get_cube_signature(cube)}
    
    moves = ["U", "U'", "U2", "D", "D'", "D2", "F", "F'", "F2",
             "B", "B'", "B2", "L", "L'", "L2", "R", "R'", "R2"]

    while q:
        current_cube, current_path = q.popleft()
        
        for move in moves:
            next_cube = current_cube.copy()
            next_cube(move)
            
            sig = get_cube_signature(next_cube)
            if sig not in visited:
                
                if count_solved_f2l_pairs(next_cube, f2l_pairs_definition) >= 2:
                    solution_path = (current_path + " " + move).strip()
                    solution_moves = solution_path.split(' ')
                    print(f"The minimum number of moves to solve two F2L pairs is: {len(solution_moves)}")
                    print(f"The sequence of moves is: {solution_path}")
                    return

                visited.add(sig)
                new_path = (current_path + " " + move).strip()
                q.append((next_cube, new_path))

if __name__ == '__main__':
    solve_f2l()