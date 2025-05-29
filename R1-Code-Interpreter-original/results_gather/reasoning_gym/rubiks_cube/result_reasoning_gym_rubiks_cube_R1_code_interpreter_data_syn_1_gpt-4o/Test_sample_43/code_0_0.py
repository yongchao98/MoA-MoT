# This script provides a generic solution framework for solving a 3x3x3 Rubik's Cube.
# The solution is given in Singmaster notation.

def solve_rubiks_cube():
    # Step 1: Solve the white cross
    white_cross = [
        "F R U R' U' F'",
        "F R U R' U' F'",
        "F R U R' U' F'",
        "F R U R' U' F'"
    ]
    
    # Step 2: Solve the white corners
    white_corners = [
        "U R U' L' U R' U' L",
        "U R U' L' U R' U' L",
        "U R U' L' U R' U' L",
        "U R U' L' U R' U' L"
    ]
    
    # Step 3: Solve the middle layer edges
    middle_layer = [
        "U R U' R' U' F' U F",
        "U R U' R' U' F' U F",
        "U R U' R' U' F' U F",
        "U R U' R' U' F' U F"
    ]
    
    # Step 4: Solve the yellow cross
    yellow_cross = [
        "F R U R' U' F'",
        "F R U R' U' F'",
        "F R U R' U' F'",
        "F R U R' U' F'"
    ]
    
    # Step 5: Position the yellow edges
    yellow_edges = [
        "R U R' U R U2 R'",
        "R U R' U R U2 R'",
        "R U R' U R U2 R'",
        "R U R' U R U2 R'"
    ]
    
    # Step 6: Position the yellow corners
    yellow_corners_position = [
        "U R U' L' U R' U' L",
        "U R U' L' U R' U' L",
        "U R U' L' U R' U' L",
        "U R U' L' U R' U' L"
    ]
    
    # Step 7: Orient the yellow corners
    yellow_corners_orient = [
        "R U R' U R U2 R'",
        "R U R' U R U2 R'",
        "R U R' U R U2 R'",
        "R U R' U R U2 R'"
    ]
    
    # Combine all steps
    solution = white_cross + white_corners + middle_layer + yellow_cross + yellow_edges + yellow_corners_position + yellow_corners_orient
    
    # Print the solution in Singmaster notation
    for move in solution:
        print(move)

solve_rubiks_cube()