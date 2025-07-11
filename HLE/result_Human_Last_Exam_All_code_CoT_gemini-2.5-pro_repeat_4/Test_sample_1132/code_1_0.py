import copy

def solve():
    """
    Simulates a sequence of Rubik's cube moves on a given jumbled state
    and prints the final state of the front face.
    """
    
    # 1. Represent the Cube's initial state
    # F=White, U=Orange, R=Blue, L=Green, B=Yellow, D=Red
    cube = {
        'front': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']], # White face
        'up':    [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']], # Orange face
        'right': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']], # Blue face
        'back':  [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']], # Yellow face
        'left':  [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']], # Green face
        'down':  [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]  # Red face
    }

    # 2. Define Rotation helper functions
    def rotate_face_cw(face):
        """Rotates a 3x3 matrix 90 degrees clockwise."""
        return [[face[2][0], face[1][0], face[0][0]], 
                [face[2][1], face[1][1], face[0][1]], 
                [face[2][2], face[1][2], face[0][2]]]

    def rotate_face_ccw(face):
        """Rotates a 3x3 matrix 90 degrees counter-clockwise."""
        return [[face[0][2], face[1][2], face[2][2]], 
                [face[0][1], face[1][1], face[2][1]], 
                [face[0][0], face[1][0], face[2][0]]]

    # 3. Implement and Execute the 5-step algorithm
    
    # Move 1: R (Right face clockwise)
    # Cycle: Front -> Up -> Back -> Down -> Front
    cube['right'] = rotate_face_cw(cube['right'])
    temp_col = [row[2] for row in cube['front']]
    for i in range(3): cube['front'][i][2] = cube['down'][i][2]
    for i in range(3): cube['down'][i][2]  = cube['back'][2 - i][0]
    for i in range(3): cube['back'][2 - i][0] = cube['up'][i][2]
    for i in range(3): cube['up'][i][2]    = temp_col[i]

    # Move 2: U (Up face clockwise)
    # Cycle: Front -> Left -> Back -> Right -> Front
    cube['up'] = rotate_face_cw(cube['up'])
    temp_row = cube['front'][0][:]
    cube['front'][0] = cube['right'][0]
    cube['right'][0] = cube['back'][0]
    cube['back'][0]  = cube['left'][0]
    cube['left'][0]  = temp_row

    # Move 3: F (Front face clockwise)
    # Cycle: Up -> Right -> Down -> Left -> Up
    cube['front'] = rotate_face_cw(cube['front'])
    temp_row = cube['up'][2][:]
    # Left's right column (reversed) -> Up's bottom row
    cube['up'][2] = [cube['left'][2][2], cube['left'][1][2], cube['left'][0][2]]
    # Down's top row -> Left's right column
    for i in range(3): cube['left'][i][2] = cube['down'][0][i]
    # Right's left column (reversed) -> Down's top row
    cube['down'][0] = [cube['right'][2][0], cube['right'][1][0], cube['right'][0][0]]
    # Original Up's bottom row -> Right's left column
    for i in range(3): cube['right'][i][0] = temp_row[i]

    # Move 4: L' (Left face counter-clockwise)
    # Cycle: Front -> Down -> Back -> Up -> Front
    cube['left'] = rotate_face_ccw(cube['left'])
    temp_col = [row[0] for row in cube['front']]
    for i in range(3): cube['front'][i][0] = cube['up'][i][0]
    for i in range(3): cube['up'][i][0]    = cube['back'][2 - i][2]
    for i in range(3): cube['back'][2 - i][2] = cube['down'][i][0]
    for i in range(3): cube['down'][i][0]  = temp_col[i]
    
    # Move 5: D (Down face clockwise)
    # Cycle: Front -> Right -> Back -> Left -> Front
    cube['down'] = rotate_face_cw(cube['down'])
    temp_row = cube['front'][2][:]
    cube['front'][2] = cube['left'][2]
    cube['left'][2]  = cube['back'][2]
    cube['back'][2]  = cube['right'][2]
    cube['right'][2] = temp_row

    # 4. Output the final state of the white face
    final_face = cube['front']
    print("The final state of the white face is:")
    print(f"[{final_face[0]},")
    print(f" {final_face[1]},")
    print(f" {final_face[2]}]")

solve()