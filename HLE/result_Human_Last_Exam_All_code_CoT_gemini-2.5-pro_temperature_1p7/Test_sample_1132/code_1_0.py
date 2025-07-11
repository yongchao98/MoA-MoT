import copy

def solve_rubiks_face():
    """
    Calculates the state of the white face of a Rubik's cube after a sequence of moves.
    """

    # Helper functions for rotating a 3x3 matrix representing a face
    def rotate_cw(face):
        return [[face[2 - j][i] for j in range(3)] for i in range(3)]

    def rotate_ccw(face):
        return [[face[j][2 - i] for j in range(3)] for i in range(3)]

    def rotate_180(face):
        return rotate_cw(rotate_cw(face))

    # --- Step 1: Standardize the Initial Cube State ---
    # The initial layouts are given from different viewing angles.
    # We must normalize them to the standard orientation (White=Front, Orange=Up).

    white_data = [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']]
    orange_data = [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']]
    blue_data = [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']]
    yellow_data = [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']]
    green_data = [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']]
    red_data = [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]

    faces = {
        'F': white_data,                   # Standard view: White front, Orange top. OK.
        'U': rotate_180(orange_data),      # Orange face: viewed with Yellow top -> 180 deg rot.
        'R': rotate_ccw(blue_data),        # Blue face: viewed with Orange top -> CCW rot.
        'B': rotate_ccw(yellow_data),      # Yellow face: viewed with Blue top -> CCW rot.
        'L': rotate_180(green_data),       # Green face: viewed with Yellow top -> 180 deg rot.
        'D': rotate_ccw(red_data)          # Red face: viewed with Yellow top -> CCW rot.
    }
    
    # --- Step 2: Execute the 5-step Algorithm ---

    # 1. R (Right face clockwise)
    faces['R'] = rotate_cw(faces['R'])
    temp_col = [row[2] for row in faces['F']]
    for i in range(3): faces['F'][i][2] = faces['D'][i][2]
    for i in range(3): faces['D'][i][2] = faces['B'][2 - i][0]
    for i in range(3): faces['B'][2 - i][0] = faces['U'][i][2]
    for i in range(3): faces['U'][i][2] = temp_col[i]

    # 2. U (Up face clockwise)
    faces['U'] = rotate_cw(faces['U'])
    temp_row = [c for c in faces['F'][0]]
    faces['F'][0] = [c for c in faces['R'][0]]
    faces['R'][0] = [c for c in faces['B'][0]]
    faces['B'][0] = [c for c in faces['L'][0]]
    faces['L'][0] = temp_row
    
    # 3. F (Front face clockwise)
    faces['F'] = rotate_cw(faces['F'])
    temp_U_bottom = [c for c in faces['U'][2]]
    temp_R_left = [row[0] for row in faces['R']]
    temp_D_top = [c for c in faces['D'][0]]
    temp_L_right = [row[2] for row in faces['L']]
    for i in range(3): faces['R'][i][0] = temp_U_bottom[i]
    for i in range(3): faces['D'][0][i] = temp_R_left[2 - i]
    for i in range(3): faces['L'][i][2] = temp_D_top[2-i]
    for i in range(3): faces['U'][2][i] = temp_L_right[2 - i]

    # 4. L' (Left face counter-clockwise)
    faces['L'] = rotate_ccw(faces['L'])
    temp_col = [row[0] for row in faces['F']]
    for i in range(3): faces['F'][i][0] = faces['D'][i][0]
    for i in range(3): faces['D'][i][0] = faces['B'][2 - i][2]
    for i in range(3): faces['B'][2 - i][2] = faces['U'][i][0]
    for i in range(3): faces['U'][i][0] = temp_col[i]

    # 5. D (Down face clockwise)
    faces['D'] = rotate_cw(faces['D'])
    temp_row = [c for c in faces['F'][2]]
    faces['F'][2] = [c for c in faces['R'][2]]
    faces['R'][2] = [c for c in faces['B'][2]]
    faces['B'][2] = [c for c in faces['L'][2]]
    faces['L'][2] = temp_row

    # --- Step 3: Print the final state of the white face ---
    print("The final state of the white face is:")
    final_face = faces['F']
    print("[", end="")
    for i, row in enumerate(final_face):
        if i > 0:
            print(" ", end="")
        print(row, end="")
        if i < len(final_face) - 1:
            print(",")
    print("]")


solve_rubiks_face()