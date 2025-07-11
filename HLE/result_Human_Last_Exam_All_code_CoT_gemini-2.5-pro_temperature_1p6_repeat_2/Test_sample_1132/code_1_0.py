import copy

def rotate_face_cw(face):
    """Rotates a 3x3 face 90 degrees clockwise."""
    return [list(row) for row in zip(*face[::-1])]

def rotate_face_ccw(face):
    """Rotates a 3x3 face 90 degrees counter-clockwise."""
    transposed = [list(row) for row in zip(*face)]
    return transposed[::-1]

def move_R(state):
    """Performs a clockwise R move."""
    new_state = copy.deepcopy(state)
    new_state['R'] = rotate_face_cw(state['R'])
    
    for i in range(3):
        # U -> F
        new_state['F'][i][2] = state['U'][i][2]
        # F -> D
        new_state['D'][i][2] = state['F'][i][2]
        # D -> B (reversed)
        new_state['B'][2-i][0] = state['D'][i][2]
        # B -> U (reversed)
        new_state['U'][i][2] = state['B'][2-i][0]
    return new_state

def move_U(state):
    """Performs a clockwise U move as defined by the problem's example (equivalent to U')."""
    new_state = copy.deepcopy(state)
    # The example implies a U' move, so we rotate CCW
    new_state['U'] = rotate_face_ccw(state['U'])

    # Cycle: F <- R <- B <- L <- F
    temp_F_top = state['F'][0][:]
    new_state['F'][0] = state['R'][0][:]
    new_state['L'][0] = temp_F_top
    new_state['B'][0] = state['L'][0][:]
    new_state['R'][0] = state['B'][0][:]
    return new_state
    
def move_F(state):
    """Performs a clockwise F move."""
    new_state = copy.deepcopy(state)
    new_state['F'] = rotate_face_cw(state['F'])
    
    # temp store old values to avoid conflicts
    temp_U_bot = state['U'][2][:]
    temp_R_left = [row[0] for row in state['R']]
    temp_D_top = state['D'][0][:]
    temp_L_right = [row[2] for row in state['L']]
    
    for i in range(3):
        # U -> R
        new_state['R'][i][0] = temp_U_bot[i]
        # R -> D (reversed)
        new_state['D'][0][i] = temp_R_left[2-i]
        # D -> L (reversed)
        new_state['L'][i][2] = temp_D_top[2-i]
        # L -> U
        new_state['U'][2][i] = temp_L_right[i]
    return new_state

def move_L_prime(state):
    """Performs a counter-clockwise L' move."""
    new_state = copy.deepcopy(state)
    new_state['L'] = rotate_face_ccw(state['L'])

    # Cycle: U -> B -> D -> F -> U
    temp_U_left = [row[0] for row in state['U']]
    temp_B_right = [row[2] for row in state['B']]
    temp_D_left = [row[0] for row in state['D']]
    temp_F_left = [row[0] for row in state['F']]
    
    for i in range(3):
        # B(rev) <- U
        new_state['B'][2-i][2] = temp_U_left[i]
        # D <- B(rev)
        new_state['D'][i][0] = temp_B_right[2-i]
        # F <- D
        new_state['F'][i][0] = temp_D_left[i]
        # U <- F
        new_state['U'][i][0] = temp_F_left[i]
    return new_state
    
def move_D(state):
    """Performs a clockwise D move."""
    new_state = copy.deepcopy(state)
    new_state['D'] = rotate_face_cw(state['D'])
    
    # Cycle: F -> L -> B -> R -> F
    temp_F_bot = state['F'][2][:]
    new_state['L'][2] = temp_F_bot
    new_state['B'][2] = state['L'][2][:]
    new_state['R'][2] = state['B'][2][:]
    new_state['F'][2] = state['R'][2][:]
    return new_state

# Initial state of the jumbled cube
cube = {
    'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White face
    'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange face
    'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue face
    'B': [['Y','W','R'],['B','Y','O'],['G','O','B']], # Yellow face
    'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green face
    'D': [['B','Y','Y'],['R','R','G'],['W','G','O']]  # Red face
}

# Apply the 5-step algorithm
state_1 = move_R(cube)
state_2 = move_U(state_1)
state_3 = move_F(state_2)
state_4 = move_L_prime(state_3)
final_state = move_D(state_4)

# Print the final state of the white (Front) face
final_white_face = final_state['F']
print(str(final_white_face))
