import copy

def rotate_face_clockwise(face):
    """Rotates a 3x3 face 90 degrees clockwise."""
    return [list(row) for row in zip(*face[::-1])]

def rotate_face_counter_clockwise(face):
    """Rotates a 3x3 face 90 degrees counter-clockwise."""
    return [list(row) for row in zip(*face)][::-1]

# Helper to get a column from a matrix
def get_col(matrix, c):
    return [row[c] for row in matrix]

# Helper to set a column in a matrix
def set_col(matrix, c, col_data):
    for r in range(len(matrix)):
        matrix[r][c] = col_data[r]

def apply_move(cube, move):
    """Applies a given move to the cube state."""
    state = copy.deepcopy(cube)

    if move == 'R':
        state['R'] = rotate_face_clockwise(state['R'])
        temp = get_col(state['F'], 2)
        set_col(state['F'], 2, get_col(state['D'], 2))
        set_col(state['D'], 2, list(reversed(get_col(state['B'], 0))))
        set_col(state['B'], 0, list(reversed(get_col(state['U'], 2))))
        set_col(state['U'], 2, temp)
    elif move == 'U':
        state['U'] = rotate_face_clockwise(state['U'])
        temp = state['F'][0]
        state['F'][0] = state['L'][0]
        state['L'][0] = state['B'][0]
        state['B'][0] = state['R'][0]
        state['R'][0] = temp
    elif move == 'F':
        state['F'] = rotate_face_clockwise(state['F'])
        temp = state['U'][2]
        state['U'][2] = list(reversed(get_col(state['L'], 2)))
        set_col(state['L'], 2, state['D'][0])
        state['D'][0] = list(reversed(get_col(state['R'], 0)))
        set_col(state['R'], 0, temp)
    elif move == "L'": # L prime (counter-clockwise)
        state['L'] = rotate_face_counter_clockwise(state['L'])
        temp = get_col(state['F'], 0)
        set_col(state['F'], 0, get_col(state['D'], 0))
        set_col(state['D'], 0, list(reversed(get_col(state['B'], 2))))
        set_col(state['B'], 2, list(reversed(get_col(state['U'], 0))))
        set_col(state['U'], 0, temp)
    elif move == 'D':
        state['D'] = rotate_face_clockwise(state['D'])
        temp = state['F'][2]
        state['F'][2] = state['R'][2]
        state['R'][2] = state['B'][2]
        state['B'][2] = state['L'][2]
        state['L'][2] = temp

    return state

def main():
    # Initial layouts as provided in the problem
    white_face_raw = [['R','Y','Y'],['G','W','W'],['G','O','O']]
    orange_face_raw = [['R','Y','W'],['B','O','Y'],['Y','R','O']]
    blue_face_raw = [['G','G','G'],['R','B','B'],['B','B','O']]
    yellow_face_raw = [['Y','W','R'],['B','Y','O'],['G','O','B']]
    green_face_raw = [['R','W','W'],['R','G','O'],['W','W','B']]
    red_face_raw = [['B','Y','Y'],['R','R','G'],['W','G','O']]
    
    # Correcting orientations to match standard F/U/R/B/L/D view
    # F(W), U(O), R(B) are already in standard orientation
    # B(Y) needs 90 deg counter-clockwise rotation
    # L(G) needs 90 deg clockwise rotation
    # D(R) needs 180 deg rotation
    cube = {
        'F': white_face_raw,
        'U': orange_face_raw,
        'R': blue_face_raw,
        'B': rotate_face_counter_clockwise(yellow_face_raw),
        'L': rotate_face_clockwise(green_face_raw),
        'D': rotate_face_clockwise(rotate_face_clockwise(red_face_raw))
    }
    
    # Algorithm to execute
    algorithm = ['R', 'U', 'F', "L'", 'D']
    
    # Apply each move in the sequence
    for move in algorithm:
        cube = apply_move(cube, move)
        
    # Final state of the white (front) face
    final_front_face = cube['F']
    print(final_front_face)

if __name__ == "__main__":
    main()
<<<E>>>