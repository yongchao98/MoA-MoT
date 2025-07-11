import copy

def rotate_face_cw(face):
    """Rotates a 3x3 face clockwise."""
    return [[face[2][0], face[1][0], face[0][0]],
            [face[2][1], face[1][1], face[0][1]],
            [face[2][2], face[1][2], face[0][2]]]

def rotate_face_ccw(face):
    """Rotates a 3x3 face counter-clockwise."""
    return [[face[0][2], face[1][2], face[2][2]],
            [face[0][1], face[1][1], face[2][1]],
            [face[0][0], face[1][0], face[2][0]]]

def apply_move(cube_state, move):
    """Applies a single move to the cube state."""
    c = copy.deepcopy(cube_state)
    new_c = copy.deepcopy(cube_state)

    if move == 'R':
        # Rotate Right face CW
        new_c['R'] = rotate_face_cw(c['R'])
        # Move adjacent sides: U -> F -> D -> B(rev) -> U
        for i in range(3):
            new_c['U'][i][2] = c['B'][2-i][0]
            new_c['F'][i][2] = c['U'][i][2]
            new_c['D'][i][2] = c['F'][i][2]
            new_c['B'][2-i][0] = c['D'][i][2]
            
    elif move == 'U':
        # Rotate Up face CW
        new_c['U'] = rotate_face_cw(c['U'])
        # Move adjacent sides: F -> L -> B -> R -> F
        new_c['F'][0] = c['L'][0]
        new_c['L'][0] = c['B'][0]
        new_c['B'][0] = c['R'][0]
        new_c['R'][0] = c['F'][0]
        
    elif move == 'F':
        # Rotate Front face CW
        new_c['F'] = rotate_face_cw(c['F'])
        # Move adjacent sides: U -> R -> D(rev) -> L(rev) -> U(rev)
        for i in range(3):
            # U.bottom -> R.left
            new_c['R'][i][0] = c['U'][2][i]
            # R.left -> D.top (reversed)
            new_c['D'][0][i] = c['R'][2-i][0]
            # D.top -> L.right (reversed)
            new_c['L'][i][2] = c['D'][0][2-i]
            # L.right -> U.bottom (reversed)
            new_c['U'][2][i] = c['L'][2-i][2]

    elif move == "L'":
        # Rotate Left face CCW
        new_c['L'] = rotate_face_ccw(c['L'])
        # Move adjacent sides (CCW): F -> U -> B(rev) -> D -> F
        for i in range(3):
            new_c['F'][i][0] = c['U'][i][0]
            new_c['U'][i][0] = c['B'][2-i][2]
            new_c['B'][2-i][2] = c['D'][i][0]
            new_c['D'][i][0] = c['F'][i][0]
            
    elif move == 'D':
        # Rotate Down face CW
        new_c['D'] = rotate_face_cw(c['D'])
        # Move adjacent sides: F -> R -> B -> L -> F
        new_c['F'][2] = c['R'][2]
        new_c['R'][2] = c['B'][2]
        new_c['B'][2] = c['L'][2]
        new_c['L'][2] = c['F'][2]
        
    return new_c

def solve():
    # Initial state of the jumbled cube
    # F:White, U:Orange, R:Blue, B:Yellow, L:Green, D:Red
    cube = {
        'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White
        'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange
        'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue
        'B': [['Y','W','R'],['B','Y','O'],['G','O','B']], # Yellow
        'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green
        'D': [['B','Y','Y'],['R','R','G'],['W','G','O']]  # Red
    }
    
    # The 5-step algorithm
    algorithm = ['R', 'U', 'F', "L'", 'D']
    
    # Apply each move in the sequence
    current_cube_state = cube
    for move in algorithm:
        current_cube_state = apply_move(current_cube_state, move)
        
    final_front_face = current_cube_state['F']
    
    # Print the final state of the front face as an "equation"
    print(f"[{final_front_face[0][0]},{final_front_face[0][1]},{final_front_face[0][2]}]")
    print(f"[{final_front_face[1][0]},{final_front_face[1][1]},{final_front_face[1][2]}]")
    print(f"[{final_front_face[2][0]},{final_front_face[2][1]},{final_front_face[2][2]}]")

solve()
<<<E>>>