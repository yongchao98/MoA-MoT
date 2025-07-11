import copy

def rotate_cw(m):
    """Rotates a 3x3 matrix 90 degrees clockwise."""
    return [list(row) for row in zip(*m[::-1])]

def rotate_ccw(m):
    """Rotates a 3x3 matrix 90 degrees counter-clockwise."""
    return [list(row) for row in zip(*m)][::-1]

def apply_move(state, move_char):
    """Applies a single move to the cube state."""
    s = copy.deepcopy(state)
    face = move_char[0]
    is_prime = len(move_char) > 1

    if is_prime:
        s[face] = rotate_ccw(s[face])
    else:
        s[face] = rotate_cw(s[face])
    
    # Adjacent strip movements
    if face == 'U':
        b_r, r_r, f_r, l_r = state['B'][0], state['R'][0], state['F'][0], state['L'][0]
        if not is_prime: s['R'][0], s['B'][0], s['L'][0], s['F'][0] = f_r, r_r, b_r, l_r
        else: s['L'][0], s['B'][0], s['R'][0], s['F'][0] = f_r, l_r, b_r, r_r
    
    elif face == 'D':
        f_r, l_r, b_r, r_r = state['F'][2], state['L'][2], state['B'][2], state['R'][2]
        if not is_prime: s['L'][2], s['B'][2], s['R'][2], s['F'][2] = f_r, l_r, b_r, r_r
        else: s['R'][2], s['B'][2], s['L'][2], s['F'][2] = f_r, r_r, b_r, l_r
        
    elif face == 'F':
        u_r, r_c, d_r, l_c = state['U'][2], [r[0] for r in state['R']], state['D'][0], [r[2] for r in state['L']]
        if not is_prime:
            for i in range(3): s['R'][i][0] = u_r[i]
            for i in range(3): s['D'][0][i] = r_c[2-i]
            for i in range(3): s['L'][i][2] = d_r[2-i]
            for i in range(3): s['U'][2][i] = l_c[i]
        else:
            for i in range(3): s['L'][i][2] = u_r[i]
            for i in range(3): s['U'][2][i] = l_c[2-i]
            for i in range(3): s['R'][i][0] = d_r[i]
            for i in range(3): s['D'][0][i] = r_c[2-i]
            
    elif face == 'R':
        u_c, b_c, d_c, f_c = [r[2] for r in state['U']], [r[0] for r in state['B']], [r[2] for r in state['D']], [r[2] for r in state['F']]
        if not is_prime:
            for i in range(3): s['B'][i][0] = u_c[2-i]
            for i in range(3): s['D'][i][2] = b_c[2-i]
            for i in range(3): s['F'][i][2] = d_c[i]
            for i in range(3): s['U'][i][2] = f_c[i]
        else:
            for i in range(3): s['F'][i][2] = u_c[i]
            for i in range(3): s['U'][i][2] = d_c[i]
            for i in range(3): s['D'][i][2] = f_c[2-i]
            for i in range(3): s['B'][i][0] = b_c[2-i]

    elif face == 'L':
        u_c, f_c, d_c, b_c = [r[0] for r in state['U']], [r[0] for r in state['F']], [r[0] for r in state['D']], [r[2] for r in state['B']]
        if not is_prime:
            for i in range(3): s['F'][i][0] = u_c[i]
            for i in range(3): s['D'][i][0] = f_c[i]
            for i in range(3): s['B'][i][2] = d_c[2-i]
            for i in range(3): s['U'][i][0] = b_c[2-i]
        else:
            for i in range(3): s['U'][i][0] = f_c[i]
            for i in range(3): s['F'][i][0] = d_c[i]
            for i in range(3): s['D'][i][0] = b_c[2-i]
            for i in range(3): s['B'][i][2] = u_c[2-i]
            
    return s

def solve():
    # Initial state of the jumbled cube
    # F=Front(White), U=Up(Orange), R=Right(Blue)
    # B=Back(Yellow), L=Left(Green), D=Down(Red)
    cube = {
        'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White face
        'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange face
        'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue face
        'B': [['Y','W','R'],['B','Y','O'],['G','O','B']], # Yellow face
        'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green face
        'D': [['B','Y','Y'],['R','R','G'],['W','G','O']]  # Red face
    }

    algorithm = ["R", "U", "F", "L'", "D"]
    
    current_cube_state = cube
    for move in algorithm:
        current_cube_state = apply_move(current_cube_state, move)
        
    final_white_face = current_cube_state['F']

    print("The final state of the white face is:")
    print(f"[[{final_white_face[0][0]},{final_white_face[0][1]},{final_white_face[0][2]}],[{final_white_face[1][0]},{final_white_face[1][1]},{final_white_face[1][2]}],[{final_white_face[2][0]},{final_white_face[2][1]},{final_white_face[2][2]}]]")

solve()