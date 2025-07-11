import copy

def rotate_matrix_cw(m):
    """Rotates a 3x3 matrix 90 degrees clockwise."""
    return [list(row) for row in zip(*m[::-1])]

def rotate_matrix_ccw(m):
    """Rotates a 3x3 matrix 90 degrees counter-clockwise."""
    return [list(row) for row in zip(*m)][::-1]

def solve():
    """
    Solves the Rubik's cube problem by simulating the moves.
    """
    # Step 1: Define initial face layouts as given in the problem.
    # The 'top' and 'right' descriptions for each face mean we need to reorient
    # some of them to fit a standard F/U/R view (Front=White, Up=Orange, Right=Blue).
    
    # White face (Front): U=O, R=B -> Standard view, no change.
    face_F = [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']]
    
    # Orange face (Up): U=Y, R=B. This means if Orange is Front, Yellow is Up.
    # Tilting back (x-axis rotation) makes Orange Up and White Front. This matches the standard view. No change.
    face_U = [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']]

    # Blue face (Right): U=O, R=Y. If Blue is Front, Orange is Up.
    # Rotating (y-axis rotation) makes White Front and Orange Up. Blue moves to Right. This matches. No change.
    face_R = [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']]

    # Yellow face (Back): U=B, R=O. If Yellow is Front, Blue is Up.
    # To get to our standard view (W front, O up), a cube rotation is needed which
    # results in rotating the yellow face matrix 90 degrees counter-clockwise.
    face_B_given = [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']]
    face_B = rotate_matrix_ccw(face_B_given)
    
    # Green face (Left): U=Y, R=O. If Green is Front, Yellow is Up.
    # A cube rotation is needed which results in rotating the green face matrix 90 degrees clockwise.
    face_G_given = [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']]
    face_L = rotate_matrix_cw(face_G_given)

    # Red face (Down): U=Y, R=G. If Red is Front, Yellow is Up.
    # A cube rotation is needed which results in rotating the red face matrix 90 degrees clockwise.
    face_R_given = [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]
    face_D = rotate_matrix_cw(face_R_given)

    faces = {'F': face_F, 'U': face_U, 'R': face_R, 'B': face_B, 'L': face_L, 'D': face_D}

    # Step 2: Define the moves based on standard cube notation
    
    moves = ["R", "U", "F", "L'", "D"]

    for move in moves:
        if move == "R":
            faces['R'] = rotate_matrix_cw(faces['R'])
            temp = [faces['U'][i][2] for i in range(3)]
            for i in range(3): faces['U'][i][2] = faces['F'][i][2]
            for i in range(3): faces['F'][i][2] = faces['D'][i][2]
            for i in range(3): faces['D'][i][2] = faces['B'][i][2]
            for i in range(3): faces['B'][i][2] = temp
        
        elif move == "U":
            faces['U'] = rotate_matrix_cw(faces['U'])
            temp = faces['F'][0]
            faces['F'][0] = faces['R'][0]
            faces['R'][0] = faces['B'][0]
            faces['B'][0] = faces['L'][0]
            faces['L'][0] = temp

        elif move == "F":
            faces['F'] = rotate_matrix_cw(faces['F'])
            temp = faces['U'][2]
            faces['U'][2] = [faces['L'][2-i][2] for i in range(3)]
            for i in range(3): faces['L'][i][2] = faces['D'][0][i]
            faces['D'][0] = [faces['R'][2-i][0] for i in range(3)]
            for i in range(3): faces['R'][i][0] = temp[i]

        elif move == "L'":
            faces['L'] = rotate_matrix_ccw(faces['L'])
            temp = [faces['U'][i][0] for i in range(3)]
            for i in range(3): faces['U'][i][0] = faces['B'][i][0]
            for i in range(3): faces['B'][i][0] = faces['D'][i][0]
            for i in range(3): faces['D'][i][0] = faces['F'][i][0]
            for i in range(3): faces['F'][i][0] = temp[i]
            
        elif move == "D":
            faces['D'] = rotate_matrix_cw(faces['D'])
            temp = faces['F'][2]
            faces['F'][2] = faces['L'][2]
            faces['L'][2] = faces['B'][2]
            faces['B'][2] = faces['R'][2]
            faces['R'][2] = temp

    # Step 3: Print the final state of the front (white) face
    final_face = faces['F']
    for row in final_face:
        print(row)

solve()
<<<E>>>