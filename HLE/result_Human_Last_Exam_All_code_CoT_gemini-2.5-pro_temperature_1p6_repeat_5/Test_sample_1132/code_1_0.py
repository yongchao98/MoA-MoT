import copy

def main():
    """
    Solves the Rubik's cube rotation problem.
    """

    # Helper function to rotate a 3x3 matrix
    def rotate_face(face, clockwise=True):
        if clockwise:
            # Transpose and reverse each row for clockwise rotation
            return [list(row) for row in zip(*face[::-1])]
        else:
            # Transpose and reverse main diagonal for counter-clockwise
            return [list(row)[::-1] for row in zip(*face)]

    # Initial layouts as described in the problem
    initial_layouts = {
        'white': {'matrix': [['R','Y','Y'],['G','W','W'],['G','O','O']], 'top': 'O', 'right': 'B'},
        'orange': {'matrix': [['R','Y','W'],['B','O','Y'],['Y','R','O']], 'top': 'Y', 'right': 'B'},
        'blue': {'matrix': [['G','G','G'],['R','B','B'],['B','B','O']], 'top': 'O', 'right': 'Y'},
        'yellow': {'matrix': [['Y','W','R'],['B','Y','O'],['G','O','B']], 'top': 'B', 'right': 'O'},
        'green': {'matrix': [['R','W','W'],['R','G','O'],['W','W','B']], 'top': 'Y', 'right': 'O'},
        'red': {'matrix': [['B','Y','Y'],['R','R','G'],['W','G','O']], 'top': 'Y', 'right': 'G'}
    }

    # Standard orientation: F(W), U(O), R(B), B(Y), L(G), D(R)
    # The views for yellow, green, and red need to be corrected to this standard.
    # Yellow (Back): Std top=O, right=G. Given top=B, right=O. This is a 90 deg CW view rotation. Correct by rotating matrix 90 deg CCW.
    # Green (Left): Std top=O, right=W. Given top=Y, right=O. This is also a 90 deg CW view rotation. Correct by rotating matrix 90 deg CCW.
    # Red (Down): Std top=W, right=B. Given top=Y, right=G. This is a 180 deg view rotation. Correct by rotating matrix 180 deg.

    faces = {
        'F': initial_layouts['white']['matrix'],  # White
        'U': initial_layouts['orange']['matrix'], # Orange
        'R': initial_layouts['blue']['matrix'],   # Blue
        'B': rotate_face(initial_layouts['yellow']['matrix'], clockwise=False), # Yellow corrected
        'L': rotate_face(initial_layouts['green']['matrix'], clockwise=False),  # Green corrected
        'D': rotate_face(rotate_face(initial_layouts['red']['matrix'])), # Red corrected (2x CW = 180)
    }

    # --- Move Definitions ---

    def move_R():
        faces['R'] = rotate_face(faces['R'])
        temp_col = [row[2] for row in faces['F']]
        for i in range(3): faces['F'][i][2] = faces['D'][i][2]
        for i in range(3): faces['D'][i][2] = faces['B'][2-i][0]
        for i in range(3): faces['B'][2-i][0] = faces['U'][i][2]
        for i in range(3): faces['U'][i][2] = temp_col[i]

    def move_U():
        faces['U'] = rotate_face(faces['U'])
        temp_row = faces['F'][0]
        faces['F'][0] = faces['R'][0]
        faces['R'][0] = faces['B'][0]
        faces['B'][0] = faces['L'][0]
        faces['L'][0] = temp_row

    def move_F():
        faces['F'] = rotate_face(faces['F'])
        temp_row = faces['U'][2]
        for i in range(3): faces['U'][2][i] = faces['L'][2-i][2]
        for i in range(3): faces['L'][i][2] = faces['D'][0][i]
        for i in range(3): faces['D'][0][i] = faces['R'][2-i][0]
        for i in range(3): faces['R'][i][0] = temp_row[i]

    def move_L_prime():
        faces['L'] = rotate_face(faces['L'], clockwise=False)
        temp_col = [row[0] for row in faces['F']]
        for i in range(3): faces['F'][i][0] = faces['U'][i][0]
        for i in range(3): faces['U'][i][0] = faces['B'][2-i][2]
        for i in range(3): faces['B'][2-i][2] = faces['D'][i][0]
        for i in range(3): faces['D'][i][0] = temp_col[i]

    def move_D():
        faces['D'] = rotate_face(faces['D'])
        temp_row = faces['F'][2]
        faces['F'][2] = faces['L'][2]
        faces['L'][2] = faces['B'][2]
        faces['B'][2] = faces['R'][2]
        faces['R'][2] = temp_row

    # Execute the 5-step algorithm
    move_R()
    move_U()
    move_F()
    move_L_prime()
    move_D()

    # Print the final state of the white (Front) face
    final_face = faces['F']
    print(f"Final White Face State:")
    print(f"[[{final_face[0][0]},{final_face[0][1]},{final_face[0][2]}],[{final_face[1][0]},{final_face[1][1]},{final_face[1][2]}],[{final_face[2][0]},{final_face[2][1]},{final_face[2][2]}]]")

if __name__ == '__main__':
    main()