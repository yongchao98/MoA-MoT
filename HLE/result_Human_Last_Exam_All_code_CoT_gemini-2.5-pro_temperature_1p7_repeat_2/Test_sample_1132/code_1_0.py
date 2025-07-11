import copy

def solve_rubiks_face():
    """
    Simulates a sequence of moves on a Rubik's cube and prints the final state of the front face.
    """

    # Helper functions for rotating a 3x3 face
    def rotate_face_cw(face):
        """Rotates a 3x3 matrix clockwise."""
        return [list(row) for row in zip(*face[::-1])]

    def rotate_face_ccw(face):
        """Rotates a 3x3 matrix counter-clockwise."""
        return [list(row) for row in zip(*[row for row in face])][::-1]

    # --- Move Functions (operate on the cube dictionary in-place) ---

    def move_R(cube):
        """Performs a clockwise R (Right) move."""
        cube['R'] = rotate_face_cw(cube['R'])
        temp_col = [cube['F'][i][2] for i in range(3)]
        for i in range(3): cube['F'][i][2] = cube['D'][i][2]      # D -> F
        for i in range(3): cube['D'][i][2] = cube['B'][2-i][0]    # B -> D
        for i in range(3): cube['B'][2-i][0] = cube['U'][i][2]    # U -> B
        for i in range(3): cube['U'][i][2] = temp_col[i]         # F_old -> U

    def move_U(cube):
        """Performs a clockwise U (Up) move."""
        cube['U'] = rotate_face_cw(cube['U'])
        temp_row = cube['F'][0][:]
        cube['F'][0] = cube['R'][0][:]                           # R -> F
        cube['R'][0] = cube['B'][0][:]                           # B -> R
        cube['B'][0] = cube['L'][0][:]                           # L -> B
        cube['L'][0] = temp_row                                  # F_old -> L

    def move_F(cube):
        """Performs a clockwise F (Front) move."""
        cube['F'] = rotate_face_cw(cube['F'])
        temp_row = cube['U'][2][:]
        for i in range(3): cube['U'][2][i] = cube['L'][2-i][2]    # L -> U
        for i in range(3): cube['L'][i][2] = cube['D'][0][i]      # D -> L
        for i in range(3): cube['D'][0][i] = cube['R'][2-i][0]    # R -> D
        for i in range(3): cube['R'][i][0] = temp_row[i]         # U_old -> R

    def move_L_prime(cube):
        """Performs a counter-clockwise L' (Left) move."""
        cube['L'] = rotate_face_ccw(cube['L'])
        temp_col = [cube['F'][i][0] for i in range(3)]
        for i in range(3): cube['F'][i][0] = cube['U'][i][0]      # U -> F
        for i in range(3): cube['U'][i][0] = cube['B'][2-i][2]    # B -> U
        for i in range(3): cube['B'][2-i][2] = cube['D'][i][0]    # D -> B
        for i in range(3): cube['D'][i][0] = temp_col[i]         # F_old -> D

    def move_D(cube):
        """Performs a clockwise D (Down) move."""
        cube['D'] = rotate_face_cw(cube['D'])
        temp_row = cube['F'][2][:]
        cube['F'][2] = cube['L'][2][:]                           # L -> F
        cube['L'][2] = cube['B'][2][:]                           # B -> L
        cube['B'][2] = cube['R'][2][:]                           # R -> B
        cube['R'][2] = temp_row                                  # F_old -> R

    # Initial jumbled state of the cube
    cube = {
        'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White face
        'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange face
        'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue face
        'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green face
        'B': [['Y','W','R'],['B','Y','O'],['G','O','B']], # Yellow face
        'D': [['B','Y','Y'],['R','R','G'],['W','G','O']]  # Red face
    }

    # Apply the 5-step algorithm in sequence
    move_R(cube)
    move_U(cube)
    move_F(cube)
    move_L_prime(cube)
    move_D(cube)

    # Print the final state of the white face
    final_white_face = cube['F']
    print(final_white_face)

solve_rubiks_face()