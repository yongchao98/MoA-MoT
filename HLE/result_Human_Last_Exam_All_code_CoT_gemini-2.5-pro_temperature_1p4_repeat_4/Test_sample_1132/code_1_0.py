import copy

def solve_rubiks_puzzle():
    """
    This script simulates a 5-step algorithm on a jumbled Rubik's cube
    and prints the final state of the white (Front) face.
    """

    # Helper functions for rotating the 3x3 matrix of a single face
    def rotate_face_cw(face):
        """Rotates a 3x3 face clockwise"""
        return [list(row) for row in zip(*face[::-1])]

    def rotate_face_ccw(face):
        """Rotates a 3x3 face counter-clockwise"""
        # This is equivalent to three clockwise rotations
        new_face = face
        for _ in range(3):
            new_face = rotate_face_cw(new_face)
        return new_face

    # --- Move Functions ---
    # Each function takes the cube state `c` and modifies it in place.

    def move_R(c):
        """Performs a clockwise R move."""
        c['R'] = rotate_face_cw(c['R'])
        
        # Store the adjacent strips before modifying them
        f_col2 = [c['F'][i][2] for i in range(3)]
        u_col2 = [c['U'][i][2] for i in range(3)]
        b_col0 = [c['B'][i][0] for i in range(3)]
        d_col2 = [c['D'][i][2] for i in range(3)]

        # Cycle the strips: F -> U -> B(rev) -> D(rev) -> F
        for i in range(3): c['F'][i][2] = d_col2[i]
        for i in range(3): c['U'][i][2] = f_col2[i]
        for i in range(3): c['B'][i][0] = u_col2[2-i]
        for i in range(3): c['D'][i][2] = b_col0[2-i]

    def move_U(c):
        """Performs a clockwise U move."""
        c['U'] = rotate_face_cw(c['U'])
        
        f_row0, r_row0, b_row0, l_row0 = c['F'][0][:], c['R'][0][:], c['B'][0][:], c['L'][0][:]
        
        # Cycle the strips: F -> R -> B -> L -> F
        c['F'][0], c['R'][0], c['B'][0], c['L'][0] = r_row0, b_row0, l_row0, f_row0
    
    def move_F(c):
        """Performs a clockwise F move."""
        c['F'] = rotate_face_cw(c['F'])

        u_row2 = c['U'][2][:]
        r_col0 = [c['R'][i][0] for i in range(3)]
        d_row0 = c['D'][0][:]
        l_col2 = [c['L'][i][2] for i in range(3)]

        # Cycle: U -> R -> D(rev) -> L(rev) -> U(rev)
        for i in range(3): c['R'][i][0] = u_row2[i]          # U bottom row -> R left col
        c['D'][0] = r_col0[::-1]                              # R left col -> D top row (rev)
        for i in range(3): c['L'][i][2] = d_row0[::-1][i]    # D top row -> L right col (rev)
        c['U'][2] = l_col2[::-1]                              # L right col -> U bottom row (rev)

    def move_L_prime(c):
        """Performs a counter-clockwise L' move."""
        c['L'] = rotate_face_ccw(c['L'])

        f_col0 = [c['F'][i][0] for i in range(3)]
        u_col0 = [c['U'][i][0] for i in range(3)]
        b_col2 = [c['B'][i][2] for i in range(3)]
        d_col0 = [c['D'][i][0] for i in range(3)]

        # Cycle (L' is opposite of L): F -> D -> B(rev) -> U -> F
        for i in range(3): c['F'][i][0] = d_col0[i]
        for i in range(3): c['D'][i][0] = b_col2[2-i]
        for i in range(3): c['B'][i][2] = u_col0[2-i]
        for i in range(3): c['U'][i][0] = f_col0[i]

    def move_D(c):
        """Performs a clockwise D move."""
        c['D'] = rotate_face_cw(c['D'])

        f_row2, r_row2, b_row2, l_row2 = c['F'][2][:], c['R'][2][:], c['B'][2][:], c['L'][2][:]

        # Cycle: F -> L -> B -> R -> F
        c['F'][2], c['L'][2], c['B'][2], c['R'][2] = l_row2, b_row2, r_row2, f_row2

    # --- Main Execution ---

    # Initial Jumbled Cube State
    initial_cube = {
        'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White face
        'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange face (Up)
        'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue face (Right)
        'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green face (Left)
        'D': [['B','Y','Y'],['R','R','G'],['W','G','O']], # Red face (Down)
        'B': [['Y','W','R'],['B','Y','O'],['G','O','B']]  # Yellow face (Back)
    }

    # Create a mutable copy to perform the simulation on
    cube_state = copy.deepcopy(initial_cube)
    
    # The 5-step algorithm
    algorithm = [move_R, move_U, move_F, move_L_prime, move_D]

    # Apply each move in the algorithm to the cube state
    for move in algorithm:
        move(cube_state)
    
    # Print the final state of the white face
    final_white_face = cube_state['F']
    print(final_white_face)

solve_rubiks_puzzle()