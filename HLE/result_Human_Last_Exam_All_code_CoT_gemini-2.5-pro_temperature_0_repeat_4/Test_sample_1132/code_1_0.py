import copy

def solve_rubiks_face():
    """
    Simulates a 5-step algorithm on a Rubik's cube and prints the final state of the white face.
    """

    # Step 1 & 2: Define the initial state of the cube based on the problem description.
    # F=White, U=Orange, R=Blue, L=Green, B=Yellow, D=Red
    state = {
        'F': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']], # White
        'U': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']], # Orange
        'R': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']], # Blue
        'B': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']], # Yellow
        'L': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']], # Green
        'D': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]  # Red
    }

    # Helper function to rotate a face's 3x3 matrix
    def rotate_face_cw(face):
        return [[face[2][0], face[1][0], face[0][0]],
                [face[2][1], face[1][1], face[0][1]],
                [face[2][2], face[1][2], face[0][2]]]

    def rotate_face_ccw(face):
        return [[face[0][2], face[1][2], face[2][2]],
                [face[0][1], face[1][1], face[2][1]],
                [face[0][0], face[1][0], face[2][0]]]

    # Step 3: Implement the move functions
    def move(current_state, move_type):
        s = copy.deepcopy(current_state)
        
        if move_type == 'R':
            s['R'] = rotate_face_cw(s['R'])
            # Cycle: F -> U -> B(rev) -> D -> F
            # Store values before overwriting
            temp_f_col = [s['F'][i][2] for i in range(3)]
            temp_u_col = [s['U'][i][2] for i in range(3)]
            temp_b_col = [s['B'][i][0] for i in range(3)]
            temp_d_col = [s['D'][i][2] for i in range(3)]
            # Apply rotations
            for i in range(3):
                s['U'][i][2] = temp_f_col[i]
                s['B'][i][0] = temp_u_col[2-i]
                s['D'][i][2] = temp_b_col[2-i]
                s['F'][i][2] = temp_d_col[i]

        elif move_type == 'U':
            s['U'] = rotate_face_cw(s['U'])
            # Cycle: F -> R -> B -> L -> F
            temp_f_row = s['F'][0]
            s['F'][0] = s['R'][0]
            s['R'][0] = s['B'][0]
            s['B'][0] = s['L'][0]
            s['L'][0] = temp_f_row

        elif move_type == 'F':
            s['F'] = rotate_face_cw(s['F'])
            # Cycle: U -> R -> D -> L -> U
            temp_u_row = s['U'][2]
            # U bottom gets L right reversed
            for i in range(3): s['U'][2][i] = s['L'][2-i][2]
            # L right gets D top
            for i in range(3): s['L'][i][2] = s['D'][0][i]
            # D top gets R left reversed
            for i in range(3): s['D'][0][i] = s['R'][2-i][0]
            # R left gets original U bottom
            for i in range(3): s['R'][i][0] = temp_u_row[i]

        elif move_type == "L'":
            s['L'] = rotate_face_ccw(s['L'])
            # CCW Cycle: F -> D -> B(rev) -> U -> F
            temp_f_col = [s['F'][i][0] for i in range(3)]
            # F left gets D left
            for i in range(3): s['F'][i][0] = s['D'][i][0]
            # D left gets B right reversed
            for i in range(3): s['D'][i][0] = s['B'][2-i][2]
            # B right gets U left reversed
            for i in range(3): s['B'][i][2] = s['U'][2-i][0]
            # U left gets original F left
            for i in range(3): s['U'][i][0] = temp_f_col[i]

        elif move_type == 'D':
            s['D'] = rotate_face_cw(s['D'])
            # CW Cycle: F -> L -> B -> R -> F
            temp_f_row = s['F'][2]
            s['F'][2] = s['L'][2]
            s['L'][2] = s['B'][2]
            s['B'][2] = s['R'][2]
            s['R'][2] = temp_f_row
            
        return s

    # Step 4: Execute the algorithm
    algorithm = ['R', 'U', 'F', "L'", 'D']
    final_state = state
    for m in algorithm:
        final_state = move(final_state, m)

    # Step 5: Output the final state of the front face
    final_front_face = final_state['F']
    print(final_front_face)

solve_rubiks_face()