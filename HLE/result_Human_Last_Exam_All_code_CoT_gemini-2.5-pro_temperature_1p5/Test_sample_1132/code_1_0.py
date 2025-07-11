import copy

class RubiksCube:
    """
    A class to represent and manipulate a 3x3 Rubik's Cube.
    The cube state is stored as a dictionary of 6 faces,
    each represented by a 3x3 list of lists (matrix).
    """
    def __init__(self, initial_state):
        self.faces = copy.deepcopy(initial_state)

    def get_face(self, face_key):
        return self.faces[face_key]

    def get_col(self, face_key, col_idx):
        return [self.faces[face_key][i][col_idx] for i in range(3)]

    def set_col(self, face_key, col_idx, new_col):
        for i in range(3):
            self.faces[face_key][i][col_idx] = new_col[i]

    def rotate_face_cw(self, face_key):
        face = self.faces[face_key]
        self.faces[face_key] = [[face[2][0], face[1][0], face[0][0]],
                                [face[2][1], face[1][1], face[0][1]],
                                [face[2][2], face[1][2], face[0][2]]]

    def rotate_face_ccw(self, face_key):
        face = self.faces[face_key]
        self.faces[face_key] = [[face[0][2], face[1][2], face[2][2]],
                                [face[0][1], face[1][1], face[2][1]],
                                [face[0][0], face[1][0], face[2][0]]]

    # Standard Singmaster Moves
    def R(self):
        self.rotate_face_cw('R')
        temp_f = self.get_col('F', 2)
        temp_u = self.get_col('U', 2)
        temp_b = self.get_col('B', 0)
        temp_d = self.get_col('D', 2)
        # Cycle: F -> U -> B(rev) -> D(rev) -> F
        self.set_col('U', 2, temp_f)
        self.set_col('B', 0, [temp_u[2], temp_u[1], temp_u[0]])
        self.set_col('D', 2, [temp_b[2], temp_b[1], temp_b[0]])
        self.set_col('F', 2, temp_d)

    def U(self):
        self.rotate_face_cw('U')
        # Cycle: F -> R -> B -> L -> F
        temp = self.faces['F'][0][:]
        self.faces['F'][0] = self.faces['R'][0]
        self.faces['R'][0] = self.faces['B'][0]
        self.faces['B'][0] = self.faces['L'][0]
        self.faces['L'][0] = temp

    def F(self):
        self.rotate_face_cw('F')
        # Cycle: U(rev) -> L -> D -> R -> U
        temp_u = self.faces['U'][2][:]
        temp_r = self.get_col('R', 0)
        temp_d = self.faces['D'][0][:]
        temp_l = self.get_col('L', 2)
        self.set_col('R', 0, temp_u)
        self.faces['D'][0] = [temp_r[2], temp_r[1], temp_r[0]]
        self.set_col('L', 2, temp_d)
        self.faces['U'][2] = [temp_l[2], temp_l[1], temp_l[0]]

    def L_prime(self):
        self.rotate_face_ccw('L')
        # Cycle for L is F -> U -> B(rev) -> D, so L' is F -> D -> B(rev) -> U
        temp_f = self.get_col('F', 0)
        temp_d = self.get_col('D', 0)
        temp_b = self.get_col('B', 2)
        temp_u = self.get_col('U', 0)
        self.set_col('D', 0, temp_f)
        self.set_col('B', 2, [temp_d[2], temp_d[1], temp_d[0]])
        self.set_col('U', 0, [temp_b[2], temp_b[1], temp_b[0]])
        self.set_col('F', 0, temp_u)

    def D(self):
        self.rotate_face_cw('D')
        # Cycle: F -> L -> B -> R -> F
        temp = self.faces['F'][2][:]
        self.faces['F'][2] = self.faces['L'][2]
        self.faces['L'][2] = self.faces['B'][2]
        self.faces['B'][2] = self.faces['R'][2]
        self.faces['R'][2] = temp


# --- Main Execution ---

# 1. Initialize cube state based on problem description and orientation adjustments
initial_state = {
    'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White (ok)
    'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange (ok)
    'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue (ok)
    # B(Y) rotate CCW: from top=B, right=O to top=O, right=G
    'B': [['R','O','B'],['W','Y','O'],['Y','B','G']],
    # L(G) rotate CW: from top=Y, right=O to top=O, right=W
    'L': [['W','R','R'],['W','G','W'],['B','O','W']],
    # D(R) rotate 180: from top=Y, right=G to top=W, right=B
    'D': [['O','G','W'],['G','R','R'],['Y','Y','B']]
}

# 2. Create cube object
cube = RubiksCube(initial_state)

# 3. Apply the 5-step algorithm
cube.R()
cube.U()
cube.F()
cube.L_prime()
cube.D()

# 4. Get the final state of the white (Front) face
final_white_face = cube.get_face('F')

# 5. Print the result
for row in final_white_face:
    print(row)
