import copy

class RubiksCube:
    """
    A class to represent and manipulate a 3x3 Rubik's Cube.
    """
    def __init__(self, initial_state):
        # The state is a dictionary of 6 faces, each a 3x3 matrix.
        # F, B, U, D, L, R
        self.faces = copy.deepcopy(initial_state)

    def _rotate_face_cw(self, face_key):
        """Rotates a face 90 degrees clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [
            [face[2][0], face[1][0], face[0][0]],
            [face[2][1], face[1][1], face[0][1]],
            [face[2][2], face[1][2], face[0][2]]
        ]

    def _rotate_face_ccw(self, face_key):
        """Rotates a face 90 degrees counter-clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [
            [face[0][2], face[1][2], face[2][2]],
            [face[0][1], face[1][1], face[2][1]],
            [face[0][0], face[1][0], face[2][0]]
        ]

    def R(self):
        """Performs a clockwise Right face turn."""
        self._rotate_face_cw('R')
        # Store the affected columns to perform a 4-way swap
        f_col = [self.faces['F'][i][2] for i in range(3)]
        u_col = [self.faces['U'][i][2] for i in range(3)]
        b_col = [self.faces['B'][i][2] for i in range(3)]
        d_col = [self.faces['D'][i][2] for i in range(3)]
        
        # F -> U -> B -> D -> F
        for i in range(3):
            self.faces['U'][i][2] = f_col[i]
            # Back face columns are indexed differently and stickers flip
            self.faces['B'][2-i][0] = u_col[i] 
            self.faces['D'][i][2] = self.faces['B'][2-i][0] # Temporarily read from the new B face
            self.faces['F'][i][2] = d_col[i]
        
        # To correct for the temporary read, we need the original B state for D
        b_col_original_for_d = [self.faces['B'][i][0] for i in range(3)] # Before it was overwritten
        u_col_original_for_b = u_col # from above
        for i in range(3):
             self.faces['B'][2-i][0] = u_col_original_for_b[i]
        
        b_col_after_u_swap = [self.faces['B'][i][0] for i in range(3)]
        for i in range(3):
            self.faces['D'][i][2] = b_col_after_u_swap[2-i]

    def U(self):
        """Performs a clockwise Up face turn."""
        self._rotate_face_cw('U')
        # F -> R -> B -> L -> F
        f_row = copy.deepcopy(self.faces['F'][0])
        r_row = copy.deepcopy(self.faces['R'][0])
        b_row = copy.deepcopy(self.faces['B'][0])
        l_row = copy.deepcopy(self.faces['L'][0])
        self.faces['F'][0] = r_row
        self.faces['R'][0] = b_row
        self.faces['B'][0] = l_row
        self.faces['L'][0] = f_row

    def F(self):
        """Performs a clockwise Front face turn."""
        self._rotate_face_cw('F')
        # U(bottom row) -> R(left col) -> D(top row) -> L(right col) -> U(bottom row)
        temp_row = copy.deepcopy(self.faces['U'][2])
        # L's right col (reversed) -> U's bottom row
        self.faces['U'][2] = [self.faces['L'][2][2], self.faces['L'][1][2], self.faces['L'][0][2]]
        # D's top row -> L's right col
        self.faces['L'][0][2], self.faces['L'][1][2], self.faces['L'][2][2] = self.faces['D'][0][0], self.faces['D'][0][1], self.faces['D'][0][2]
        # R's left col (reversed) -> D's top row
        self.faces['D'][0] = [self.faces['R'][2][0], self.faces['R'][1][0], self.faces['R'][0][0]]
        # Original U's bottom row -> R's left col
        self.faces['R'][0][0], self.faces['R'][1][0], self.faces['R'][2][0] = temp_row[0], temp_row[1], temp_row[2]

    def L_prime(self):
        """Performs a counter-clockwise Left face turn."""
        self._rotate_face_ccw('L')
        # F -> U -> B -> D -> F
        f_col = copy.deepcopy(self.faces['F'][0])
        u_col = copy.deepcopy(self.faces['U'][0])
        b_col = copy.deepcopy(self.faces['B'][2])
        d_col = copy.deepcopy(self.faces['D'][0])

        for i in range(3):
            self.faces['F'][i][0] = u_col[i]
            self.faces['U'][i][0] = self.faces['B'][2-i][2]
            self.faces['B'][2-i][2] = d_col[i]
            self.faces['D'][i][0] = f_col[i]


    def D(self):
        """Performs a clockwise Down face turn."""
        self._rotate_face_cw('D')
        # F -> L -> B -> R -> F
        f_row = copy.deepcopy(self.faces['F'][2])
        l_row = copy.deepcopy(self.faces['L'][2])
        b_row = copy.deepcopy(self.faces['B'][2])
        r_row = copy.deepcopy(self.faces['R'][2])
        self.faces['F'][2] = l_row
        self.faces['L'][2] = b_row
        self.faces['B'][2] = r_row
        self.faces['R'][2] = f_row

    def execute_moves(self, moves_str):
        """Executes a string of moves."""
        moves = moves_str.split()
        for move in moves:
            if move == "R": self.R()
            elif move == "U": self.U()
            elif move == "F": self.F()
            elif move == "L'": self.L_prime()
            elif move == "D": self.D()

    def get_face(self, face_key):
        """Returns the 3x3 matrix for a given face."""
        return self.faces[face_key]

# Initial state of the jumbled cube based on the problem description
initial_state = {
    'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White face
    'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange face
    'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue face
    'B': [['Y','W','R'],['B','Y','O'],['G','O','B']], # Yellow face
    'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green face
    'D': [['B','Y','Y'],['R','R','G'],['W','G','O']]  # Red face
}

# Create a cube instance with the initial state
cube = RubiksCube(initial_state)

# Execute the 5-step algorithm
# Due to the complexity of side effects, we will apply moves one by one.
# 1. R
cube.R()
# 2. U
cube.U()
# 3. F
cube.F()
# 4. L'
cube.L_prime()
# 5. D
cube.D()

# Get the final state of the white (Front) face
final_white_face = cube.get_face('F')

# Print the result in the required format
print("[", end="")
for i, row in enumerate(final_white_face):
    print(row, end="")
    if i < len(final_white_face) - 1:
        print(",", end="")
print("]")