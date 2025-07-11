import copy

class RubiksCube:
    def __init__(self, initial_state):
        # Using deepcopy to ensure the initial state is not modified elsewhere
        self.faces = copy.deepcopy(initial_state)

    def _rotate_face_cw(self, face):
        """Rotates a 3x3 face 90 degrees clockwise."""
        return [[face[2 - j][i] for j in range(3)] for i in range(3)]

    def _rotate_face_ccw(self, face):
        """Rotates a 3x3 face 90 degrees counter-clockwise."""
        return [[face[j][2 - i] for j in range(3)] for i in range(3)]

    def move_R(self):
        """Performs a 90-degree clockwise rotation of the Right face."""
        old_faces = copy.deepcopy(self.faces)
        self.faces['R'] = self._rotate_face_cw(old_faces['R'])
        
        # Cycle of adjacent columns: F -> U -> B -> D -> F
        for i in range(3):
            self.faces['U'][i][2] = old_faces['F'][i][2]
            self.faces['B'][i][2] = old_faces['U'][i][2]
            self.faces['D'][i][2] = old_faces['B'][i][2]
            self.faces['F'][i][2] = old_faces['D'][i][2]

    def move_U(self):
        """Performs a 90-degree clockwise rotation of the Up face."""
        old_faces = copy.deepcopy(self.faces)
        self.faces['U'] = self._rotate_face_cw(old_faces['U'])
        
        # Cycle of adjacent rows: F -> R -> B -> L -> F
        self.faces['R'][0] = old_faces['F'][0]
        self.faces['B'][0] = old_faces['R'][0]
        self.faces['L'][0] = old_faces['B'][0]
        self.faces['F'][0] = old_faces['L'][0]

    def move_F(self):
        """Performs a 90-degree clockwise rotation of the Front face."""
        old_faces = copy.deepcopy(self.faces)
        self.faces['F'] = self._rotate_face_cw(old_faces['F'])
        
        # Store adjacent strips
        temp_U_bottom = old_faces['U'][2]
        temp_R_left = [old_faces['R'][i][0] for i in range(3)]
        temp_D_top = old_faces['D'][0]
        temp_L_right = [old_faces['L'][i][2] for i in range(3)]

        # U -> R (no reversal)
        for i in range(3): self.faces['R'][i][0] = temp_U_bottom[i]
        # R -> D (reversed)
        for i in range(3): self.faces['D'][0][i] = temp_R_left[2 - i]
        # D -> L (reversed)
        for i in range(3): self.faces['L'][i][2] = temp_D_top[2 - i]
        # L -> U (reversed)
        for i in range(3): self.faces['U'][2][i] = temp_L_right[2 - i]

    def move_L(self):
        """Performs a 90-degree clockwise rotation of the Left face."""
        old_faces = copy.deepcopy(self.faces)
        self.faces['L'] = self._rotate_face_cw(old_faces['L'])

        # Cycle: F -> D -> B(rev) -> U(rev) -> F
        temp_F_left = [old_faces['F'][i][0] for i in range(3)]
        temp_D_left = [old_faces['D'][i][0] for i in range(3)]
        temp_B_right = [old_faces['B'][i][2] for i in range(3)]
        temp_U_left = [old_faces['U'][i][0] for i in range(3)]

        for i in range(3): self.faces['D'][i][0] = temp_F_left[i]
        for i in range(3): self.faces['B'][i][2] = temp_D_left[2 - i]
        for i in range(3): self.faces['U'][i][0] = temp_B_right[2 - i]
        for i in range(3): self.faces['F'][i][0] = temp_U_left[i]

    def move_L_prime(self):
        """Performs a 90-degree counter-clockwise L-face turn."""
        self.move_L()
        self.move_L()
        self.move_L()

    def move_D(self):
        """Performs a 90-degree clockwise rotation of the Down face."""
        old_faces = copy.deepcopy(self.faces)
        self.faces['D'] = self._rotate_face_cw(old_faces['D'])

        # Cycle of adjacent rows: F -> L -> B -> R -> F
        self.faces['L'][2] = old_faces['F'][2]
        self.faces['B'][2] = old_faces['L'][2]
        self.faces['R'][2] = old_faces['B'][2]
        self.faces['F'][2] = old_faces['R'][2]

# Initial state of the jumbled cube
initial_state = {
    'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White face
    'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange face (Top)
    'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue face (Right)
    'B': [['Y','W','R'],['B','Y','O'],['G','O','B']], # Yellow face (Back)
    'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green face (Left)
    'D': [['B','Y','Y'],['R','R','G'],['W','G','O']]  # Red face (Bottom)
}

# Create a cube instance
cube = RubiksCube(initial_state)

# Apply the 5-step algorithm
cube.move_R()
cube.move_U()
cube.move_F()
cube.move_L_prime()
cube.move_D()

# Print the final state of the white (Front) face
final_white_face = cube.faces['F']
for row in final_white_face:
    print(row)